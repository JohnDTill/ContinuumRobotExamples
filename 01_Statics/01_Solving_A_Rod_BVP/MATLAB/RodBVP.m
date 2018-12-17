function RodBVP
    %Parameters
    E = 200e9;
    G = 80e9;
    r = 0.001;
    rho = 8000;
    g = [9.81; 0; 0];
    L = 0.5;
    A = pi*r^2;
    I = pi*r^4/4;
    J = 2*I;
    Kse = diag([G*A, G*A, E*A]);
    Kbt = diag([E*I, E*I, G*J]);
		
    %Boundary Conditions
    p0 = [0;0;0];
    R0 = eye(3);
    pL = [0; -0.1*L; 0.8*L];
    RL = eye(3);

    %Main Simulation
    init_guess = zeros(6,1);
    global y; %Forward declaration for future scoping rule changes
    fsolve(@RodShootingMethod, init_guess); %Use convex optimization to solve ICs
    
    %Visualization
    plot3(y(:,1),y(:,2),y(:,3));  title('Rod BVP Solution');  axis([-L/2 L/2 -L/2 L/2 0 L])
    grid on;  daspect([1 1 1]);  xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
    
    %Subfunctions
    function residual = RodShootingMethod(guess) %Optimization objective function
        n0 = guess(1:3);                  %Update guessed initial conditions
        m0 = guess(4:6);
        y0 = [p0; reshape(R0,9,1); n0; m0];
        
        [s,y] = ode45(@RodODE,[0 L],y0);  %Numerically solve the resulting IVP
        
        pL_shot = y(end,1:3)';            %Calculate distal constraint violation
        RL_shot = reshape(y(end,4:12),3,3);
        position_error = pL_shot - pL;
        rotation_error = inv_hat( RL_shot'*RL - RL_shot*RL' );
        residual = [position_error; rotation_error];
    end

    function ys = RodODE(s,y)
        R = reshape(y(4:12),3,3);
        n = y(13:15);
        m = y(16:18);

        v = Kse^-1*R.'*n + [0;0;1];
        u = Kbt^-1*R.'*m;

        ps = R*v;
        Rs = R*hat(u);
        ns = -rho*A*g;
        ms = -hat(ps)*n;

        ys = [ps; reshape(Rs,9,1); ns; ms];
    end

    function skew_symmetric_matrix = hat(y)
        skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                                  y(3)   0   -y(1) ;
                                 -y(2)  y(1)   0  ];
    end

    function R3 = inv_hat(skew)
        R3 = [skew(3,2); skew(1,3); skew(2,1)];
    end
end