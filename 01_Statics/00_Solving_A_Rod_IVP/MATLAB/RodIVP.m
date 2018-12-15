function RodIVP
    %Independent Parameters
    E = 200e9;        %Young's modulus
    G = 80e9;         %Shear modulus
    r = 0.001;        %Cross-sectional radius
    rho = 8000;       %Density
    g = [9.81; 0; 0]; %Gravitational acceleration
    L = 0.5;          %Length (before strain)

    %Dependent Parameters
    A = pi*r^2;   %Cross-sectional area
    I = pi*r^4/4; %Area moment of inertia
    J = 2*I;      %Polar moment of inertia

    Kse = diag([G*A, G*A, E*A]); %Stiffness matrices
    Kbt = diag([E*I, E*I, G*J]);

    %Measured base force and moment
    n0 = [0; 1; 0];
    m0 = [0; 0; 0];

    %Arbitrary base frame assignment
    p0 = [0;0;0];
    R0 = eye(3);

    %Numerical Integration
    y0 = [p0; reshape(R0,9,1); n0; m0]; %Combine states into single state vector
    [s,y] = ode45(@RodODE,[0 L],y0);    %Solve IVP with numerical integration
    
    %Visualization
    plot3( y(:,1), y(:,2), y(:,3) )
    axis([-L/2 L/2  -L/2 L/2  0 L])
    daspect([1 1 1])
    grid on
    title('Rod IVP Solution')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')

    %Subfunctions
    function ys = RodODE(s,y) %State vector derivative function
        %Unpack state vector
        R = reshape(y(4:12),3,3);
        n = y(13:15);
        m = y(16:18);

        %Constitutive equation
        v = Kse^-1*R.'*n + [0;0;1];
        u = Kbt^-1*R.'*m;

        %Static Cosserat rod equations - system of nonlinear ODEs
        ps = R*v;
        Rs = R*hat(u);
        ns = -rho*A*g;
        ms = -cross(ps,n);

        %Pack state vector derivative
        ys = [ps; reshape(Rs,9,1); ns; ms];
    end

    function skew_symmetric_matrix = hat(y)
        skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                                  y(3)   0   -y(1) ;
                                 -y(2)  y(1)   0  ];
    end
end