function CantileverRod
    %Parameters - Section 2.1
    L = 0.4;
    N = 40;
    E = 207e9;
    r = 0.0012;
    rho = 8000;
    C = zeros(3);
    g = 10*[-9.81;0;0]; %10x greater for effect
    Bse = zeros(3);
    Bbt = zeros(3);
    del_t = 0.005;
    STEPS = 50;
    F_tip = [0;0;0];
    M_tip = [0;0;0];
    vstar = [0;0;1];
    %Boundary Conditions - Section 2.4
    p0 = [0;0;0];
    h0 = [1;0;0;0];
    q0 = [0;0;0];
    w0 = [0;0;0];
    %Dependent Parameter Calculations
    A = pi*r^2;
    G = E/( 2*(1+0.3) );
    J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);
    Kse = diag([G*A, G*A, E*A]);
    Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]);
    ds = L/(N-1);
    %BDF2 Coefficients
    c0 = 1.5/del_t;
    c1 = -2/del_t;
    c2 = 0.5/del_t;
    %Expressions extracted from simulation loop
    Kse_plus_c0_Bse_inv = (Kse+c0*Bse)^-1;
    Kbt_plus_c0_Bbt_inv = (Kbt+c0*Bbt)^-1;
    Kse_vstar = Kse*vstar;
    rhoA = rho*A;
    rhoAg = rho*A*g;
    rhoJ = rho*J;
 
    %Initialize to straight configuration
    % y and z are general variables as in Eq (12)
    % y:=[p;h;n;m;q;w] and z:=[v;u]
    y = [zeros(2,N); linspace(0,L,N); zeros(16,N)];
    z = [zeros(2,N); ones(1,N); zeros(3,N)];
    y_prev = y;
    z_prev = z;
 
    %Main Simulation Loop - Section 2.4
    visualize();
    G = zeros(6,1); %Shooting method initial guess
    for i = 2 : STEPS
        %Set history terms - Eq (5)
        yh = c1*y+c2*y_prev;
        zh = c1*z+c2*z_prev;
        y_prev = y;
        z_prev = z;
        %Midpoints are linearly interpolated for RK4
        yh_int = 0.5*(yh(:,1:end-1) + yh(:,2:end));
        zh_int = 0.5*(zh(:,1:end-1) + zh(:,2:end));
        %Shooting method solver call
        G = fsolve(@getResidual, G);
 
        visualize();
    end
 
    %Function Definitions
    function visualize()
        plot(y(3,:),y(1,:));
        title('Cantilever Rod');
        xlabel('z (m)');
        ylabel('x (m)');
        axis([0 1.1*L  -0.55*L 0.55*L]);
        grid on;
        daspect([1 1 1]);
        drawnow;
    end
 
    function E = getResidual(G)
        %Reaction force and moment are guessed
        n0 = G(1:3);
        m0 = G(4:6);
        y(:,1) = [p0; h0; n0; m0; q0; w0];
        %Fourth-Order Runge-Kutta Integration
        for j = 1 : N-1
            yj = y(:,j);
            yhj_int = yh_int(:,j);
            [k1,z(:,j)]=ODE(yj,yh(:,j),zh(:,j));
            [k2,~]=ODE(yj+k1*ds/2,yhj_int,zh_int(:,j));
            [k3,~]=ODE(yj+k2*ds/2,yhj_int,zh_int(:,j));
            [k4,~]=ODE(yj+k3*ds,yh(:,j),zh(:,j+1));
            y(:,j+1) = yj + ds*(k1 + 2*(k2+k3) + k4)/6;
            %y(:,j+1) = yj + ds*k1; %Euler's Method
        end
        %Cantilever boundary conditions
        nL = y(8:10,N);
        mL = y(11:13,N);
        E = [F_tip-nL;  M_tip-mL];
    end
 
    function [ys,z] = ODE(y,yh,zh)
        h = y(4:7);
        n = y(8:10);
        m = y(11:13);
        q = y(14:16);
        w = y(17:19);
        vh = zh(1:3);
        uh = zh(4:6);
         
        %Quaternion to Rotation - Eq (10)
        h1=h(1);
        h2=h(2);
        h3=h(3);
        h4=h(4);
        R = eye(3) + 2/(h'*h) * ...
            [-h3^2-h4^2  , h2*h3-h4*h1,  h2*h4+h3*h1;
            h2*h3+h4*h1, -h2^2-h4^2 ,  h3*h4-h2*h1;
            h2*h4-h3*h1, h3*h4+h2*h1, -h2^2-h3^2  ];
         
        %Solved Constutitve Law - Eq (6)
        v=Kse_plus_c0_Bse_inv*(R'*n+Kse_vstar-Bse*vh);
        u=Kbt_plus_c0_Bbt_inv*(R'*m-Bbt*uh);
        z=[v;u];
         
        %Time Derivatives - Eq (5)
        yt = c0*y + yh;
        zt = c0*z + zh;
         
        vt = zt(1:3);
        ut = zt(4:6);
        qt = yt(14:16);
        wt = yt(17:19);
         
        %Weight and Square-Law-Drag - Eq (3)
        f = rhoAg - R*C*q.*abs(q);
         
        %Rod State Derivatives - Eq (7)
        ps = R*v;
        ns = rhoA*R*(cross(w,q) + qt) - f;
        ms = R*(cross(w,rhoJ*w) +rhoJ*wt)-cross(ps,n);
        qs = vt - cross(u,q) + cross(w,v);
        ws = ut - cross(u,w);
         
        %Quaternion Derivative - Eq (9)
        hs = [ 0  , -u(1), -u(2), -u(3);
            u(1),   0  ,  u(3), -u(2);
            u(2), -u(3),   0  ,  u(1);
            u(3),  u(2), -u(1),   0  ] * h/2;
         
        ys = [ps;hs;ns;ms;qs;ws];
    end
end