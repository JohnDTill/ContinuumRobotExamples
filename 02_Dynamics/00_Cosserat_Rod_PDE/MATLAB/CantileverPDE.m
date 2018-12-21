function CantileverPDE
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];
    global p R j n m v u q w vt ut qt wt vh uh qh wh %Make vars available in whole program
    %Parameters
    L = 0.4;         %Length (before strain)
    N = 40;          %Spatial resolution
    E = 207e9;       %Young's modulus
    r = 0.0012;      %Cross-section radius
    rho = 8000;      %Density
    g = [-9.81;0;0]; %Gravity
    Bse = zeros(3);  %Material damping coefficients - shear and extension
    Bbt = zeros(3);  %Material damping coefficients - bending and torsion
    C = zeros(3);    %Square-law-drag damping coefficients
    dt = 0.002;      %Time step
    alpha = -0.48;   %BDF-alpha parameter
    STEPS = 50;      %Number of timesteps to completion
    vstar = @(s)[0;0;1]; %Value of v when static and absent loading
    ustar = @(s)[0;0;0]; %Precurvature
    %Boundary Conditions
    for i = 1 : STEPS
        p{i,1} = [0;0;0]; %Clamped base
        R{i,1} = eye(3);
        q{i,1} = [0;0;0];
        w{i,1} = [0;0;0];
    end
    nL = 0.2*g;  %Start with a weight hung at the tip
    mL = [0;0;0];

    %Dependent Parameter Calculations
    A = pi*r^2;                                 %Cross-sectional area
    J = diag([pi*r^4/4  pi*r^4/4  pi*r^4/2]);   %Inertia
    G = E/( 2*(1+0.3) );                        %Shear modulus
    Kse = diag([G*A, G*A, E*A]);                %Stiffness matrix - shear and extension
    Kbt = diag([E*J(1,1), E*J(2,2), G*J(3,3)]); %Stiffness matrix - bending and torsion
    ds = L/(N-1);                               %Grid distance (before strain)
    c0 = (1.5 + alpha) / ( dt*(1+alpha) );      %BDF-alpha coefficients
    c1 = -2/dt;
    c2 = (0.5 + alpha) / ( dt*(1+alpha) );
    d1 = alpha / (1+alpha);
    
    %Main Simulation
    i = 1;
    fsolve(@staticIVP, zeros(6,1)); %Solve static BVP w/ shooting method
    applyStaticBDFalpha();
    visualize();
    nL = [0;0;0]; %Tip weight is released
    
    for i = 2 : STEPS
        fsolve(@dynamicIVP, [n{i-1,1}; m{i-1,1}]); %Solve semi-discretized PDE w/ shooting
        applyDynamicBDFalpha();
        visualize();
    end

    %Function Definitions
    function applyStaticBDFalpha()
        for j = 1 : N-1
            vh{i+1,j} = (c1+c2)*v{i,j};
            uh{i+1,j} = (c1+c2)*u{i,j};
            qh{i+1,j} = [0;0;0];
            wh{i+1,j} = [0;0;0];
            q{i,j} = [0;0;0];
            w{i,j} = [0;0;0];
        end
    end

    function applyDynamicBDFalpha()
        for j = 1 : N-1        
            vh{i+1,j} = c1*v{i,j} + c2*v{i-1,j} + d1*vt{i,j};
            uh{i+1,j} = c1*u{i,j} + c2*u{i-1,j} + d1*ut{i,j};
            qh{i+1,j} = c1*q{i,j} + c2*q{i-1,j} + d1*qt{i,j};
            wh{i+1,j} = c1*w{i,j} + c2*w{i-1,j} + d1*wt{i,j};
        end
    end

    function E = staticIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-1
            [ps, Rs, ns, ms, v{i,j}, u{i,j}] = staticODE(p{i,j},R{i,j},n{i,j},m{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
        end
        E = [ n{i,N} - nL;  m{i,N} - mL ];
    end

    function E = dynamicIVP(G)
        n{i,1} = G(1:3);
        m{i,1} = G(4:6);

        %Euler's method
        for j = 1 : N-1
            [ps, Rs, ns, ms, qs, ws, v{i,j}, u{i,j}, ...
                vt{i,j}, ut{i,j}, qt{i,j}, wt{i,j}] = ...
                dynamicODE(p{i,j},R{i,j},n{i,j},m{i,j},q{i,j},w{i,j});
            p{i,j+1} = p{i,j} + ds*ps;
            R{i,j+1} = R{i,j} + ds*Rs;
            n{i,j+1} = n{i,j} + ds*ns;
            m{i,j+1} = m{i,j} + ds*ms;
            q{i,j+1} = q{i,j} + ds*qs;
            w{i,j+1} = w{i,j} + ds*ws;
        end
        E = [n{i,N} - nL;  m{i,N} - mL];
    end

    function [ps, Rs, ns, ms,  v, u] = staticODE(p,R,n,m)
        v = Kse\R'*n + vstar(ds*(j-1));
        u = Kbt\R'*m + ustar(ds*(j-1));

        ps = R*v;
        Rs = R*hat(u);
        ns = -rho*A*g;
        ms = -hat(ps)*n;
    end

    function [ps,Rs,ns,ms,qs,ws, v,u,vt,ut,qt,wt] = dynamicODE(p,R,n,m,q,w)
        v = (Kse + c0*Bse)\(R'*n + Kse*vstar(ds*(j-1)) - Bse*vh{i,j});
        u = (Kbt + c0*Bbt)\(R'*m + Kbt*ustar(ds*(j-1)) - Bbt*uh{i,j});
        vt = c0*v + vh{i,j};
        ut = c0*u + uh{i,j};
        qt = c0*q + qh{i,j};
        wt = c0*w + wh{i,j};
        f = -R*C*q.*abs(q) + rho*A*g;

        ps = R*v;
        Rs = R*hat(u);
        ns = rho*A*R*(hat(w)*q + qt) - f;
        ms = rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
    end

    function visualize()
        for j = 1 : N,  x(j) = p{i,j}(1);  z(j) = p{i,j}(3);  end
        plot(z,x);  axis([0 1.1*L  -0.55*L 0.55*L]);  daspect([1 1 1]);
        title('Cantilever Rod');  xlabel('z (m)');  ylabel('x (m)');
        grid on;  drawnow;  pause(0.05);
    end
end