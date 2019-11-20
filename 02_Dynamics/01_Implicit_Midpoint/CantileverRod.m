function CantileverRod %Section 3.1.4
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0]; global y z
    %Parameters
    L = 0.4;  E = 207e9;  r = 0.0012;  rho = 8000;  g = [-9.81;0;0];
    Bse = zeros(3);  Bbt = zeros(3);  C = zeros(3);  vstar = [0;0;1];
    N = 40;  dt = 0.002;  STEPS = 50;  M_tip = [0;0;0];  F_tip = 0.2*g; %Start with tip load
    %Boundary Conditions
    p0 = [0;0;0];  R0 = eye(3);  q0 = [0;0;0];  w0 = [0;0;0];
    %Dependent Parameter Calculations
    A = pi*r^2;  G = E/( 2*(1+0.3) );  ds = L/(N-1);
    Ixx = pi*r^4/4;  Iyy = Ixx;  Izz = 2*Ixx;  J = diag([Ixx  Iyy  Izz]);
    Kse = diag([G*A, G*A, E*A]);   Kbt = diag([E*Ixx, E*Iyy, G*Izz]);
    %Main Simulation
    i = 1;  G = fsolve(@getResidual, zeros(6,1)); %Solve static BVP
    z(:,N) = solveStaticConstitutiveLaw(y(:,N));
    ys(:,N) = f(y(:,end),zeros(24,1),z(:,end),zeros(6,1)); %Distal ys needed
    F_tip = [0;0;0]; %Tip weight is released
    for i = 2 : STEPS
        y_old = y;  ys_old = ys;  z_old = z;  visualize();
        G = fsolve(@getResidual, G); %Solve semi-discretized BVP
        ys(:,N) = implicitMidptODE(y(:,end),y_old(:,end),ys_old(:,end),z_old(:,end));
    end
    %Function Definitions
    function visualize()
        plot(y(3,:),y(1,:));  title('Cantilever Rod');  xlabel('z (m)');  ylabel('x (m)');
        axis([0 1.1*L  -0.55*L 0.55*L]);  grid on;  daspect([1 1 1]);  drawnow;
    end
    function E = getResidual(G)
        n0 = G(1:3);  m0 = G(4:6); %Reaction force and moment are guessed
        y(:,1) = [p0; reshape(R0,9,1); n0; m0; q0; w0];
        for j = 1 : N-1 %Euler Integration
            if i == 1 %First time step is static
                z(:,j) = solveStaticConstitutiveLaw(y(:,j));
                ys(:,j) = f(y(:,j),zeros(24,1),z(:,j),zeros(6,1));
            else      %Next time steps use PDE semi-discretization
                [ys(:,j),z(:,j)]=implicitMidptODE(y(:,j),y_old(:,j),ys_old(:,j),z_old(:,j));
            end
            y(:,j+1) = y(:,j) + ds*ys(:,j); %Euler's Method
        end
        nL = y(13:15,N);  mL = y(16:18,N);  E = [F_tip-nL;  M_tip-mL];
    end
    function [ys,z] = implicitMidptODE(y,y_old,ys_old,z_old) %Equation (3.7)
        z = solveDynamicConstitutiveLaw(y,y_old,z_old);
        ys = -ys_old + 2*f( (y+y_old)/2, (y-y_old)/dt, (z+z_old)/2, (z-z_old)/dt );
    end
    function z = solveStaticConstitutiveLaw(y)
        R = reshape(y(4:12),3,3);  n = y(13:15);  m = y(16:18);
        v = Kse\R'*n + vstar;  u = Kbt\R'*m;  z = [v;u];
    end
    function z = solveDynamicConstitutiveLaw(y,y_old,z_old) %Equation (3.8)
        Rbar = (reshape(y(4:12),3,3) + reshape(y_old(4:12),3,3))/2;
        nbar = (y(13:15) + y_old(13:15))/2;  mbar = (y(16:18) + y_old(16:18))/2;
        v_prev = z_old(1:3);  u_prev = z_old(4:6);
        v = (Kse/2+Bse/dt) \ ((-Kse/2+Bse/dt)*v_prev + Rbar'*nbar + Kse*vstar);
        u = (Kbt/2+Bbt/dt) \ ((-Kbt/2+Bbt/dt)*u_prev + Rbar'*mbar);
        z = [v;u];
    end
    function ys = f(y,yt,z,zt) %Equation (3.5)
        R = reshape(y(4:12),3,3);  n = y(13:15);  m = y(16:18);
        q = y(19:21);  w = y(22:24);  qt = yt(19:21);  wt = yt(22:24);
        v = z(1:3);  u = z(4:6);  vt = zt(1:3);  ut = zt(4:6);
        ps = R*v;
        Rs = R*hat(u);
        ns = R*(rho*A*(hat(w)*q + qt) + C*q.*abs(q)) - rho*A*g;
        ms = rho*R*(hat(w)*J*w + J*wt) - hat(ps)*n;
        qs = vt - hat(u)*q + hat(w)*v;
        ws = ut - hat(u)*w;
        ys = [ps; reshape(Rs,9,1); ns; ms; qs; ws];
    end
end