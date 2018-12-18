function CSG_Inverse_Kinematics
    %Properties
    E=200e9;  G=80e9;  rad=0.001;  rho=8000;  g=[0;0;-9.81];  ee_mass=0.1;  R0=eye(3);
    A=pi*rad^2;  I=pi*rad^4/4;  J=2*I;  Kse=diag([G*A,G*A,E*A]);  Kbt=diag([E*I,E*I,G*J]);

    %Given End Effector Pose and Wrench
    F = ee_mass*g;   M = [0;0;0];
    pE = [0; 0; 0.4];  bend = 10*pi/180;  RE = [ cos(bend)  0  sin(bend) ;
                                                     0      1      0     ;
                                                -sin(bend)  0  cos(bend)];
    %Hole Pattern (vector operations)
    scrib_R = 0.087;  alpha1 = 100*pi/180;  alpha2 = 120*pi/180 - alpha1;  i = 1:6;
    theta_B = (-alpha2 + (i-mod(i,2))*alpha2 + (i-1-mod(i-1,2))*alpha1)/2;
    theta_E = (-alpha1 + (i-mod(i,2))*alpha1 + (i-1-mod(i-1,2))*alpha2)/2;
    p0 = scrib_R*[cos(theta_B); sin(theta_B); zeros(1,6)];
    r  = scrib_R*[cos(theta_E); sin(theta_E); zeros(1,6)];
    
    hat=@(y)[0,-y(3),y(2);y(3),0,-y(1);-y(2),y(1),0];  inv_hat_xy=@(y)[y(3,2);y(1,3)];
    
    %Main Simulation
    init_guess = [zeros(30,1); pE(3)*ones(6,1)]; %Initial guess
    global p
    fsolve(@CSG_BVP_Function, init_guess); %Solve CSG BVP with shooting method
    
    %Visualization
    plot3( p{1}(1,:), p{1}(2,:), p{1}(3,:), 'b' );  hold on;
    plot3([p{6}(1,end) p{1}(1,end)],[p{6}(2,end) p{1}(2,end)],[p{6}(3,end) p{1}(3,end)],'r')
    for i = 2 : 6
        plot3(p{i}(1,:),p{i}(2,:),p{i}(3,:),'b');  ee_line = [p{i-1}(:,end), p{i}(:,end)];
        plot3(ee_line(1,:),ee_line(2,:),ee_line(3,:),'r')
    end
    hold off;  axis([-0.25 0.25  -0.25 0.25  0 0.5]);  daspect([1 1 1]);
    title('CSG BVP Solution'); grid on; xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
    
    %Subfunctions
    function E = CSG_BVP_Function(G)
        EF = F;
        EM = M; %Begin summing forces

        for i = 1 : 6 %Loop over each compliant robot link
            n0_range = 1+5*(i-1) : 3+5*(i-1);   m0_range = 4+5*(i-1) : 5*i;
            n0 = G(n0_range);  m0 = [G(m0_range); 0];  L = G(30+i);
            y0 = [p0(:,i); reshape(R0,9,1); n0; m0];

            [~,y] = ode45(@RodODE,[0 L],y0); %Numerically integrate this rod
            p{i} = y(:,1:3)'; %Store centerlines in cells for plotting

            pL_shot = y(end,1:3)';
            RL_shot = reshape(y(end,4:12),3,3);
            nL = y(end,13:15)';
            mL = y(end,16:18)';

            Ep = pL_shot - (pE + RE*r(:,i)); %Calculate geometric error
            ER = inv_hat_xy( RL_shot'*RE - RL_shot*RE' );
            geometric_error_range = 1+5*(i-1) : 5*i;
            E(geometric_error_range) = [Ep; ER];
            
            EF = EF - nL;
            EM = EM - mL - cross( RE*r(:,i), nL ); %Continue summing forces
        end
        
        E(31:33) = EF;
        E(34:36) = EM; %Force and moment summation are complete
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
end