function TendonRobotDisplacement
    %Independent Parameters
    E = 200e9;
    shear_modulus = 80e9;
    rad = 0.001;
    rho = 8000;
    g = [9.81; 0; 0];
    L = 0.5;
    num_tendons = 4;
    
    tendon_offset = 0.01506;
    rot = @(theta) tendon_offset*[cos(theta); sin(theta); 0];
    r = {rot(0), rot(pi/2), rot(pi), rot(pi*3/2)}; %Tendons 90 deg apart
    
    %Boundary conditions
    p0 = [0;0;0];
    R0 = eye(3);

    q=[10e-3; 0; -10e-3; 0];
    l_star = [L; L; L; L];
    C=[1e-4; 1e-4; 1e-4; 1e-4];
    tau = [0; 0; 0; 0];

    %Dependent parameter calculations
    area = pi*rad^2;
    I = pi*rad^4/4;
    J = 2*I;
    Kse = diag([shear_modulus*area, shear_modulus*area, E*area]);
    Kbt = diag([E*I, E*I, shear_modulus*J]);
    
   
    %Solve with shooting method
    init_guess = zeros(6+num_tendons,1);
    global Y; %Forward declaration for future scoping rule changes

    fsolve(@shootingFunction, init_guess);

    %export data
    centerline=Y(:,1:3).';
    %dlmwrite('centerline.dat', centerline, 'Delimiter', '\t');
    
    tendonlines = zeros(3*num_tendons, size(Y,1));
    for i = 1 : size(Y,1)
        p = Y(i,1:3).';
        R = reshape(Y(i,4:12).',3,3);
        for j = 1 : num_tendons
            tendonlines(3*j-2 : 3*j, i) = p + R*r{j};
        end
    end
    dlmwrite('tendonlines.dat', tendonlines, 'Delimiter', '\t');
    
    num_disks = 9;
    disks = zeros(3,4*num_disks);
    for i = 1 : num_disks
        j = round(size(Y,1) * i / num_disks);
        p = Y(j,1:3).';
        R = reshape(Y(j,4:12).',3,3);
        disks(1:3, 4*i-3:4*i-1) = R;
        disks(1:3, 4*i) = p;
    end
    dlmwrite('disks.dat', disks, 'Delimiter', '\t');
    

    %plot TendonRobot with Matlab
    Visualize(centerline,tendonlines,num_tendons, disks,num_disks)

    disp(['Simulation finished.', newline, 'For visualization, move ', ...
         '''.dat'' files to Blender subfolder and run Blender script.',...
         newline])
    
    %Subfunctions
    function y_s = cosseratTendonRobotOde(s,y)
        %Unpack state vector
        R = reshape(y(4:12),3,3);
        v = y(13:15);
        u = y(16:18);
        
        %Setup tendon linear system
        a = zeros(3,1);
        b = zeros(3,1);
        A = zeros(3,3);
        G = zeros(3,3);
        H = zeros(3,3);
        pb_s_norm=y(19:19+num_tendons-1);

        for i = 1 : num_tendons
            pb_si = cross(u,r{i}) + v;
            pb_s_norm(i) = norm(pb_si);
            A_i = -hat(pb_si)^2 * (tau(i)/pb_s_norm(i)^3);
            G_i = -A_i * hat(r{i});
            a_i = A_i * cross(u,pb_si);
            
            a = a + a_i;
            b = b + cross(r{i}, a_i);
            A = A + A_i;
            G = G + G_i;
            H = H + hat(r{i})*G_i;
        end
        
        K = [A+ Kse,     G   ;
                G.' ,  H + Kbt];
            
        nb = Kse*(v - [0;0;1]);
        mb = Kbt*u;
        
        rhs = [-cross(u,nb) - R.'*rho*area*g - a;
               -cross(u,mb) - cross(v,nb) - b];
        
        %Calculate ODE terms
        p_s = R*v;
        R_s = R*hat(u);
        vs_and_us = K\rhs;
        
        %Pack state vector derivative
        y_s = [p_s; reshape(R_s,9,1); vs_and_us;pb_s_norm];
    end

    function distal_error = shootingFunction(guess)
        n0 = guess(1:3);
        v0 = Kse\n0 + [0;0;1];
        u0 = guess(4:6);
        
        tau=max(guess(7:7+num_tendons-1),0); %拉力,和0比取最大
        slack=-min(guess(7:7+num_tendons-1),0);

        y0 = [p0; reshape(R0,9,1); v0; u0; q]; %初始y0
        
        [~,Y] = ode45(@cosseratTendonRobotOde, linspace(0,L), y0);
        
        %Find the internal forces in the backbone prior to the final plate
        vL = Y(end,13:15).';
        uL = Y(end,16:18).';
        
        nb = Kse*(vL - [0;0;1]);
        mb = Kbt*uL;
        
        %Find the equilibrium error at the tip, considering tendon forces
        force_error = -nb;
        moment_error = -mb;
        for i = 1 : num_tendons
            pb_si = cross(uL,r{i}) + vL;
            Fb_i = -tau(i)*pb_si/norm(pb_si);
            force_error = force_error + Fb_i;
            moment_error = moment_error + cross(r{i}, Fb_i);
        end
        
        %Find the length violation error
        integrated_lengths=Y(end,19:19+num_tendons-1)';
        stretch = l_star.*(C.*tau);
        length_error = integrated_lengths + slack - (l_star + stretch);
        distal_error = [force_error; moment_error;length_error];
    end

    function skew_symmetric_matrix = hat(y)
        skew_symmetric_matrix = [  0   -y(3)  y(2) ;
                                  y(3)   0   -y(1) ;
                                 -y(2)  y(1)   0  ];
    end
end
