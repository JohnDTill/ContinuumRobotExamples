%plot 3d tendonRobot
%
function [] = Visualize(centerline,tendonlines,num_tendons, disks,num_disks)
L=0.8;
p1=plot3(centerline(1,:),centerline(2,:),centerline(3,:));
p1.LineWidth=3;
p1.Color='b';

axis([-0.3*L 0.3*L  -0.3*L 0.3*L 0 1.05*L]);
title('');  xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
daspect([1 1 1]); %Lock aspect ratio
hold on;

%plot tendons%
for j = 1 : num_tendons
    plot3(tendonlines(3*j-2, :),tendonlines(3*j-2+1,:)...
        ,tendonlines(3*j,:),'Color','#000000','LineWidth', ...
        1);
end

%-plot disks---%
for i = 1 : num_disks
    R=disks(1:3, 4*i-3:4*i-1);
    P=disks(1:3, 4*i);
    u1=P+R*[0;0;-0.001/2];
    u2=P+R*[0;0;0.001/2];
    [cycleX,cycleY,cycleZ] = myCylinder(u1,u2,0.025);
    plot3(cycleX(1,:),cycleY(1,:),cycleZ(1,:),'black');
    %fill3(cycleX(1,:),cycleY(1,:),cycleZ(1,:),'yellow') %bottom
    %fill3(cycleX(2,:),cycleY(2,:),cycleZ(2,:),'yellow') %top
    %surf(cycleX,cycleY,cycleZ,'facecolor','yellow','EdgeColor','none') %Cylindrical surface

end
grid on;  drawnow;
%pause(0.05);
hold off;

%% function
    function [X,Y,Z] = myCylinder(u1,u2,r)
        
        n=u2-u1;
        theta=(0:2*pi/100:2*pi)'; %theta  0 to 2*pi
        a=cross(n,[1 0 0]);       
        if ~any(a) 
            a=cross(n,[0 1 0]);
        end
        b=cross(n,a); 
        a=a/norm(a); 
        b=b/norm(b);  

        
        c1=u1(1)*ones(size(theta,1),1);
        c2=u1(2)*ones(size(theta,1),1);
        c3=u1(3)*ones(size(theta,1),1);
        x1=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta); %圆上各点的x坐标
        y1=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta); %圆上各点的y坐标
        z1=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta); %圆上各点的z坐标

        
        c1=u2(1)*ones(size(theta,1),1);
        c2=u2(2)*ones(size(theta,1),1);
        c3=u2(3)*ones(size(theta,1),1);
        x2=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta); 
        y2=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta); 
        z2=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta); 

        X(1,:)=x1(:);
        X(2,:)=x2(:);
        Y(1,:)=y1(:);
        Y(2,:)=y2(:);
        Z(1,:)=z1(:);
        Z(2,:)=z2(:);
    end
end

