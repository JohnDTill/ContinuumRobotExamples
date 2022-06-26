%plot 3d tendonRobot
%
function [] = Visualize(centerline,tendonlines,num_tendons, disks,num_disks)
    
    p1=plot3(centerline(1,:),centerline(2,:),centerline(3,:));
    p1.LineWidth=3;
    p1.Color='b';
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
        if i==(num_disks*0.5) %驱动盘单独画
            fill3(cycleX(1,:),cycleY(1,:),cycleZ(1,:),'green') %底面
            fill3(cycleX(2,:),cycleY(2,:),cycleZ(2,:),'green') %顶面
            surf(cycleX,cycleY,cycleZ,'facecolor','yellow','EdgeColor','none')    %圆柱表面
        else
            fill3(cycleX(1,:),cycleY(1,:),cycleZ(1,:),'yellow') %底面
            fill3(cycleX(2,:),cycleY(2,:),cycleZ(2,:),'yellow') %顶面
            surf(cycleX,cycleY,cycleZ,'facecolor','yellow','EdgeColor','none')    %圆柱表面
        end
    
    end
    % title('Cantilever Rod');  xlabel('z (m)');  ylabel('x (m)');
    %axis([-1.2*L 1.2*L  -1.2*L 1.2*L 0 2.1*L])
    title('');  xlabel('x (m)');  ylabel('y (m)');  zlabel('z (m)');
    daspect([1 1 1]); %锁定纵横比
    grid on;  drawnow;
    %pause(0.05);
    hold off;
end

