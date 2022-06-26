function [X,Y,Z] = myCylinder(u1,u2,r)
%通过空间中两点和半径得到圆柱底和顶的圆点
n=u2-u1;
theta=(0:2*pi/100:2*pi)'; %theta角从0到2*pi
a=cross(n,[1 0 0]);       %n与i叉乘，求取a向量
if ~any(a) %如果a为零向量，将n与j叉乘   %any检测矩阵中是否有非零元素，如果有，则返回1，否则，返回0
    a=cross(n,[0 1 0]);
end
b=cross(n,a); %求取b向量
a=a/norm(a);  %单位化a向量
b=b/norm(b);  %单位化b向量
%a,b是满足既垂直于n，又互相垂直的任意单位向量

%底面圆方程
c1=u1(1)*ones(size(theta,1),1);
c2=u1(2)*ones(size(theta,1),1);
c3=u1(3)*ones(size(theta,1),1);
x1=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta); %圆上各点的x坐标
y1=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta); %圆上各点的y坐标
z1=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta); %圆上各点的z坐标

%顶面圆方程
c1=u2(1)*ones(size(theta,1),1);
c2=u2(2)*ones(size(theta,1),1);
c3=u2(3)*ones(size(theta,1),1);
x2=c1+r*a(1)*cos(theta)+r*b(1)*sin(theta); %圆上各点的x坐标
y2=c2+r*a(2)*cos(theta)+r*b(2)*sin(theta); %圆上各点的y坐标
z2=c3+r*a(3)*cos(theta)+r*b(3)*sin(theta); %圆上各点的z坐标

X(1,:)=x1(:);
X(2,:)=x2(:);
Y(1,:)=y1(:);
Y(2,:)=y2(:);
Z(1,:)=z1(:);
Z(2,:)=z2(:);
end

