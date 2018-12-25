%参数设置和模型建立
phi=pi/4;%初始纬度
omega=7.2916*10^-5;%地球自转角速度
g=9.84;%重力加速度
%ac=100*10^-6*g;
ar=50*10^-6*g;%加速度计初始偏差和随机偏差
%pc=0.02*pi/180/3600;   
pr=0.01*pi/180/3600;%陀螺仪常值漂移和随机漂移 
deg=1/180*pi;%一度对应的的弧度（初始失准角）

xr0=[0.1 0.1 deg deg deg zeros(1,5)]';%状态初值，真实值
xe0=zeros(1,10)';%状态估计初值
Q=diag([ar^2 ar^2 pr^2 pr^2 pr^2 zeros(1,5)]);%系统噪声协方差矩阵
P0=diag([0.1^2 0.1^2 deg^2 deg^2 deg^2 zeros(1,5)]);%系统误差协方差矩阵
R=diag([0.1^2 0.1^2]);%观测噪声协方差矩阵

Ou=omega*sin(phi);On=omega*cos(phi);
F=[0 2*Ou 0 -g 0;
    -2*Ou 0 g 0 0;
    0 0 0 Ou -On;
    0 0 -Ou 0 0;
    0 0 On 0 0];
T=eye(5);%初始这个矩阵可以取为单位阵，即认为两个坐标系重合，之后需要使用估计出来的失准角进行修正
A=[F T;
    zeros(5,5) zeros(5,5)];%10*10
B=eye(10);%10*10
C=[eye(2) zeros(2,8)];%2*10
D=zeros(2,10);%这里虽然设了B和D实际上都没有用，只是为了能够使用ss函数才设成D=zeros(2,10)
G=ss(A,B,C,D);%创建状态空间模型
t=1;
Gd=c2d(G,t);%c2d()函数的作用是将s域（频域）的表达式转化成z域（z变换）的表达式（零阶保持器）
PHI=Gd.a;%返回z域的A矩阵
N=1000;
X=zeros(10,N);%估计输出

%开始计算了
xrk_1=xr0;%状态的真实值
xk_1=xe0;%状态估计值
Pk_1=P0;%系统误差协方差矩阵
for i=1:N
    Z=C*xrk_1+(mvnrnd([0;0],R,1))';%当前时刻观测值Z
    Pkk=PHI*Pk_1*PHI'+Q;%一步预测协方差矩阵
    %K=Pkk*C'*inv(C*Pkk*C'+R);速度更慢准确性更低
    K=(Pkk*C')/(C*Pkk*C'+R);%滤波增益矩阵，计算卡尔曼增益K
    Pk=(eye(10)-K*C)*Pkk;%协方差矩阵更新
    xk=PHI*xk_1;%状态一步预测
    xk=xk+K*(Z-C*xk);%状态估计/状态更新
    %xk=xk;
    X(:,i)=xk;
    xk_1=xk;
    Pk_1=Pk;
    %根据失准角修正T
   xrk=A*xrk_1+(mvnrnd(zeros(10,1),Q,1))';
   xrk_1=xrk;
end

t=1:N;
%plot(t,X(1,:),'-r',t,X(2,:),'-k'); 
%title('1');
figure;
plot(t,X(3,:)*3600*180/pi,'-k');
title('东向失准角');
figure;
plot(t,X(4,:)*3600*180/pi,'-k');%化成秒（角）输出
title('北向失准角');
figure;
plot(t,X(5,:)*3600*180/pi,'-k');
title('方位失准角');
