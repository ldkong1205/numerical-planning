clear all;
format long;

global aa alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T mu_p mu_v qP qM qDp qDm qDDp qDDm m Pinfty Minfty gamma m n znn_W k;
load mvn_SYSdata;
load mvn_SNDdata;

Time=0;                             %Time清零
t_last=0;
t_time=0;
[z,n]=size(qAll);
q_Time=zeros(z,1);

[Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q0);          %计算雅可比矩阵
ix0=Ppx(7);iy0=Ppy(7);iz0=Ppz(7);
for jj=1:length(t),
   qjj=qAll(jj,:)';%simulated end effector position
   [Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qjj); 
   j0px(jj,1)=Ppx(1);
   j1px(jj,1)=Ppx(2);
   j2px(jj,1)=Ppx(3);
   j3px(jj,1)=Ppx(4);
   j4px(jj,1)=Ppx(5);
   j5px(jj,1)=Ppx(6);
   posx(jj,1)=Ppx(7);%joint6
   j0py(jj,1)=Ppy(1);
   j1py(jj,1)=Ppy(2);
   j2py(jj,1)=Ppy(3);
   j3py(jj,1)=Ppy(4);
   j4py(jj,1)=Ppy(5);
   j5py(jj,1)=Ppy(6);
   posy(jj,1)=Ppy(7);%joint6
   j0pz(jj,1)=Ppz(1);
   j1pz(jj,1)=Ppz(2);
   j2pz(jj,1)=Ppz(3);
   j3pz(jj,1)=Ppz(4);
   j4pz(jj,1)=Ppz(5);
   j5pz(jj,1)=Ppz(6);
   posz(jj,1)=Ppz(7);%joint6
   %-------------------------------------------------
   dqjj=dqAll(jj,:)';%simulated end effector velocity
  [J,DJ]=ZYNjdj(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qjj,dqjj); 
   %-----------------
   dpos=J*dqjj;
   dposx(jj,1)=dpos(1);
   dposy(jj,1)=dpos(2);
   dposz(jj,1)=dpos(3);
   
   ddqjj=ddqAll(jj,:)';
   ddpos=J*ddqjj+DJ*dqjj;
   ddposx(jj,1)=ddpos(1);
   ddposy(jj,1)=ddpos(2);%simulated end-effector acceleration
   ddposz(jj,1)=ddpos(3);

t_time=t(jj)-Time*T;                    %每个周期内的相对于周期起点的运行时间，对该变量初始化计算出处于本周期中的哪一个时刻
if t_time>T;
   q_Time(jj,1)=1;
   Time=Time+1;
   t_time=t_time-T;
end
   
% phi_sin=2*pi*sin(0.5*pi*t_time/T);
% phi=phi_sin*sin(0.5*pi*t_time/T);
% phiDot=phi_sin*pi*cos(0.5*pi*t_time/T)/T;
% phiDotDot=pi^3*cos(pi*t_time/T)/T^2;
% num=4;   
% rx(jj,1)=aa*(cos(phi)+cos(num*phi)/num-1-1/num)+ix0;      %实际的海星形状极坐标
% ry(jj,1)=aa*(sin(phi)-sin(num*phi)/num)+iy0;
% rz(jj,1)=0+iz0;
% drx(jj,1)=-aa*((2*pi^2*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T))*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T + (2*pi^2*cos((pi*t_time)/(2*T))*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T)))/T);
% dry(jj,1)=aa*((2*pi^2*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T))*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T - (2*pi^2*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T)))/T);
% drz(jj,1)=0;
% ddrx(jj,1)=-aa*((pi^3*cos((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (pi^3*cos((pi*t_time)/(2*T))^2*sin(2*pi*num*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2 + (4*pi^4*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (4*pi^4*num*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2)/T^2);
% ddry(jj,1)=aa*((pi^3*cos((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))^2)/T^2 + (pi^3*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2 - (4*pi^4*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (4*pi^4*num*cos((pi*t_time)/(2*T))^2*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2);
% ddrz(jj,1)=0;

phi_sin=2*pi*sin(0.5*pi*t(jj)/T);
phi=phi_sin*sin(0.5*pi*t(jj)/T);
phiDot=phi_sin*pi*cos(0.5*pi*t(jj)/T)/T;
fac=exp(cos(phi))-2*cos(4*phi)-sin(phi/12)^5;%factor
rx(jj,1)=aa*sin(phi)*fac+ix0;%*T/2/pi;
ry(jj,1)=aa*(cos(phi)*fac-exp(1)+2)+iy0;%T/2/pi;-exp(1)+65/32+16*10^(-4)
rz(jj,1)=iz0;
drx(jj,1)=-aa*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*((5*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T + (2*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T) - (2*pi^2*aa*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T;
dry(jj,1)=(2*pi^2*aa*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T - aa*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*((5*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T + (2*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T);
drz(jj,1)=0;
ddrx(jj,1)=aa*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*((8*pi^3*cos((pi*t(jj))/(2*T))^2*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 - (8*pi^3*sin((pi*t(jj))/(2*T))^2*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (5*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5)/(36*T^2) - (5*pi^3*cos((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(12*T^2) - (pi^3*cos((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (5*pi^3*sin((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(12*T^2) + (pi^3*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (128*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*cos(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 - (5*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)^2*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^3)/(9*T^2) + (4*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2)^2)/T^2 - (4*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*cos(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2) - (pi^3*aa*cos((pi*t(jj))/(2*T))^2*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2 + (pi^3*aa*sin((pi*t(jj))/(2*T))^2*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2 - (4*pi^2*aa*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*((5*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T + (2*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T))/T + (4*pi^4*aa*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2;
ddry(jj,1)=aa*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*((8*pi^3*cos((pi*t(jj))/(2*T))^2*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 - (8*pi^3*sin((pi*t(jj))/(2*T))^2*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (5*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5)/(36*T^2) - (5*pi^3*cos((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(12*T^2) - (pi^3*cos((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (5*pi^3*sin((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(12*T^2) + (pi^3*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2 + (128*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*cos(8*pi*sin((pi*t(jj))/(2*T))^2))/T^2 - (5*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)^2*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^3)/(9*T^2) + (4*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2)^2)/T^2 - (4*pi^4*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*cos(2*pi*sin((pi*t(jj))/(2*T))^2))/T^2) + (pi^3*aa*cos((pi*t(jj))/(2*T))^2*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2 - (pi^3*aa*sin((pi*t(jj))/(2*T))^2*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2 + (4*pi^4*aa*cos((pi*t(jj))/(2*T))^2*sin((pi*t(jj))/(2*T))^2*cos(2*pi*sin((pi*t(jj))/(2*T))^2)*(sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^5 - exp(cos(2*pi*sin(1/2/T*pi*t(jj))^2)) + 2*cos(8*pi*sin((pi*t(jj))/(2*T))^2)))/T^2 + (4*pi^2*aa*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(2*pi*sin((pi*t(jj))/(2*T))^2)*((5*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*cos((pi*sin((pi*t(jj))/(2*T))^2)/6)*sin((pi*sin((pi*t(jj))/(2*T))^2)/6)^4)/(6*T) - (16*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*sin(8*pi*sin((pi*t(jj))/(2*T))^2))/T + (2*pi^2*cos((pi*t(jj))/(2*T))*sin((pi*t(jj))/(2*T))*exp(cos(2*pi*sin((pi*t(jj))/(2*T))^2))*sin(2*pi*sin((pi*t(jj))/(2*T))^2))/T))/T;
ddrz(jj,1)=0;

Q=[eye(n,n), J'; J, zeros(m,m)];
znn_d=[drx(jj,1);dry(jj,1);drz(jj,1)];
znn_u=[-k*(qjj-q0); znn_d];
error_function=Q*uAll(jj,:)'-znn_u;
error(jj)=norm(error_function,2);
t(jj)

end
%--Errors-
erposx=rx-posx;
erposy=ry-posy;
erposz=rz-posz;
averposx=mean(erposx);
averposy=mean(erposy);
averposz=mean(erposz);
erdposx=drx-dposx;
erdposy=dry-dposy;
erdposz=drz-dposz;
erddposx=ddrx-ddposx;
erddposy=ddry-ddposy;
erddposz=ddrz-ddposz;


save mvn_SRDdata1 t qAll dqAll ddqAll uAll posx posy posz dposx dposy dposz ddposx ddposy ddposz;
save mvn_SRDdata2 erposx erposy erposz erdposx erdposy  erdposz erddposx erddposy erddposz error averposx averposy averposz;
save mvn_SRDdata3 j0px j1px j2px j3px j4px j5px j0py j1py j2py j3py j4py j5py j0pz j1pz j2pz j3pz j4pz j5pz rx ry rz;

%以下为数据筛选和输出
%等间隔数据筛选
[m,n]=size(qAll);
step_long=0.000005;
jjj=1;
k=1;
q_step(1,:)=qAll(1,:);
for iii=2:m
    deta=(posx(iii)-posx(k))*(posx(iii)-posx(k))+(posy(iii)-posy(k))*(posy(iii)-posy(k))+(posz(iii)-posz(k))*(posz(iii)-posz(k));
    if deta>step_long||q_Time(iii)==1
        jjj=jjj+1;
        k=iii;
        q_step(jjj,:)=qAll(iii,:);
        t_step(jjj)=t(iii);
        q_deta(jjj,1)=deta;
        q_finish(jjj,1)=q_Time(iii);
    end
end
q_deta(1,1)=step_long;
%数据输出，由于Jaco建模和实际的Jaco有初始0度偏角以及旋转方向的的差异，在这里首先对偏角和旋转方向进行修正，再保存数据
q_tran=[2*pi*ones(jjj,1),(1/2*pi)*ones(jjj,1),(3/2*pi)*ones(jjj,1),0*ones(jjj,1),(pi)*ones(jjj,1),0*ones(jjj,1),0*ones(jjj,1),0*ones(jjj,1)];
q_Real=[-q_step(:,1),q_step(:,2),q_step(:,3),q_step(:,4),q_step(:,5),-q_step(:,6), q_deta, q_finish];
q_Real=q_Real+q_tran;

save mvn_Datarun q_Real;


