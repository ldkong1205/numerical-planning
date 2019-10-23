function doty=mvn_ZZJLJrebi1_net(t,y)

global aa alpha_i a_i d_i D sa ca s2a c2a pen_long q0 T m gamma nu n Time znn_W k tt errorr K;

persistent ix0 iy0 iz0 ix1 iy1 iz1 t_last;


   q=y(1:n);                         %??????????
   dq=y(n+1:n+n);                        %????Theta dot??????????????Theta dot
   yy=y(n+1:n+n+m);


if t==0
    [x,y,z]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q0);  %????????????
        
    ix0=x(7);iy0=y(7);iz0=z(7);
end

t_time=t-Time*T;                    %??????????????????????????????????????????????????????????????????????????????????

if t_time>T;
   Time=Time+1;
   t_time=t_time-T;
end



num=4;                                  %4??????4+1=5????
phi_sin=2*pi*sin(0.5*pi*t_time/T);
phi=phi_sin*sin(0.5*pi*t_time/T);
phiDot=phi_sin*pi*cos(0.5*pi*t_time/T)/T;

%??????????????1????????r
rx=aa*(cos(phi)+cos(num*phi)/num-1-1/num);      %????????????????????
ry=aa*(sin(phi)-sin(num*phi)/num);
rz=0;
%??????????????2????????r_dot
drx=-aa*((2*pi^2*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T))*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T + (2*pi^2*cos((pi*t_time)/(2*T))*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T)))/T);
dry=aa*((2*pi^2*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T))*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T - (2*pi^2*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))*sin((pi*t_time)/(2*T)))/T);
drz=0;
%??????????????????r_dotdot
ddrx=-aa*((pi^3*cos((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (pi^3*cos((pi*t_time)/(2*T))^2*sin(2*pi*num*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2 + (4*pi^4*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (4*pi^4*num*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2)/T^2);
ddry=aa*((pi^3*cos((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*sin((pi*t_time)/(2*T))^2*cos(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 - (pi^3*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*cos((pi*t_time)/(2*T))^2)/T^2 + (pi^3*cos(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2 - (4*pi^4*cos((pi*t_time)/(2*T))^2*sin((pi*t_time)/(2*T))^2*sin(2*pi*sin((pi*t_time)/(2*T))^2))/T^2 + (4*pi^4*num*cos((pi*t_time)/(2*T))^2*sin(2*pi*num*sin((pi*t_time)/(2*T))^2)*sin((pi*t_time)/(2*T))^2)/T^2);
ddrz=0;

dr=[drx;dry;drz];           %????????????????r_dot
ddr=[ddrx;ddry;ddrz];
znn_d=dr;
[J,DJ]=ZYNjdj(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q,dq); 


znn_W=eye(n,n);
znn_Q=[znn_W, J'; J, zeros(m,m)];
znn_dotQ=[zeros(n,n), DJ'; DJ, zeros(m,m)];


% C=[J,zeros(m,1)];
k=6;                                                %????????????c????????gamma????????0????????????????????????
% zz=k*(q-q0);
znn_q=k*(q-q0);                                         %????????????c
znn_u=[-znn_q; znn_d];
znn_dotu=[zeros(n,1); ddr];


%%%%%%%%%%%%%%%%%
%VP-CDNN??????????????
% dotyy=-znn_dotQ*yy-gamma*exp(t)*AFMlinear(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*exp(t)*AFpower(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*exp(t)*AFsigmoid(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*exp(t)*AFsinh(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*exp(t)*AFMpowersigmoid(znn_Q*yy-znn_u)+znn_dotu;

%%%%%%%%%%%%%%%%%
%FP-CDNN??????????????
% dotyy=-znn_dotQ*yy-gamma*AFMlinear(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*AFpower(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*AFsigmoid(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*AFsinh(znn_Q*yy-znn_u)+znn_dotu;
% dotyy=-znn_dotQ*yy-gamma*AFMpowersigmoid(znn_Q*yy-znn_u)+znn_dotu;

%%%%%%%%%%%%%%%%%????
delt_D=1.0*[sin(t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(3*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             cos(2*t/T*pi) cos(2*t/T*pi) -cos(4*t/T*pi) sin(t/T*pi) -sin(2*t/T*pi) -cos(3*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(3*t/T*pi) -sin(4*t/T*pi) cos(5*t/T*pi) -cos(3*t/T*pi) cos(4*t/T*pi) -sin((5*t/T*pi)) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(2*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi);...
             sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(t/T*pi);...
             -sin(2*t/T*pi) cos(2*t/T*pi) -sin(t/T*pi) cos(t/T*pi) -sin(t/T*pi) cos(2*t/T*pi) -sin(2*t/T*pi) cos(t/T*pi) -sin(t/T*pi)];
delt_s=1.5*[2*sin(t/T*pi);4*cos(2*t/T*pi);-5*sin(3*t/T*pi); cos(3*t/T*pi); -sin(3*t/T*pi); 3*cos(t/T*pi);-sin(t/T*pi);-sin(2*t/T*pi);cos(t/T*pi)];

% delt_D=0.01*rand(n+m,n+m);
% delt_s=t*ones(9,1);
% delt_s=[1;2;3;4;5;6;7;8;9];
% delt_s=0.1*rand(n+m,1);



%IE-RNN
dotyy=-(znn_dotQ)*yy-gamma*AFMpowersigmoid(znn_Q*yy-znn_u)+znn_dotu-nu*(znn_Q*yy-znn_u)*t+delt_s;

%Z-RNN
%dotyy=-(znn_dotQ)*yy-gamma*AFMpowersigmoid(znn_Q*yy-znn_u)+znn_dotu+delt_s;


t

doty=[dq;dotyy];
