clear;
format long;
load mvn_SYSdata;
load mvn_SRDdata1;
load mvn_SRDdata2;
load mvn_SRDdata3;

%==plotting angles==
figure;
plot(t,qAll(:,1),'m--',t,qAll(:,2),'g.-',t,qAll(:,3),'k.',t,qAll(:,4),'r:',t,qAll(:,5),'k--',t,qAll(:,6),'b-','linewidth',1.5);hold on;%axis([0 10 -2.5 2.5]);
% plot(t,0.0349,t,-1.7453);
title('Simulated theta');
xlabel('Time (Second)');
legend('q1','q2','q3','q4','q5','q6');
%%%%%%%%%%%%%%%%%%%%%%%%
% qplim1=qP(1)*ones(size(t));
% qplim2=qP(2)*ones(size(t));
% qplim3=qP(3)*ones(size(t));
% qplim4=qP(4)*ones(size(t));
% qplim5=qP(5)*ones(size(t));
% qplim6=qP(6)*ones(size(t));
% 
% qmlim1=qM(1)*ones(size(t));
% qmlim2=qM(2)*ones(size(t));
% qmlim3=qM(3)*ones(size(t));
% qmlim4=qM(4)*ones(size(t));
% qmlim5=qM(5)*ones(size(t));
% qmlim6=qM(6)*ones(size(t));
% %%%%%%%%%%%%%%%%%%%%
% 
% figure;
% plot(t,qAll(:,1),'linewidth',1.5);hold on;
% plot(t,qplim1,'k-',t,qmlim1,'k-','linewidth',1.5);
% ylim([-4 4]);
% xlim([0 15]);
% title('Simulated theta q1');
% xlabel('Time (Second)');
% 
% figure;
% plot(t,qAll(:,2),'linewidth',1.5);hold on;
% plot(t,qplim2,'k-',t,qmlim2,'k-','linewidth',1.5);
% ylim([-4 2]);
% xlim([0 15]);
% title('Simulated theta q2');
% xlabel('Time (Second)');
% 
% figure;
% plot(t,qAll(:,3),'linewidth',1.5);hold on;
% plot(t,qplim3,'k-',t,qmlim3,'k-','linewidth',1.5);
% ylim([-2 4]);
% xlim([0 15]);
% title('Simulated theta q3');
% xlabel('Time (Second)');
% 
% figure;
% plot(t,qAll(:,4),'linewidth',1.5);hold on;
% plot(t,qplim4,'k-',t,qmlim4,'k-','linewidth',1.5);
% ylim([-2.5 4]);
% xlim([0 15]);
% title('Simulated theta q4');
% xlabel('Time (Second)');
% 
% figure;
% plot(t,qAll(:,5),'linewidth',1.5);hold on;
% plot(t,qplim5,'k-',t,qmlim5,'k-','linewidth',1.5);
% ylim([-2 0.5]);
% xlim([0 15]);
% title('Simulated theta q5');
% xlabel('Time (Second)');
% 
% figure;
% plot(t,qAll(:,6),'linewidth',1.5);hold on;
% plot(t,qplim6,'k-',t,qmlim6,'k-','linewidth',1.5);
% ylim([-4 4]);
% xlim([0 15]);
% title('Simulated theta q6');
% xlabel('Time (Second)');
% 
% 
% 
% 
% 
% 
% 
% 
figure;
plot(t,dqAll(:,1),'m--',t,dqAll(:,2),'g.-',t,dqAll(:,3),'k.',t,dqAll(:,4),'r:',t,dqAll(:,5),'k--',t,dqAll(:,6),'b-','linewidth',1.5);hold on;
title('Simulated dotTheta');
xlabel('Time (Second)');
legend('dq1','dq2','dq3','dq4','dq5','dq6');

figure;
plot(t,ddqAll(:,1),'m--',t,ddqAll(:,2),'g.-',t,ddqAll(:,3),'k.',t,ddqAll(:,4),'r:',t,ddqAll(:,5),'k--',t,ddqAll(:,6),'b-','linewidth',1.5);hold on;
title('Simulated dotdotTheta');
xlabel('Time (Second)');
legend('ddq1','ddq2','ddq3','ddq4','ddq5','ddq6');

% figure
% plot(t,uAll(:,1:7));hold on;
% xlabel('Time (Second)');
% title('compare last figure: the same or not?');

% figure;
% plot(t,tauAll);hold on;
% title('Simulated Torque');
% xlabel('Time (Second)');
% legend('\tau_1','\tau_2','\tau_3','\tau_4','\tau_5','\tau_6');

%==plotting 3d trajectories of joints==
figure;
grid on;
plot3(j2px,j2py,j2pz);hold on;
plot3(j3px,j3py,j3pz);hold on;
plot3(j4px,j4py,j4pz);hold on;
plot3(j5px,j5py,j5pz);hold on;
plot3(posx,posy,posz);hold on;
xlabel('Time (Second)');
title('joints Position');
legend('X','Y','Z');

%==plotting end effector==
figure;
plot(t,posx,'m--',t,posy,'c.-',t,posz,'k-');hold on;
%axis([0 10 -0.5 1]);
xlabel('Time (Second)');
title('End Effector Position');
legend('X','Y','Z');

figure;
plot(t,dposx,t,dposy,t,dposz);
%axis([0 10 -0.2 0.2]);
xlabel('Time (Second)');
title('End Effector Velocity');
legend('dX','dY','dZ');

figure;
plot(t,ddposx,t,ddposy,t,ddposz);
%axis([0 10 -0.2 0.2]);
xlabel('Time (Second)');
title('End Effector Acceleration');
legend('ddX','ddY','ddZ');

figure;
plot(t,erposx,'m--',t,erposy,'b-',t,erposz,'k-.','linewidth',2.5);
title('Position Error');
xlabel('Time (Second)');
legend('ex','ey','ez');

figure;

plot(t,erdposx,t,erdposy,t,erdposz);
title('Velocity Error');
xlabel('Time (Second)');
legend('edx','edy','edz');

figure;
plot(t,erddposx,t,erddposy,t,erddposz);
%axis([0 10 -15*10^-5 15*10^-5]);
title('Acceleration Error');
legend('eddx','eddy','eddz');
xlabel('Time (Second)');

figure;
plot3(rx,ry,rz,'k--','linewidth',2);
hold on;
for jj=1:8:length(t)
plot3(posx(jj),posy(jj),posz(jj),'md-','linewidth',2);
grid on;
hold on;
end
legend('Expected Path','Actual Trajectories');
xlabel('x');
ylabel('y');
zlabel('z');


u_initial=u0;
[Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qAll(1,:));
pos_initial=[Ppx Ppy Ppz]
[Ppx,Ppy,Ppz]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qAll(length(t),:));
pos_final=[Ppx Ppy Ppz]
EndeffectorPe=pos_final(6,:)-pos_initial(6,:);
%save position status information to "txt" file
fid=fopen('joint_position_information.txt','a');
fprintf(fid,'u_initial: %g\n',u_initial);
% fprintf(fid,'initial_joint_position: %g\n',pos_initial);
% fprintf(fid,'final_joint_position: %g\n',pos_final);
fprintf(fid,'End-effector_position_error: %g\n',EndeffectorPe);
fclose(fid);
%save join drift information to "txt" file
joindrift=qAll(length(t),:)-qAll(1,:)
joindriftnorm=norm(joindrift)
fid=fopen('join_angle_drift_information.txt','a');
fprintf(fid,'u_initial: %g\n',u_initial);
fprintf(fid,'joindrift: %g\n',joindrift);
fprintf(fid,'joindriftnorm: %g\n',joindriftnorm);
fclose(fid);

% figure;
% plot(t,qAll(:,1),'m--','linewidth',1.5);hold on;
% legend('q1');
% figure;
% plot(t,qAll(:,2),'g.-','linewidth',1.5);hold on;
% legend('q2');
% figure;
% plot(t,qAll(:,3),'k.','linewidth',1.5);hold on;
% legend('q3');
% figure;
% plot(t,qAll(:,4),'r:','linewidth',1.5);hold on;
% legend('q4');
% figure;
% plot(t,qAll(:,5),'k--','linewidth',1.5);hold on;
% legend('q5');
% figure;
% plot(t,qAll(:,6),'b-','linewidth',1.5);hold on;
% legend('q6');

%%????????????
figure(30);  
plot(t,error,'linewidth',2.5);
title('function error')
xlabel('Time(Second)');
ylabel('error');
hold on;