clear;
format long;
load mvn_SYSdata;
load mvn_INITdata;
 global aa q0 T mu_p mu_v qP qM qDp qDm qDDp qDDm m Pinfty Minfty gamma myInf n tt errorr;

interval=0.005;                         %????????
jj=0;
epsilon=0;
tt=0; errorr=0;
for ii=1:length(t),
    if(t(ii,1)>=epsilon)                %??????????????????????????????????????
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        yn(jj,:)=y(ii,:);
        %%%%%%%%%%%%%%%%%
        doty=mvn_ZZJLJrebi1_net(tn(jj,1),yn(jj,:)');     %????????????????????????????????
        ddqAll(jj,:)=doty((n+1):(n+n),1)';              %????u_dot????n??????
        %%%%%%%%%%%%%%%%%
        epsilon=jj*interval;        
    elseif(ii==length(t))
        jj=jj+1;
        tn(jj,1)=t(ii,1);
        yn(jj,:)=y(ii,:);
        %%%%%%%%%%%%%%%%%
        doty=mvn_ZZJLJrebi1_net(tn(jj,1),yn(jj,:));
        ddqAll(jj,:)=doty((n+1):(n+n),1)';
        %%%%%%%%%%%%%%%%%
        epsilon=jj*interval;
    end
end
clear t y;
t=tn;
y=yn;
size(t)
size(y)
clear tn yn;
qAll=y(:,1:n);                                      %??????????Theta
size(qAll)
uAll=y(:,(n+1):(n+n+m));                            %??????????u,????????????????
size(uAll)
dqAll=uAll(:,1:n);                                  %????u????n??
size(dqAll)
size(ddqAll)
save mvn_SNDdata t qAll dqAll uAll ddqAll;

% ????????
% 
% [m,n]=size(qAll);
% q_tran=[2*pi*ones(m,1),1/2*pi*ones(m,1),3/2*pi*ones(m,1),0*ones(m,1),pi*ones(m,1),0*ones(m,1),];
% q_Real=[-qAll(:,1),qAll(:,2),qAll(:,3),qAll(:,4),qAll(:,5),-qAll(:,6)];
% q_Real=q_Real+q_tran;
% 
% save mvn_Datarun q_Real t;