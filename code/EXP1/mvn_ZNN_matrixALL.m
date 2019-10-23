function znn_ALL=mvn_ZNN_matrixALL(t,y)

global n m alpha_i a_i d_i D sa ca s2a c2a pen_long  

q=y(1:n);

J=ZZJLJrebiJacobian(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q);  
znn_W=eye(n,n);
znn_Q=[znn_W, J'; J, zeros(m,m)];
znn_ALL=[znn_W,zeros(n,n+m); zeros(n+m,n),znn_Q];