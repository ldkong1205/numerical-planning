
clear all;

t=0:0.1:10;
x=-t;
plot(t,x,'go-');
hold on;
plot(t,x+1,'ks-.');
legend 'velocity-level resolution' 'acceleration-level resolution';

