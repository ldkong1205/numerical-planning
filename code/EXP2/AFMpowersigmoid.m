function y=AFMpowersigmoid(E, xi, p)
if nargin==1, xi=1; p=1;
elseif nargin==2, p=3;
end
y=(1+exp(-xi))/(1-exp(-xi))*(1-exp(-xi*E))./(1+exp(-xi*E));
i=find(abs(E)>=1);
y(i)=E(i).^p;