function du=ZYNU0(t,u)
global mu J0 DJ0 taug0 tauc0 H0 ddr0 dq0;
u=[u(1);u(2);u(3)];
ddq=H0*(J0'*u-tauc0-taug0);
du=(ddr0-DJ0*dq0-J0*ddq)*mu;