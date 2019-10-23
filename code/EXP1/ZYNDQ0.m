function dotDQ0=ZYNDQ0(t,DQ0)
global dr0 J0 mu;
dotDQ0=-mu*(J0'*J0*DQ0-J0'*dr0);