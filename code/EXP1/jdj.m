function [J,DJ]=jdj(q,dq) 
%ian planar manipulator 2x6
% J is the 2x6 jocobian matrix of a 6dof planar, unit-link-length robot.
%q/theta is a 6x1 column vector of angles in radians meansuring ccw from
%the previous link (from horizontal for the first link).

len=[0.2755, 0.4100, 0.2073, 0.07433, 0.07433, 0.13];
len1=len(1);
len2=len(2);
len3=len(3);
len4=len(4);
len5=len(5);
len6=len(6);

s6=len6*sin(q(1)+q(2)+q(3)+q(4)+q(5)+q(6));
s5=len5*sin(q(1)+q(2)+q(3)+q(4)+q(5));
s4=len4*sin(q(1)+q(2)+q(3)+q(4));
s3=len3*sin(q(1)+q(2)+q(3));
s2=len2*sin(q(1)+q(2));
s1=len1*sin(q(1));

c6=len6*cos(q(1)+q(2)+q(3)+q(4)+q(5)+q(6));
c5=len5*cos(q(1)+q(2)+q(3)+q(4)+q(5));
c4=len4*cos(q(1)+q(2)+q(3)+q(4));
c3=len3*cos(q(1)+q(2)+q(3));
c2=len2*cos(q(1)+q(2));
c1=len1*cos(q(1));

d6=dq(1)+dq(2)+dq(3)+dq(4)+dq(5)+dq(6);
d5=dq(1)+dq(2)+dq(3)+dq(4)+dq(5);
d4=dq(1)+dq(2)+dq(3)+dq(4);
d3=dq(1)+dq(2)+dq(3);
d2=dq(1)+dq(2);
d1=dq(1);
F=[c1+c2+c3+c4+c5+c6,c2+c3+c4+c5+c6,c3+c4+c5+c6,c4+c5+c6,c5+c6,c6;...
    s1+s2+s3+s4+s5+s6,s2+s3+s4+s5+s6,s3+s4+s5+s6,s4+s5+s6,s5+s6,s6];
J=[-s1-s2-s3-s4-s5-s6,-s2-s3-s4-s5-s6,-s3-s4-s5-s6,-s4-s5-s6,-s5-s6,-s6;...
    c1+c2+c3+c4+c5+c6,c2+c3+c4+c5+c6,c3+c4+c5+c6,c4+c5+c6,c5+c6,c6];
DJ=[-c1*d1-c2*d2-c3*d3-c4*d4-c5*d5-c6*d6,-c2*d2-c3*d3-c4*d4-c5*d5-c6*d6,-c3*d3-c4*d4-c5*d5-c6*d6,-c4*d4-c5*d5-c6*d6,-c5*d5-c6*d6,-c6*d6;...
    -s1*d1-s2*d2-s3*d3-s4*d4-s5*d5-s6*d6,-s2*d2-s3*d3-s4*d4-s5*d5-s6*d6,-s3*d3-s4*d4-s5*d5-s6*d6,-s4*d4-s5*d5-s6*d6,-s5*d5-s6*d6,-s6*d6];