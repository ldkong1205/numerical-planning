%old message:
%1.tau_gravity=ZYNrne(q,zeros(1,6),zeros(1,6),[0;0;9.81])
%2.tau_coriolis=ZYNrne(q,dq,zeros(1,6),[0;0;0])
%3.inertia=[ZYNrne(q,zeros(1,6),[1,0,0,0,0,0],[0;0;0]);...
%        ZYNrne(q,zeros(1,6),[0,1,0,0,0,0],[0;0;0]);...
%        ZYNrne(q,zeros(1,6),[0,0,1,0,0,0],[0;0;0]);...
%        ZYNrne(q,zeros(1,6),[0,0,0,1,0,0],[0;0;0]);...
%        ZYNrne(q,zeros(1,6),[0,0,0,0,1,0],[0;0;0]);...
%        ZYNrne(q,zeros(1,6),[0,0,0,0,0,1],[0;0;0])];
%new usage below (VERSION of April 16th 2001):
%[taug,tauc,inerH]=ZYNtaugcI(a2,a3,d3,d4,d6,q,dq)
function [taug,tauc,inerH]=ZYNtaugcI(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q,qd)%new add
    %common data below
    keepqd=qd;
    z0=[0;0;1];fext=zeros(6,1);%force/moments on end of arm 
    alpha_r=[1.5707963267949,0,-1.5707963267949,1.5707963267949,-1.5707963267949,0];
    linkR=[zeros(1,3);-0.3638,0.006,0.2275;-0.0203,-0.0141,0.070;0,0.019,0;zeros(1,3);0,0,0.032];
    linkM=[0;17.4;4.8;0.82;0.34;0.09];
    linkI=[0,0.35,0,0,0,0;0.13,0.524,0.539,0,0,0;0.066,0.086,0.0125,0,0,0;...
            1.8e-3,1.3e-3,1.8e-3,0,0,0;0.3e-3,0.4e-3,0.3e-3,0,0,0;0.15e-3,0.15e-3,0.04e-3,0,0,0];
    linkJm=[200e-6;200e-6;200e-6;33e-6;33e-6;33e-6];
    linkG=[-62.6111;107.815;-53.7063;76.0364;71.923;76.686];
    mys=sin(q);
    s1=mys(1);s2=mys(2);s3=mys(3);s4=mys(4);s5=mys(5);s6=mys(6);
    myc=cos(q);
    c1=myc(1);c2=myc(2);c3=myc(3);c4=myc(4);c5=myc(5);c6=myc(6);
    
    alpha_4 = alpha_i(4);alpha_5 = alpha_i(5);
    d_1 = d_i(1);d_2 = d_i(2);d_3 = d_i(3);d_4 = d_i(4);
    d_5 = d_i(5);d_6 = d_i(6);
%     sp = sin(phi);
%     cp = cos(phi);
    sa4 = sin(alpha_4);ca4 = cos(alpha_4);
    sa5 = sin(alpha_5);ca5 = cos(alpha_5);
 
    sp1=sin(0.02);cp1=cos(0.02);
    sp2=sin(0.05);cp2=cos(0.05);    
    
    AW1 = [1, 0, 0, 0;...
    0, cp1, sp1, 0;...
    0, -sp1, cp1, 0;...
    0, 0, 0, 1];

    AW2 = [cp2, 0, -sp2, 0;...
    0, 1, 0, 0;...
    sp2, 0, cp2, 0;...
    0, 0, 0, 1];

    AW1=AW1*AW2;
    
    A1 = [c1, 0, s1, 0;s1, 0, -c1, 0;0, 1, 0, d_1;0, 0, 0, 1];
    A2 = [c2, s2, 0, D(2)*c2;s2, -c2, 0, D(2)*s2;0, 0, -1, 0;0, 0, 0, 1];
    A3 = [c3, 0, s3, 0;s3, 0, -c3, 0;0, 1, 0, d_3;0, 0, 0, 1];
    A4 = [c4, -ca4*s4, sa4*s4, 0;s4, ca4*c4, -sa4*c4, 0;0, sa4, ca4, d_4;0, 0, 0, 1];
    A5 = [c5, -ca5*s5, sa5*s5, 0;s5, ca5*c5, -sa5*c5, 0;0, sa5, ca5, d_5;0, 0, 0, 1];
    A6 = [c6, s6, 0, 0;s6, -c6, 0, 0;0, 0, -1, d_6;0, 0, 0, 1];%new add
    Rm{1}=[A1(1:3,1:3)];
    Rm{2}=[A2(1:3,1:3)];
    Rm{3}=[A3(1:3,1:3)];
    Rm{4}=[A4(1:3,1:3)];
    Rm{5}=[A5(1:3,1:3)];
    Rm{6}=[A6(1:3,1:3)];
    pstarm=[0,0,0;0,0,0;0,0*sin(alpha_r(3)),0*cos(alpha_r(3));...
            0,0*sin(alpha_r(4)),0*cos(alpha_r(4));0,0,0;...
            0,0*sin(alpha_r(6)),0*cos(alpha_r(6))]';%new add
    
    %calculate inertial matrix H
    for ji=1:6,
        grav=[0;0;0];qd=zeros(1,6);
        qdd=zeros(1,6);qdd(ji)=1;
        w=zeros(3,1);wd=zeros(3,1);vd=grav;v=zeros(3,1);%seemingly v useless
    	for j=1:6,%the forward recursion
    		R=Rm{j}';
    		pstar=pstarm(:,j);
    		r=linkR(j,:)';
    		wd=R*(wd+z0*qdd(j)+cross(w,z0*qd(j)));
    		w=R*(w+z0*qd(j));
            vd=cross(wd,pstar)+cross(w,cross(w,pstar))+R*vd;
            vhat=cross(wd,r)+cross(w,cross(w,r))+vd;
            F=linkM(j)*vhat;
            LinkIj=diag(linkI(j,1:3));
    		N=LinkIj*wd+cross(w,LinkIj*w);
            if j==1
                Fm=F;
                Nm=N;
            else
                Fm=[Fm F];
                Nm=[Nm N];
            end;		
    	end
    	f=fext(1:3);nn=fext(4:6);
    	for j=6:-1:1,%the backward recursion
            pstar=pstarm(:,j);
    		if j==6,
    			R=eye(3,3);
    		else
    			R=Rm{j+1};
    		end
            r=linkR(j,:)';
    		nn=R*(nn+cross(R'*pstar,f))+cross(pstar+r,Fm(:,j))+Nm(:,j);
    		f=R*f+Fm(:,j);
    		R=Rm{j};
    		inerH(ji,j)=nn'*(R'*z0)+linkJm(j)*qdd(j)*(linkG(j)^2)+0;
        end
    end
    
    %calculate taug
    grav=[0;0;9.81];qd=zeros(1,6);qdd=zeros(1,6);
    w=zeros(3,1);wd=zeros(3,1);vd=grav;v=zeros(3,1);%seemingly v useless
	for j=1:6,%the forward recursion
		R=Rm{j}';
		pstar=pstarm(:,j);
		r=linkR(j,:)';
		wd=R*(wd+z0*qdd(j)+cross(w,z0*qd(j)));
		w=R*(w+z0*qd(j));
        vd=cross(wd,pstar)+cross(w,cross(w,pstar))+R*vd;
        vhat=cross(wd,r)+cross(w,cross(w,r))+vd;
        F=linkM(j)*vhat;
        LinkIj=diag(linkI(j,1:3));
		N=LinkIj*wd+cross(w,LinkIj*w);
        if j==1
            Fm=F;
            Nm=N;
        else
            Fm=[Fm F];
            Nm=[Nm N];
        end;		
	end
	f=fext(1:3);nn=fext(4:6);%the backward recursion
	for j=6:-1:1,
        pstar=pstarm(:,j);
		if j==6,
			R=eye(3,3);
		else
			R=Rm{j+1};
		end
        r=linkR(j,:)';
		nn=R*(nn+cross(R'*pstar,f))+cross(pstar+r,Fm(:,j))+Nm(:,j);
		f=R*f+Fm(:,j);
		R=Rm{j};
		taug(1,j)=nn'*(R'*z0)+linkJm(j)*qdd(j)*(linkG(j)^2)+0;
    end
    
    %calculate tauc
    grav=[0;0;0];qd=keepqd;qdd=zeros(1,6);
    w=zeros(3,1);wd=zeros(3,1);vd=grav;v=zeros(3,1);%seemingly v useless
	for j=1:6,%the forward recursion
		R=Rm{j}';
		pstar=pstarm(:,j);
		r=linkR(j,:)';
		wd=R*(wd+z0*qdd(j)+cross(w,z0*qd(j)));
		w=R*(w+z0*qd(j));
        vd=cross(wd,pstar)+cross(w,cross(w,pstar))+R*vd;
        vhat=cross(wd,r)+cross(w,cross(w,r))+vd;
        F=linkM(j)*vhat;
        LinkIj=diag(linkI(j,1:3));
		N=LinkIj*wd+cross(w,LinkIj*w);
        if j==1
            Fm=F;
            Nm=N;
        else
            Fm=[Fm F];
            Nm=[Nm N];
        end;		
	end
	f=fext(1:3);nn=fext(4:6);%the backward recursion
	for j=6:-1:1,
        pstar=pstarm(:,j);
		if j==6,
			R=eye(3,3);
		else
			R=Rm{j+1};
		end
        r=linkR(j,:)';
		nn=R*(nn+cross(R'*pstar,f))+cross(pstar+r,Fm(:,j))+Nm(:,j);
		f=R*f+Fm(:,j);
		R=Rm{j};
		tauc(1,j)=nn'*(R'*z0)+linkJm(j)*qdd(j)*(linkG(j)^2)+0;
    end