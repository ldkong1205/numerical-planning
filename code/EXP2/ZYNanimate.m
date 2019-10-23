function ZYNanimate(handle,alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q)
	n=6;
	h=handle;
    hr=h(1);
    hx=h(2);hy=h(3);hz=h(4);mag=h(5);
    [x,y,z,t]=ZYNposition(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,q);
    %compute the wrist axes
    xv=t*[mag;0;0;1];
    yv=t*[0;mag;0;1];
    zv=t*[0;0;mag;1];
    %update the line segments, wrist axis and links
    set(hx,'xdata',[t(1,4) xv(1)], 'ydata', [t(2,4) xv(2)], ...
        'zdata', [t(3,4) xv(3)]);
    set(hy,'xdata',[t(1,4) yv(1)], 'ydata', [t(2,4) yv(2)], ...
        'zdata', [t(3,4) yv(3)]);
    set(hz,'xdata',[t(1,4) zv(1)], 'ydata', [t(2,4) zv(2)], ...
        'zdata', [t(3,4) zv(3)]);
    set(hr,'xdata',x,'ydata',y,'zdata', z);
	drawnow