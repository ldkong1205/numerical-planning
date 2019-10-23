function ZYNplot(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,tg)
    %VERSION of April 15th, 2001
	np=numrows(tg);
	n=6;
    erasemode='xor';
	wrist=1;
	base=0.0;
	reach=0.1*(D(1)+D(2)+D(3)+D(4)+D(5)+D(6));
    figure(gcf);%bring to the top
    grid on
	title('PUMA560 with d6');
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
	set(gca,'drawmode','fast');%get handle of current axes, but ?
	line('xdata', [0;0], 'ydata', [0;0], 'zdata', [-reach;base], 'color', 'magenta');
	%create a line which we will subsequently modify. Set erase mode to xor for fast update
	hr=line('color', 'black', 'erasemode', erasemode,'LineWidth',2);
	hx=line('xdata', [0;0], 'ydata', [0;0], 'zdata', [0;base], ...
		'color', 'red', 'erasemode', 'xor');
	hy=line('xdata', [0;0], 'ydata', [0;0], 'zdata', [0;base], ...
		'color', 'green', 'erasemode', 'xor');
	hz=line('xdata', [0;0], 'ydata', [0;0], 'zdata', [0;base], ...
		'color', 'blue', 'erasemode', 'xor');
	mag=0.3*reach;
    %save the handles in the passed robot object
	handle=[hr hx hy hz mag];%handle
	%assignin('base',inputname(1),p560);
    %i still don't know what the above sentence means?
	for repeat=1:1,
    for p=1:np,
        ZYNanimate(handle,alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,tg(p,:));
        for ii=1:5000,iii=ii^2;end;
	end
    end