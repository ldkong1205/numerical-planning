clear;
format long;
load mvn_SYSdata;
load mvn_SRDdata1;
load mvn_SRDdata2;
load mvn_SRDdata3;

figure;
plot3(posx,posy,posz,'g');hold on;
for j=1:15:length(t),
    linksX=[j0px(j);j1px(j);j2px(j);j3px(j);j4px(j);j5px(j);posx(j)];
    linksY=[j0py(j);j1py(j);j2py(j);j3py(j);j4py(j);j5py(j);posy(j)];
    linksZ=[j0pz(j);j1pz(j);j2pz(j);j3pz(j);j4pz(j);j5pz(j);posz(j)];
    plot3(linksX,linksY,linksZ,'k');hold on;
end;
%title('Black-white PUMA560 with d6');
%xlabel('X');ylabel('Y');zlabel('Z');
text(0.25,-0.6,-0.2,'X (m)');
text(-0.1,0,-0.2,'Y (m)');
text(0,0.4,0.9,'Z (m)');
grid on;

figure;
plot3(rx,ry,rz,'k--','linewidth',2);
hold on;
for jj=1:8:length(t)
plot3(posx(jj),posy(jj),posz(jj),'md-','linewidth',2);
grid on;
hold on;
end
legend('Expected Path','Actual Trajectories');
xlabel('x');
ylabel('y');
zlabel('z');
plot3(posx,posy,posz, 'LineWidth',2);hold on;
for j=1:20:length(t),
    %if(rem(j,10)==0)
    basejoint1=line('xdata',[j0px(j);j1px(j)],'ydata',[j0py(j);j1py(j)],'zdata',[j0pz(j);j1pz(j)],...
        'color', 'yellow','erasemode','none');
    joints12=line('xdata',[j1px(j);j2px(j)],'ydata',[j1py(j);j2py(j)],'zdata',[j1pz(j);j2pz(j)],...
        'color', 'red','erasemode','none');
    joints23=line('xdata',[j2px(j);j3px(j)],'ydata',[j2py(j);j3py(j)],'zdata',[j2pz(j);j3pz(j)],...
        'color', 'cyan','erasemode','none');
    joints34=line('xdata',[j3px(j);j4px(j)],'ydata',[j3py(j);j4py(j)],'zdata',[j3pz(j);j4pz(j)],...
        'color', 'magenta','erasemode','none');
    joints45=line('xdata',[j4px(j);j5px(j)],'ydata',[j4py(j);j5py(j)],'zdata',[j4pz(j);j5pz(j)],...
        'color', 'blue','erasemode','none');
    joints56=line('xdata',[j5px(j);posx(j)],'ydata',[j5py(j);posy(j)],'zdata',[j5pz(j);posz(j)],...
        'color', 'green','erasemode','none');
    drawnow
    %end
end;
%title('colourful PUMA560 with d6');
xlabel('X');ylabel('Y');zlabel('Z');
grid on;

figure;
linksX=[j0px(j);j1px(j);j2px(1);j3px(1);j4px(1);j5px(1);posx(1)];
linksY=[j0py(j);j1py(j);j2py(1);j3py(1);j4py(1);j5py(1);posy(1)];
linksZ=[j0pz(j);j1pz(j);j2pz(1);j3pz(1);j4pz(1);j5pz(1);posz(1)];
plot3(linksX,linksY,linksZ,'r','linewidth',2);hold on;
plot3(posx,posy,posz,':','LineWidth',2);hold on;
ZYNplot(alpha_i,a_i,d_i,D,sa,ca,s2a,c2a,pen_long,qAll);
hold on;