function gf=gfunction(rm,u,rp)
for i=1:length(u),
    if (u(i)>=rp(i)) 
        gf(i,1)=rp(i);
    elseif (u(i)<=rm(i)) 
        gf(i,1)=rm(i);
    else 
        gf(i,1)=u(i);
    end;
end;