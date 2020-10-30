
L=[0 1 0; 0 0 0;0 0 0];
tincr=.05;
steps=500;
SI=0;
tic
for i=1:steps
    
    
    [SI,Gamma,~]=Calc_Gamma(L,tincr,i);
    if SI>=2
        break
    end

end
toc