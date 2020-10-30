function plot_rou_vs_t(Wk)

%  Input parameters--------------------------------------------------------
  %Wk= 1;
  gamma=1;
  epsilon= .5*(((1/(Wk^2))-1)^(.5));
 if Wk==0
     gamma=0;
     epsilon=1;
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 else
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 end 
 D = 0.5 * (L + L');
% normalizing with respect to strain rate
 L= L/norm(D);

%  Matrix flow field
   %L     = [0 1 0;0 0 0;0 0 0]; 
   
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 200;
  
% ------------------------------------------------------------------------   

SI=zeros(1,steps);
Gamma=zeros(1,steps);
  
for i=1:steps
   
 [SI(i),Gamma(i),~]=Calc_Gamma(L,tincr,i);
 %if SI>=6.5
  %   break
 %end
end

STEP= 1:steps;
if Wk==1
    plot(STEP,Gamma,'*');
    pbaspect([1 1 1])
    hold on
    plot(STEP,SI);
    pbaspect([1 1 1])
    hold on
else
    plot(STEP,SI)
    pbaspect([1 1 1])
    hold on
end

end

