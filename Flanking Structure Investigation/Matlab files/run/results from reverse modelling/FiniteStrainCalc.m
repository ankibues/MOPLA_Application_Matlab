
% calculating finite strain from the reverse modelling data 

load('Result2_nature2_10p.mat')
 K=zeros(1,1000);
for i=1:1000
   K(i)= min(real(ic(i).residual));
end

 [val,ind]=min(K);   
   ic(ind).WK;
   Wk= ic(ind).WK;    % Wk
   
  gamma=1;
  epsilon= .5*(((1/(Wk^2))-1)^(.5));
 if Wk==0
     gamma=0;
     epsilon=1;
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 else
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 end 

steps= 22;
tincr=.05;
[SI,~,~]=Calc_Gamma(L,tincr,steps); % strain intensity

viscosityratio= ic(ind).viscR; % viscosity ratio