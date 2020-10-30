function [SI,Gamma,f]=Calc_Gamma(L,tincr,steps)
  % calculating shear strain gamma for the number of steps of simple shear
  % flow.            

   t=tincr*steps;
   f =expm(L*t);
   [emat,~]=eig(f*f');
   eiganval=eig(f*f'); % this gives the principal quadric stretches of strain ellipsoid
   [~,ix]=sort(eiganval,'descend');
   ls=emat(:,ix(1));
   Y = [0;1;0];
   %X=[1;0;0];
   phi= acos(dot(ls,Y));
   %phi1= acos(dot(ls,X));
   Gamma= 2/tan(2*phi);
   %Gamma1= 2/tan(2*phi1);
    % strain_Intensity
    SI= ((log(eiganval(3)/eiganval(2)))^2 + (log(eiganval(2)/eiganval(1)))^2)^.5 ;
end