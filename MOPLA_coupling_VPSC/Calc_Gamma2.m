function [SI,Gamma,f]=Calc_Gamma2(L,tincr,steps)
  % calculating shear strain gamma for the number of steps of simple shear
  % flow.            

   t=tincr*steps;
   f =expm(L*t);
   [emat,~]=eig(f*f');
   eiganval=sqrt(eig(f*f')); % this gives the principal stretches of strain ellipsoid
   [~,ix]=sort(eiganval,'descend');
   ls=emat(:,ix(1));
   Y= [0;1;0];
   phi= -acos(dot(ls,Y));
   Gamma= -2/tan(2*phi);
 
    % strain_Intensity
    SI= ((log(eiganval(3)/eiganval(2)))^2 + (log(eiganval(2)/eiganval(1)))^2)^.5 ;
end