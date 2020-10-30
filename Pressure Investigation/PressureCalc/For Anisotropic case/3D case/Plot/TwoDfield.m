function KK= TwoDfield(pressure)
[n1,n2]= size(pressure);
% Code for converting 1/4th part of the pressure field to whole 
P1= flip(pressure,1);
K= cat(1,P1,pressure);
K(n1, :,:)=[];
P2= flip(K,2);
KK= cat(2,P2,K);
KK(:,n2,:)=[];

end