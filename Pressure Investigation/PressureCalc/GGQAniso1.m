function [T, lambda] = GGQAniso1(sub,x,ang,Nn,m)
% x is 2*1 matrix with semi axises of ellipse
% Using matlab integral function for calculation to check the results


a=x(1);
b=x(2);
phi=ang;
% Green Interaction Tensor for  inside the ellipse

f=@(psi)Abar(m,Nn,phi,psi,sub(1),sub(3)).*z(psi,sub(2)).*z(psi,sub(4)).*B(psi,a,b);
T= (a*b/pi)*integral(f,0,pi);

% Green tensor for pressure calculation

g= @(psi) C(m,phi,psi,sub(1)).*z(psi,sub(2)).*B(psi,a,b) ;
lambda = (-a*b/pi)* integral(g,0,pi);
end

function Z= z(psi,a)
zz= [cos(psi); sin(psi)];
Z= zz(a);
end
function eps= epsilon(m,phi,psi)
eps = (m*(sin(2*(phi + psi))).^2 + cos(2*(phi + psi)).^2).^(-1);
end
function P= B(psi,a,b)
P = 1./((a*cos(psi)).^2 + (b*sin(psi)).^2);
end
function A= Abar(m,Nn,phi,psi,i,j)
del= [1,0;0,1];

A = (m*epsilon(m,phi,psi)/Nn)*(del(i,j)- z(psi,i)*z(psi,j));

end
function c= C(m,phi,psi,i)
cc = [(m+1)*cos(psi)-(m-1)*cos(4*phi + 3*psi); (m+1)*sin(psi)+(m-1)*sin(4*phi + 3*psi)];
c=  (epsilon(m,phi,psi)/2).* cc(i);
end



