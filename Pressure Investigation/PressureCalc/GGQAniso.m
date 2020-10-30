function [T, lambda] = GGQAniso(sub,x,ang,Nn,m)
% x is 2*1 matrix with semi axises of ellipse
a=x(1);
b=x(2);
[psi,w]=lgwt(100,0,pi);

phi=ang;

z= [cos(psi'); sin(psi')];
del= [1,0;0,1];
epsilon = (m*(sin(2*(phi + psi'))).^2 + (cos(2*(phi + psi'))).^2).^(-1);
epsilon_rep= repmat(epsilon, 2,1)/2;
C = epsilon_rep.*[(m+1)*cos(psi')-(m-1)*cos(4*phi + 3*psi'); (m+1)*sin(psi')+(m-1)*sin(4*phi + 3*psi')];

for o=1:numel(epsilon)
    for i=1:2
        for j= 1:2
            Abar(i,j,o) = (m*epsilon(o)/Nn)*(del(i,j) -z(i,o)*z(j,o));
        end
    end
end

% Green Interaction Tensor for  inside the ellipse
p= 1./((a*cos(psi')).^2 + (b*sin(psi')).^2);
for o=1:numel(epsilon)
    f(o)= Abar(sub(1),sub(3),o)*z(sub(2),o)*z(sub(4),o)*p(o);
end
T= (a*b/pi) * sum (f.*w');

% Green tensor for pressure calculation

g= C(sub(1),:).*z(sub(2),:).*p;
lambda= (-a*b/pi)* sum (g.*w');
end


