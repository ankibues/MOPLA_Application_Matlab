function [T, lambda] = GGQAnisoExt(sub,aa,x,phi,Nn,m)
% aa is 2*1 matrix with semi axises of ellipse
% x is 2*1 matrix for x, and y coordinates 
a=aa(1);
b=aa(2);
[psi,w]=lgwt(50,0,pi);


z= [cos(psi'); sin(psi')];

del= [1,0;0,1];
epsilon = (m*(sin(2*(phi + psi'))).^2 + cos(2*(phi + psi')).^2).^(-1);
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
k= zeros(1,numel(epsilon));

for o=1:numel(epsilon)
    k(o)= real(1- (i*(dot(x,z(:,o)))/((a*z(1,o))^2 + (b*z(2,o))^2 - (dot(x,z(:,o)))^2  )^0.5));
end

B= p.* k;

for o=1:numel(epsilon)
    f(o)= Abar(sub(1),sub(3),o)*z(sub(2),o)*z(sub(4),o)*B(o);
end
T= (a*b/pi) * sum (f.*w');

% Green tensor for pressure calculation

g= C(sub(1),:).*z(sub(2),:).*B;
lambda= (-a*b/pi)* sum (g.*w');
end


