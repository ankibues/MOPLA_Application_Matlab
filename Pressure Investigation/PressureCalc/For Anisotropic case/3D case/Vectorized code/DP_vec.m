function Z= DP_vec(a,x,theta,phi,C,psi,ww1)
% a is 3*1 matrix ,3 semi-axes of an inclusion;
% x is 3*1, 3 coordinates of the external point 
% theta  is 1*N1^2 matrix containing Gaussian grid nodes corresponding to theta, 
% phi is 1*N1^2 matrix containing Gaussian grid nodes corresponding to phi
% Here, N1 is the number of nodes and weights; The total grid nodes becomes N1^2.
% C is the stiffness tensor of the embedding medium
%  psi is N2*1 matrix containing Gaussian grid node corresponding to psi
%  ww1 is N2*1 matrix containing Gaussian weights corresponding to psi
%  Here, N2 is the number of nodes and weights;

% Output Z will be 3*3*(N1^2)

[~,N]=size(theta);
Z=zeros(3,3,N);
xi= [ cos(theta).*sin(phi); sin(theta).*sin(phi);cos(phi)]; %3*N
xp= [a(1).*cos(theta).*sin(phi); a(2).*sin(theta).*sin(phi); a(3).*cos(phi)]; %3*N
A=xi./a; % 3*N     % this is corresponding to Ej/aj in the integrand
r= x-xp; % 3*N
h = Hm_vec(r,C,psi,ww1);
for i=1:3
    for j=1:3
        Z(i,j,:)=h(i,:).*A(j,:);
    end
end
Z= Z.*a(1).*a(2).*a(3);
end