function Z= DT_vec(a,x,theta,phi,C,psi,ww1)
% a is 3*1, x is 3*1, theta is 1*N, phi is 1*N, C is 3*3*3*3
% Output Z will be 3*3*3*3*N
% 
[~,N]=size(theta);
Z=zeros(3,3,3,3,N);
xi= [ cos(theta).*sin(phi); sin(theta).*sin(phi);cos(phi)]; %3*N
xp= [a(1).*cos(theta).*sin(phi); a(2).*sin(theta).*sin(phi); a(3).*cos(phi)]; %3*N
A=xi./a; % 3*N
r= x-xp; % 3*N
g = Gm_vec(r,C,psi,ww1); %here, Gm_vec should give 3*3*3*N
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
               Z(i,j,k,l,:)=squeeze(g(i,k,l,:))'.*squeeze(A(j,:));
            end
        end
    end
end
Z= Z.*a(1).*a(2).*a(3);
end