function plotTvsGP(gp)

k=pwd;
K= extractBefore(k,"Pressure Investigation");
addpath(genpath(K)); 
a= [5;3;2];   % inclusion shape
EP=[0;0;4];   % an external point 


% for Anisotropic case
[Jd, ~, ~, ~, ~] = Jnb(); 
%C= 2*Jd;  % Isotropic medium

C    = zeros(3,3,3,3); % stiffness tensor
m= 10;    % m, is the strength of planar anisotropy, m= 1 gives the isotropic case

Nn=1;

% defining stiffness tensor for material with planar anisotropy.
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    C(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    C(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end
 
%--Evaluating Gaussian grid nodes and weights----------------------
[Alpp1,Bett1,www1] = GaussGGQ(gp);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(gp);
%-------------------------------------------------------------------

% Using Non-vectorized code

tstart1= tic;
Lambda1 = LambdaExtAniso3DPress_nv(a,EP,C,theta1,phi1,Wout1,psi1,Win1);
tend1= toc(tstart1);

% Using Vectorized code

tstart2= tic;
Lambda2 = LambdaExtAniso3DPress_vec(a,EP,C,theta1,phi1,Wout1,psi1,Win1);
tend2= toc(tstart2);

% Note the elapsed time for
%  Non-vectorized code : 8.62 seconds
%  Vectorized code:  2.36 seconds
scatter(gp,tend1)
hold on
scatter(gp,tend2,'*')
hold on
end