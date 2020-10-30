
k=pwd;
K= extractBefore(k,"Pressure Investigation");
addpath(genpath(K)); 
a= [3;2;1];   % inclusion shape
EP=[0;0;4];   % an external point 
lambda = solve_eq(a, EP);

% Using quasi-analytical solutions for Isotropic case
% as a validation for the anisotropic code---------------------------------
% tic             
        G1point = Ex_Gtensor(a,EP);
        
        % LAMBDA_Ex 
        LAMBDA_Ex   = zeros(3,3);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j) = -1/3.*(G1point(ii,j,1,1)+ G1point(ii,j,2,2)+ G1point(ii,j,3,3));
                LAMBDA_Ex(j,ii) = LAMBDA_Ex(ii,j);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G1point(1,1,iii,iii)+ G1point(2,2,iii,iii)+ G1point(3,3,iii,iii));
            LAMBDA_Ex(iii,iii) = squeeze(LAMBDA_Ex(iii,iii))+t;
        end
        

%toc
%--------------------------------------------------------------------------
% for Anisotropic case
[Jd, Js, Ja, Jm, b] = Jnb(); 
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
[Alpp1,Bett1,www1] = GaussGGQ(10);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(10);
%-------------------------------------------------------------------

% Using Non-vectorized code
tic
tstart1= tic;
Lambda1 = LambdaExtAniso3DPress_nv(a,EP,C,theta1,phi1,Wout1,psi1,Win1);
tend1= toc(tstart1);
toc
% Using Vectorized code
tic
tstart2= tic;
Lambda2 = LambdaExtAniso3DPress_vec(a,EP,C,theta1,phi1,Wout1,psi1,Win1);
tend2= toc(tstart2);
toc
% Note the elapsed time for
%  Non-vectorized code : 8.62 seconds
%  Vectorized code:  2.36 seconds