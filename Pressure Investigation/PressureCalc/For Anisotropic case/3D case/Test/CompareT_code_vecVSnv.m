
k=pwd;
K= extractBefore(k,"Pressure Investigation");
addpath(genpath(K)); 
a= [3;2;1];   % inclusion shape
EP=[0;5;0];   % an external point 
lambda = solve_eq(a, EP);

% Using quasi-analytical solutions for Isotropic case
% as a validation for the anisotropic code---------------------------------
 tic             
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
        
        % S_Ex & PI_Ex
        S_Ex   = zeros(3,3,3,3);
        PI_Ex  = zeros(3,3,3,3);
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex(ii,j,k,l)  = squeeze(G1point(ii,j,k,l))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex(ii,j));
                        S_Ex(j,ii,k,l)  = S_Ex(ii,j,k,l);
                        S_Ex(j,ii,l,k)  = S_Ex(ii,j,k,l);
                        S_Ex(ii,j,l,k)  = S_Ex(ii,j,k,l);
                        %PI_Ex
                        PI_Ex(ii,j,k,l) = 1/2.*(delt(j,k).*LAMBDA_Ex(ii,l) +...
                                            delt(j,l).*LAMBDA_Ex(ii,k) - delt(ii,k)...
                                            .*LAMBDA_Ex(j,l) - delt(ii,l).*...
                                            LAMBDA_Ex(j,k));
                        PI_Ex(ii,j,l,k) = PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,k,l) = -PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,l,k) = -PI_Ex(ii,j,k,l);
                    end
                end
            end
        end
        S_Ex(1,1,1,1) = -(S_Ex(1,1,2,2)+S_Ex(1,1,3,3));
        S_Ex(2,2,2,2) = -(S_Ex(2,2,1,1)+S_Ex(2,2,3,3));
        S_Ex(3,3,3,3) = -(S_Ex(3,3,1,1)+S_Ex(3,3,2,2));

toc
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
[Alpp1,Bett1,www1] = GaussGGQ(20);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(20);
%-------------------------------------------------------------------

% Using Non-vectorized code
tic
T1 = T_ExtAniso3D(a,EP,C,theta1,phi1,Wout1,psi1,Win1);   % Green's tensor
S_Ex_Aniso1= Contract(Contract(Jd,T1),C);                % Exterior Eshelby tensor
toc

% Using Vectorized code
tic
T2 = T_ExtAniso3D_vec(a,EP,C,theta1,phi1,Wout1,psi1,Win1);  % Green's tensor
S_Ex_Aniso2= Contract(Contract(Jd,T2),C);                   % Exterior Eshelby tensor
toc

% Note the elapsed time for
%  Non-vectorized code :176.04 seconds
%  Vectorized code: 24.62 seconds
