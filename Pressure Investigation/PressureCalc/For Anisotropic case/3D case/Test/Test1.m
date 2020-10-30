

% External Eshelby Tensor Test
% test for exterior eshelby tensor for stress.
% This code compares the code for anisotropic medium with the isotropic
% quasi-analytical solutions.

k=pwd;
K= extractBefore(k,"Pressure Investigation");
addpath(genpath(K)); 
a= [3;2;1];   % inclusion shape
EP=[0;3;0];   % external coordinate in the inclusion coordinate system

lambda = solve_eq(a, EP);

%
% for Isotropic case
 tic             
        G1point = Ex_Gtensor(a,EP);
        
        % LMABDA_Ex 
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

% for Anisotropic case
[Jd, Js, Ja, Jm, b] = Jnb(); 
C= 2*Jd;

[Alpp1,Bett1,www1] = GaussGGQ(60);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(25);

tic
T = T_ExtAniso3D(a,EP,C,theta1,phi1,Wout1,psi1,Win1);
S_Ex_Aniso= Contract(Contract(Jd,T),C);
toc

%}