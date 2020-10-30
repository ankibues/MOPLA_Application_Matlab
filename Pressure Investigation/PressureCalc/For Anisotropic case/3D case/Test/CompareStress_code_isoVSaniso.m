
% Test for exterior eshelby solution for deviatoric stress  
k=pwd;
K= extractBefore(k,"Pressure Investigation");
addpath(genpath(K)); 


%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 0;0 -1 0;0 0 0]; 
   
   
%  Power law coefficients for matrix and clast
   Nm = 1;
     % since linear isotropic case is being considered
   r= 10;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
   C= 2*Jd ; % Viscous stiffness tensor of the matrix
 

% geometry of the inclusion   
   a= [3;1.5;1];  % length of the ellipsoid axis(inclusion shape)
   ang = [0; 0;(pi/180)*0];      % orientation of the inclusion
   q = Q(ang);
   
    d= q*D*q';              % Transforming D and W with respect to Ellipsoidal axis system
   
   
% Solutions using quasi-analytical solutions after Jiang (2016)
tic
% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,PI_el] = SnP(a);

% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green tensor for pressure)
  p   = zeros(3,3);
  for j=1:3
      p(j,j) = -1/3* (S_el(j,j,1,1)+ S_el(j,j,2,2)+ S_el(j,j,3,3)); 
  end
% S  
  S       = S_el;
  for k=1:3
      for l=1:3
          S(k,k,l,l) = p(k,k)+ S_el(k,k,l,l);
      end
  end
  
invS= Inverse(S,b);
h= Contract(2*Jd, invS);                        
dE       = Multiply(fdE(Nm,r,S,b,Jd),d);    %  strain rate tensor in inclusion                  
u2= dE-d; 

%----------------------------------------------------------------------
 % some exterior point
 
 ep= [0;0;3];
 
    % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,ep);
        % LMABDA_Ex 
        LAMBDA_Ex   = zeros(3,3);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j) = -1/3.*(G(ii,j,1,1)+ G(ii,j,2,2)+ G(ii,j,3,3));
                LAMBDA_Ex(j,ii) = LAMBDA_Ex(ii,j);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G(1,1,iii,iii)+ G(2,2,iii,iii)+ G(3,3,iii,iii));
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
                        S_Ex(ii,j,k,l)  = squeeze(G(ii,j,k,l))+ delt(k,l).*...
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
  
       

      
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex, invS);
            v2          = Multiply(v1, u2);
            e_Ex = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            sigma_out1 = Multiply(C, e_Ex);     % deviatoric stress at an exterior point around the inclusion using isotropic solutions
            
        
  toc  
    
     %--------------------USING ANISOTROPIC SOLUTIONS-----------------------
   % Gaussian points   
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
  
      
  
   [Jd, Js, Ja, Jm, b] = Jnb(); 



 
  
% stiffness tensor of the inclusion, 
   Ce    = 2*r*Jd; 
   Sigma = Multiply(C,D); % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);


    Cm_a    = Transform(C, q);
    Carray  = C2OneDarray(Cm_a);
    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
    S = Contract(Jd, Contract(T,Cm_a));
    
    lambda = SnpInAniso3DPress(a,Cm_a);
    invS = Inverse((S),b);
    h= Contract(Cm_a, invS);                        
    p_in = R_Multiply(lambda,h);
    H = Contract(Cm_a,(invS-Js));           % Eq.12a Jiang 2014
   
    A = Contract(FourTensorInv(H + Ce),(H + Cm_a));                     
    e= Multiply(A,d);    % e is strain rate in ellipsoid
    u1= e-d;                  
    sigmainside= Multiply(Ce,e);
    sigmaininva= inva(sigmainside);
       
    Pdev_in= contract1(p_in,u1);
    
        
 lambda = solve_eq(a, ep);  % function for calculating lamda value for each point
%
        % Exterior Fields  
 %--------------------------------------------------------------------------
% Gaussian point calculation

[Alpp1,Bett1,www1] = GaussGGQ(25);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(25);

% Exterior fields 
     
       
  tic
    
        % Lambda_External 
        
         T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta1,phi1,Wout1,psi1,Win1);
        
        % S_Ex: exterior Eshelby tensor     strain rate e_Ex
            S_Ex2 = Contract(Contract(Jd,T_Ex),Cm_a);
            hh=  Contract(S_Ex2, invS);
            e_Ex2= Multiply(hh,u1)+ d;
           
            sigma_out2= Multiply(Cm_a,e_Ex2); % deviatoric stress at an exterior point around the inclusion using anisotropic solutions
  toc
       
       
        
disp('Stress using isotropic formulation =');
disp(sigma_out1);

disp('Stress using anisotropic formulation =');
disp(sigma_out2);






