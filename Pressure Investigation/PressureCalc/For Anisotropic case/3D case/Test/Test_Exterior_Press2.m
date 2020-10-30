% Test of Exterior Pressure calculation Analytical vs Numerical 3D
% Anisotropic case !
 L     = [1 0 0;0 -1 0;0 0 0]; 
  D = 0.5 * (L + L');
   ang = [0; 0; 0];
    q = Q(ang);
    d= q*D*q'; 
a=[3;1.5;1];
 [Jd, Js, Ja, Jm, b] = Jnb(); 
 Cm= 2*Jd;
 Nm=1;
 r=.1;
 %{
 Cm1    = zeros(3,3,3,3);
m= 1;
Nn=1;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end
 %}
%-------------------------------------------------------------------------
 % Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,PI_el] = SnP(a);

% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green
% tensor for pressure)
% p
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
    p_in = R_Multiply(p,h);
    dE       = Multiply(fdE(Nm,r,S,b,Jd),d);                  
                                                                     % dE is strain rate tensor in ellipsoid
    
    u2= dE-d;
 %----------------------------------------------------------------------
 % some exterior points
 
 ep= [0;0;2];
 % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,ep);
        % LMABDA_Ex 
        LAMBDA_Ex   = zeros(3,3,1);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j,:) = -1/3.*(G(ii,j,1,1,:)+ G(ii,j,2,2,:)+ G(ii,j,3,3,:));
                LAMBDA_Ex(j,ii,:) = LAMBDA_Ex(ii,j,:);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G(1,1,iii,iii,:)+ G(2,2,iii,iii,:)+ G(3,3,iii,iii,:));
            LAMBDA_Ex(iii,iii,:) = squeeze(LAMBDA_Ex(iii,iii,:))+t;
        end
        
 
       
       
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex(:,:,1), Cm);
            uu2          = R_Multiply(uu1, invS);
            uu3          = contract1(uu2, u2);
            p_Ex    = uu3 ; 
    
     
        
        %--------------FOr anisotropic formulation.............
        
        
        %--------------------------------------------------------------------------
% Gaussian point calculation

[Alpp1,Bett1,www1] = GaussGGQ(120);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(75);

[Alpp2,Bett2,www2] = GaussGGQ(75);
theta2 = reshape(Alpp2,1,[]);
phi2   = reshape(Bett2,1,[]);
Wout2 = reshape(www2,1,[]);
[psi2,Win2] = Gauss1(15);

[Alpp3,Bett3,www3] = GaussGGQ(30);
theta3 = reshape(Alpp3,1,[]);
phi3   = reshape(Bett3,1,[]);
Wout3 = reshape(www3,1,[]);
[psi3,Win3] = Gauss1(30);

[Alpp4,Bett4,www4] = GaussGGQ(30);
theta4 = reshape(Alpp4,1,[]);
phi4   = reshape(Bett4,1,[]);
Wout4 = reshape(www4,1,[]);
[psi4,Win4] = Gauss1(30);


%-----------------------------------------------------------------



  
        % Exterior fields (p_Ex)
     lambda = solve_eq(a, ep);
        LAMBDA_Ex2= zeros(3,3,4);
        %p_Ex2    = zeros(1,4);
 
     
        tic
        % Lambda_External 
         if lambda<= .64
             LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_nv(a,ep(:,1),Cm,theta1,phi1,Wout1,psi1,Win1);
         elseif lambda>.64 && lambda<= 1.36
             LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_nv(a,ep(:,1),Cm,theta2,phi2,Wout2,psi2,Win2);
         elseif lambda> 1.36 && lambda<= 2.16
             LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_nv(a,ep(:,1),Cm,theta3,phi3,Wout3,psi3,Win3);
         else 
             LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_nv(a,ep(:,1),Cm,theta4,phi4,Wout4,psi4,Win4);
         end
         
        
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex2(:,:,1),Cm);
            uu2          = R_Multiply(uu1, invS);
            uu3          = contract1(uu2, u2);
            p_Ex2     = uu3 ; 
    toc
       
        
        