function Plot_ep_X(a)

%Plotting pressure for exterior points along X-axis 

% shape and orientation of inclusion
%a=[10;2;1]; % a1, a2 and a3 are the semi-axes with a1>=a2>=a3     
ang = [0; 0; 0];  
% these spherical angles correspond to the case when a1 lie along Z-axis, a2 along X-axis and a3 along Y-axis 
q = Q(ang);

%exterior points
x=0:.5:8;
%y=0:.5:10;
Xpts=[x;zeros(1,length(x));zeros(1,length(x))];    % points along X-axis
%Ypts=[zeros(1,length(x));y;zeros(1,length(x))];    % points along Y-axis

% note, these points are defined in coordinate system of the flow, and
% should be converted to inclusion's coordinate system for calculation of
% exterior fields.
Xpts_inc= zeros(size(Xpts));
%Ypts_inc= zeros(size(Ypts));
for i=1:length(Xpts)
    Xpts_inc(:,i)= q*Xpts(:,i);
    %Ypts_inc(:,i)= q*Ypts(:,i);
end


 ind   = (Xpts_inc(1,:)./a(1)).^2 + (Xpts_inc(2,:)./a(2)).^2 + (Xpts_inc(3,:)./a(3)).^2 > 1;
 EP=Xpts_inc(:,ind);
 
 N=length(EP);
 
 
L     = [1 0 0;0 -1 0;0 0 0];  % Macroscale flow field given by velocity gradient tensor
D = 0.5 * (L + L');
d= q*D*q';        
[Jd, Js, Ja, Jm, b] = Jnb(); % 4th order identitiy tensors
Cm= 2*Jd;   % stiffness matrix for the isotropic medium
 Nm=1;      % power law exponent
 r=10;      % viscosity ratio of the inclusion's viscosity to the embedding medium.
 
 
 %-------------------------------------------------------------------------------
 %--------------Quasi-analytical solution for isotropic medium------------
 %-------------------------------------------------------------------------------
 
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
  dE       = Multiply(fdE(Nm,r,S,b,Jd),d);  % dE is strain rate  in the inclusion              
                                                                     
u2= dE-d;
 %----------------------------------------------------------------------
 % some exterior point[ Important note: these points are defined in
 % inclusion's coordinate system ]
 P_anal= zeros(1,N);
 P_num= zeros(1,N);
 parfor kkk=1:N
     ep= EP(:,kkk);
 
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
            p_Ex    = uu3 ; %pressure at the exterior point(ep) using isotropic formulation
            P_anal(1,kkk)= p_Ex;
 %-------------------------------------------------------------------------------
 

 
        
 %--------------Calculation using code for the general Anisotropic medium.............
        
        
 %--------------------------------------------------------------------------
% Gaussian point calculation

           [Alpp1,Bett1,www1] = GaussGGQ(300);
           theta1 = reshape(Alpp1,1,[]);
           phi1   = reshape(Bett1,1,[]);
           Wout1 = reshape(www1,1,[]);
           [psi1,Win1] = Gauss1(25);

           [Alpp2,Bett2,www2] = GaussGGQ(150);
           theta2 = reshape(Alpp2,1,[]);
           phi2   = reshape(Bett2,1,[]);
           Wout2 = reshape(www2,1,[]);
           [psi2,Win2] = Gauss1(25);

           [Alpp3,Bett3,www3] = GaussGGQ(100);
           theta3 = reshape(Alpp3,1,[]);
           phi3   = reshape(Bett3,1,[]);
           Wout3 = reshape(www3,1,[]);
           [psi3,Win3] = Gauss1(25);

           [Alpp4,Bett4,www4] = GaussGGQ(100);
           theta4 = reshape(Alpp4,1,[]);
           phi4   = reshape(Bett4,1,[]);
           Wout4 = reshape(www4,1,[]);
           [psi4,Win4] = Gauss1(25);
  
           % Exterior fields (p_Ex)
  
            LAMBDA_Ex2= zeros(3,3,1);
  
 
     
            tic
        % --------------------------------------------------------------
             lambda = solve_eq(a, ep);
             if lambda<= .64
                 LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_vec(a,ep,Cm,theta1,phi1,Wout1,psi1,Win1);
             elseif lambda>.64 && lambda<= 1.36
                 LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_vec(a,ep,Cm,theta2,phi2,Wout2,psi2,Win2);
             elseif lambda> 1.36 && lambda<= 2.16
                 LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_vec(a,ep,Cm,theta3,phi3,Wout3,psi3,Win3);
             else 
                 LAMBDA_Ex2(:,:,1) = LambdaExtAniso3DPress_vec(a,ep,Cm,theta4,phi4,Wout4,psi4,Win4);
             end
            % This approach allows us to choose number of Gaussian points 
            % and grid nodes based on exterior point's position from the 
            % inclusion-medium boundary-------------------------------------
         
        
            % p_Ex: exterior field pressure
              uu1          = R_Multiply(LAMBDA_Ex2(:,:,1),Cm);
              uu2          = R_Multiply(uu1, invS);
              uu3          = contract1(uu2, u2);
              p_Ex2     = uu3 ; %pressure at an exterior point using Anisotropic formulation
              P_num(1,kkk)= p_Ex2;
              
             toc
             
 end
 XX= EP(2,:); 
 
 
 scatter(XX, P_anal,'*')
 hold on
 scatter(XX, P_num,'s')
 
 
end