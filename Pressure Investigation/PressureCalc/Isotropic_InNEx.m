%  function [e,s,w,p,e_Ex,s_Ex,w_Ex,p_Ex,pressure]=Isotropic_InNEx()
% Isotropic_InNEx.m
%
% Modeling the interior and exterior fields of a single shear zone, regarding
% a shear zone as an incompressible linear viscous ellipsoid.
%
% Scripts are based on formal solutions for interior and exterior fields of
% an incompressible linear viscous ellisoid.(Jiang, 2016 in press)
%
% (Assume both ellipsoid and matrix are isotropic.)
%--------------------------------------------------------------------------

% clear all variables, Comment Windows and figures
  clear;
  clc;
  clf;
  
% Input parameters:
% imposed far-field flow L
  L   = [1   0   0;...
         0   -1  0;...
         0   0   0]; 


% the viscosity ratio at matrix strain rate, eta_e/eta_matrix
  r     = 100/3;
% respect ratio of ellipsoid: a1 = 1, a2 = R*a1, a3 = 100R
%  R     = 20;
% semi-axes of the elliptic cylinder (a1 >> a2 > a3)
  a     = [5; 3; 1];  
% spherical angles in radian (Jiang,2007a);
  ang   = [0; 0; 0]; 
% imposed hydrostatic pressure P
%  P   = 1;

% 3D meshgrid in clast's coordinate
  xgv     = 0;           % grid vector: x'axis,a1
  ygv     = -5:0.1:5;           % grid vector: y'axis,a2
  zgv     = -1.5:0.1:1.5;    % grid vector: z'axis,a3
  [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
% Exterior points
  ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
  x_ex  = X(ind);
  y_ex  = Y(ind);
  z_ex  = Z(ind);
  ep    = cat(1,x_ex',y_ex',z_ex');
  
% total points number
  num   = numel(X);
% exterior points number
  [~,n] = size(ep);
% decompose the far-field flow L into a strain rate tensor D and a vorticity 
% tensor W, Eqn(3) in Jiang(2007a)     
  D     = 0.5 * (L + L');
  W     = 0.5 * (L - L');
% obtain the transformation matrix Q from three spherical angles, Eqns(8) 
% -(12) in Jiang(2007a)       
  q     = Q(ang);  
% generate 4th-order identity tensors   
  [Jd, Js, Ja, Jm] = FourIdentity();   
% stiffness tensor of the matrix, assuming it's isotropic here
  Cm    = 2*Jd;
% stiffness tensor of the ellipsoid, assuming it's isotropic here
  Ce    = 2*r*Jd;

  
% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,~,PI_el] = SnpIn(a,Jd,Js,Ja);

% Eshelby Tensors (S,PI,LAMDA) for Interior points
% LAMBDA
  LAMBDA   = zeros(3,3);
  for i=1:3
      LAMBDA(i,i) = -1/3* (S_el(i,i,1,1)+ S_el(i,i,2,2)+ S_el(i,i,3,3)); 
  end
% S  
  S       = S_el;
  for i=1:3
      for j=1:3
          S(i,i,j,j) = LAMBDA(i,i)+ S_el(i,i,j,j);
      end
  end
% PI  
  PI      = PI_el;
% Interior fields (e,w,p)

% describe D,W,C in the clast's coordinate system 
  D_bar  = q * D * q';
  W_bar  = q * W * q';
  Cm    = Transform(Cm, q); 
  Ce    = Transform(Ce, q);  
  S_bar  = Multiply(Cm, D_bar);
  SI     = Inva(S_bar);
  P      = SI;
% e: interior field strain-rate
  invS    = FourTensorInv(S);
  H       = Contract(Cm,(invS - Jd));   %Hill's constraint tensor
  t1      = FourTensorInv(H + Ce);
  A       = Contract(t1,(H + Cm));      %strain-rate partitioning tensor
  e       = Multiply(A, D_bar);
% s: interior field deviatoric-stress
  s       = Multiply(Ce, e);
% w: interior field vorticity
  de      = e - D_bar;                  
  t2      = Multiply(Contract(PI, invS), de);
  w       = t2 + W_bar;
% p: interior field pressure
  b1      = R_Multiply(LAMBDA, Cm);
  b2      = contract1(R_Multiply(b1, invS), de);
  p       = b2 + P;

% Exterior Fields  
%--------------------------------------------------------------------------
% G
  G           = Ex_Gtensor(a,ep);
% LMABDA_Ex 
  LAMBDA_Ex   = zeros(3,3,n);
  for i=1:3
      for j=1:3
          LAMBDA_Ex(i,j,:) = -1/3.*(G(i,j,1,1,:)+ G(i,j,2,2,:)+ G(i,j,3,3,:));
          LAMBDA_Ex(j,i,:) = LAMBDA_Ex(i,j,:);
      end
  end
  for i=1:3
      t                = 1/3.*squeeze(G(1,1,i,i,:)+ G(2,2,i,i,:)+ G(3,3,i,i,:));
      LAMBDA_Ex(i,i,:) = squeeze(LAMBDA_Ex(i,i,:))+t;
  end
% S_Ex & PI_Ex
  S_Ex   = zeros(3,3,3,3,n);
  PI_Ex  = zeros(3,3,3,3,n);
  delt   = eye(3);
  for i=1:3
      for j=i:3
          for k=1:3
              for l=k:3
                  %S_Ex
                  S_Ex(i,j,k,l,:)  = squeeze(G(i,j,k,l,:))+ delt(k,l).*...
                                     squeeze(LAMBDA_Ex(i,j,:));
                  S_Ex(j,i,k,l,:)  = S_Ex(i,j,k,l,:);
                  S_Ex(j,i,l,k,:)  = S_Ex(i,j,k,l,:);
                  S_Ex(i,j,l,k,:)  = S_Ex(i,j,k,l,:);
                  %PI_Ex
                  PI_Ex(i,j,k,l,:) = 1/2.*(delt(j,k).*LAMBDA_Ex(i,l,:) +...
                                     delt(j,l).*LAMBDA_Ex(i,k,:) - delt(i,k)...
                                     .*LAMBDA_Ex(j,l,:) - delt(i,l).*...
                                     LAMBDA_Ex(j,k,:));
                  PI_Ex(i,j,l,k,:) = PI_Ex(i,j,k,l,:);
                  PI_Ex(j,i,k,l,:) = -PI_Ex(i,j,k,l,:);
                  PI_Ex(j,i,l,k,:) = -PI_Ex(i,j,k,l,:);
              end
          end
      end
  end
  S_Ex(1,1,1,1,:) = -(S_Ex(1,1,2,2,:)+S_Ex(1,1,3,3,:));
  S_Ex(2,2,2,2,:) = -(S_Ex(2,2,1,1,:)+S_Ex(2,2,3,3,:));
  S_Ex(3,3,3,3,:) = -(S_Ex(3,3,1,1,:)+S_Ex(3,3,2,2,:));
  
% Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
  e_Ex    = zeros(3,3,n);
  s_Ex    = zeros(3,3,n);
  w_Ex    = zeros(3,3,n);
  p_Ex    = zeros(1,n);
 
 for r=1:n
% e_Ex: exterior field strain-rate in clast's coordinate
     v1          = Contract(S_Ex(:,:,:,:,r), invS);
     v2          = Multiply(v1, de);
     e_Ex(:,:,r) = v2 + D_bar;
% s_Ex: exterior field stress in clast's coordinate
     s_Ex(:,:,r) = Multiply(Cm, e_Ex(:,:,r));
    
% w_Ex: exterior field vorticity in clast's coordinate
     v3          = Contract(PI_Ex(:,:,:,:,r), invS);
     v4          = Multiply(v3, de);
     w_Ex(:,:,r) = v4 + W_bar;
% p_Ex: exterior field pressure
     u1          = R_Multiply(LAMBDA_Ex(:,:,r), Cm);
     u2          = R_Multiply(u1, invS);
     u3          = contract1(u2, de);
     p_Ex(r)     = u3 + P;
    
 end
     p1          = p / SI;
     p2          = p_Ex ./ SI;
     
 % deviatoric stress invariant plot
  
   
   pressure        = zeros(num,1);
   pressure(~ind)  = p1;
   pressure(ind)   = p2;
   pressure        = reshape(pressure,size(squeeze(Z)));
  
   contourf(squeeze(Y),squeeze(Z),pressure)
   axis equal
        

%end