% Pressure variation in an ellipsoid in Linear Isotropic Material, with Pressure
% calcutions inside the ellipsoid, for different orientation, for different r.
%
function Pressasfuncofrdifforient(r)

%  Input parameters--------------------------------------------------------
%  Matrix flow field (Velocity field for simple shear flow)
   L     = [1 0 0;0 -1 0;0 0 0]; 


%  Power law coefficients for matrix and clast( since Linear case is
%  considered, both coefficients as 1)

   Nm=1;
   Nc=1;
   
%  Ratio of ellipsoid's viscosity to the matrix viscosity  ( r taken as an input)  
   
   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
   
%  generating 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Calculating Far field Stress Invariant for normalization :
   Sigma = 2* Multiply(Jd,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   Cm    = 2*Jd;                     % stiffness tensor for the matrix
   
% initial state of the ellipsoid
                                     % initial length of the ellipsoid axis
   a= [10;5;1];
   R=5;
                                     % initial angles of the ellipsoid axis
   
 epsilonII=.5;                      % initial strain rate invariant( This is useful when considering non-linear case)

  p1= zeros(1,180);
  p2=zeros(1,180);
  
  for i=1:180
      ang = [0; 0;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
      q = Q(ang);
    
      d= q*D*q';                                      % Transforming D and W to Ellipsoidal axis system
      w= q*W*q'; 
                     
   
      
    
  % [S,p,PI] = SnpIn(a,Jd,Js,Ja); this method uses quadratures. So, since
   % it is isotropic, so we use analytical results
   
   % Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
      [S_el,~] = SnP(a);

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
        E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
        r_new= E(1,4);
                                  % new r values being stored to see the variation of viscousity ratio
        dE= E(1:3,1:3);               % dE is strain rate tensor in ellipsoid
        u2= dE-d;
                      
        Pdev_in= contract1(p_in,u2);
        Sigma_ellipsoid = 2* Multiply(Jd,dE)* r_new;
        Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);
    
        p1(1,i) = Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
        p2(1,i) = Pdev_in/Sigma_Inva_ellipsoid ; 
    
  end
  
    
    
    ii=1:180;
    plot(ii, p2)
    



end



