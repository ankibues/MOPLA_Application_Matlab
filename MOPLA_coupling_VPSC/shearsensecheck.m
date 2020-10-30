% Shear sense check in the inclusion, due to flow field partitioning !
%  Input parameters--------------------------------------------------------
%  Matrix flow field (Velocity field for simple shear flow)
   L     = [0 2 0;0 0 0;0 0 0]; 

%  Power law coefficients for matrix and clast( since Linear case is
%  considered, both coefficients as 1)
   Nm=3;
   Nc=3;
%  Ratio of ellipsoid's viscosity to the matrix viscosity   
   r= .01;
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
   Wvec= [2*W(3,2); 2*W(1,3);2*W(2,1)];  % vorticity vector
   zaxis= [0;0;1];
   dotprod1= dot(Wvec,zaxis);
%  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 
   
% Calculating Far field Stress Invariant for normalization :
   Sigma = 2* Multiply(Jd,D);   % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   strainInva = inva(D);
   Cm    = 2*Jd;                % stiffness tensor for the matrix
   
% initial state of the ellipsoid
  a= [5;3;1];                 % initial length of the ellipsoid axis    
  dotprodcheck= zeros(180,1);
for i=1:180
    ang = [0;0;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
    q = Q(ang);
    rot_z= q*zaxis;
    d= q*D*q';                                      % Transforming D and W to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);                    
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
   % [S,p,PI] = SnpIn(a,Jd,Js,Ja); this method uses quadratures. So, since
   % it is isotropic, so we use analytical results
   
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
% PI  
  PI      = PI_el;
    invS= Inverse(S,b);
    RR= R_Multiply(p,Inverse((S-Jd),b));   
    h= Contract(Cm_a/Nm, invS);                        
    p_in = R_Multiply(p,h);
    dE       = Multiply(fdE(Nm,r,S,b,Jd),d);
 
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    LL = dE + we;
    wEE= q'*we*q;
    wvec= [2*wEE(3,2); 2*wEE(1,3);2*wEE(2,1)]; % vorticity vector
    dotprod= dot(wvec,zaxis);
    dotprodcheck(i,1)= dotprod; 
  
end
  
    

