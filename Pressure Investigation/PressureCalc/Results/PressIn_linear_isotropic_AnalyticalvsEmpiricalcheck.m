% Pressure variation in an ellipsoid in Linear Isotropic Material, with Pressure
% calcutions inside the ellipsoid, for different instantatneous angles of orientation.
%

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 2 0;0 0 0;0 0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 1000;
%  number of computation steps between output sets 
   mm=10; 
   
%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1;
   r= 1000;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 
   
% for stress invariant for normalization :
   Sigma = 2* Multiply(Jd,D);                 % far field stress value( sigma= 2*Jd:E)i.e Viscosity of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   strainInva = inva(D);
   Cm    = 2*Jd;
% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                      % initial length of the ellipsoid axis
   a= [1000;5;1];
   R=5;
                                            % initial angles of the ellipsoid axis
   
   epsilonII=0.5;                                                      % initial strain rate invariant
% calculations to follow:
 Pdev_steps_alpha1= zeros(180,1);
 Pdev_steps_alpha2= zeros(steps,1);
 Pdev_steps_analytical= zeros(180,1);
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  strainininvacheck= zeros(1,steps);
  sigmaininvacheck= zeros(1,steps);
  r_new_check= zeros(steps,1);
  SI_check = zeros(steps,1);
for i=1:180
    
    ang = [0; 0;(pi/180)*i];
    q = Q(ang);
    
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    SI_check(i,1)= SI;
    % w is the vorticity of material, in ellipsoid's frame of reference 
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
    h= Contract(2*Jd, invS);                        
    p_in = R_Multiply(p,h);
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    r_new= E(1,4);
    r_new_check(i,1)= r_new;                         % new r values being stored to see the variation of viscousity ratio
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    sigmainside= 2*dE*r_new;
    strainininva = inva(dE);
    sigmaininva= inva(sigmainside);
    
    
    Pdev_in= contract1(p_in,u2);
    Sigma_ellipsoid = 2* Multiply(Jd,dE)* r_new;
    Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);
    sigmaininvacheck(1,i) = sigmaininva;   
    strainininvacheck(1,i)= strainininva;
    
    Pdev_analytical = 2*d(2,2)*(1-r_new)*((R^2)-1)/((R^2) +2*r_new*R+1);
    
    Pdev_steps_analytical(i,1)= Pdev_analytical/Sigma_Inva;
    Pdev_steps_alpha1(i,1)= Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
    Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipsoid ; 
    
    
    
    
end 

x=1:180;
plot(x,Pdev_steps_alpha1,'bo')
xlabel('Angle in degress');
ylabel('Pin/SI');
hold on
plot(x,Pdev_steps_analytical,'r')
    

