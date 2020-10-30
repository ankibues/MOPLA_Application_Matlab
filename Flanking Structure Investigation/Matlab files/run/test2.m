%function Single_marker_FS(Wk,r,a,pp)% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure

k=pwd;
 K= extractBefore(k,"Flanking Structure Investigation");
 addpath(genpath(K));
%parpool(32);
tic

 Wk= 1;
  gamma=1;
  epsilon= .5*(((1/(Wk^2))-1)^(.5));
 if Wk==0
     gamma=0;
     epsilon=1;
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 else
     L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
 end
 
%  Input parameters--------------------------------------------------------
%  Matrix flow field
  % L     = [0 1 0;0 0 0;0 0 0]; 
D = 0.5 * (L + L');
% normalizing with respect to strain rate
L= L/norm(D);


%  Power law coefficients for matrix and inclusion
   Nm=1;
   Nc=1; 
   r= .1;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Stiffness tensor of the matrix   
   Cm    = 2*Jd;
   
% initial state of the ellipsoid                                    
   a= [1;.5;.1];
   A=a;% initial length of the ellipsoid axis
   pp= 0;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
 
% Here, I define the markers in the external coordinate system

%
        
        
        
    

   
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    
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
    h= Contract(2*Jd/Nm, invS);           % Here is the Cinverse, divided by Nm             
    p_in = R_Multiply(p,h);
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    r_new= E(1,4);
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    ll= dE + we;                                    % Partitioned velocity gradient tensor

    X= 0;
    Y= .6;
    Z= .1;
    V1 = Velocity_field_calc5(Y,Z,a,invS,Jd,d,u2,w,ll);
    V2 = Velocity_field_calc_general3D(X,Y,Z,a,invS,Jd,d,u2,w,ll);
    
    