 % Internal Pressure field for an isotropic ellipse in Linear Anisotropic Material, with Pressure
% calcutions inside the ellipsoid.

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 2 ;0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 1000;
%  number of computation steps between output sets 
   mm=10; 
%  Viscousity of isotropic inclusion
   r= 10;
%  two viscosities for anisotropic medium
   m = 25;                              % strength of anisotropy
   Ns =  1;
   Nn = m*Ns;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   


% initial state of the ellipse
                                     

   a= [5;1];                                             % initial length of the ellipse axis

   ang = 0 ;                                          % initial angle of the ellipse axis
 
   
   
   
 Pdev_steps_alpha1= zeros(180,1);
 Pdev_steps_alpha2= zeros(180,1);
 SIcheck= zeros(180,1);
for i=0:180  
   ang = (pi/180)*i;
   
   q = [cos(ang), sin(ang); -sin(ang), cos(ang)];
   % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  
  Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + (sin(2*ang))^2/m);
  Cm(2,2,2,2)=  Cm(1,1,1,1);
  Cm(1,1,2,2)= -Cm(1,1,1,1);
  Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
  Cm(1,2,2,2)= -Cm(1,1,1,2);
  Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + (cos(2*ang))^2/m);
  
% stiffness tensor of the ellipsoid, 
  Ce    = 2*r*Jd; 


% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S_el,PI_el) for Interior points  
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Js,Ja,ang,Nn,m,Cm);

% Interior fields (e,w,p)

% describe D,W,C in the clast's coordinate system 
  D_bar  = q * D * q';
  W_bar  = q * W * q';
   
  Ce    = Transform2D(Ce, q);  
  S_bar  = Multiply2D(Cm, D_bar);
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization 
  P      = SI;                               
  SIcheck(i+1,1)= SI;
% e: interior field strain-rate
  invS= FourTensorInv2D(S);
  H       = Contract2D(Cm,(invS - Jd));   %Hill's constraint tensor
  t1      = FourTensorInv2D(H + Ce);
  A       = Contract2D(t1,(H + Cm));      %strain-rate partitioning tensor
  e       = Multiply2D(A, D_bar);
  
  
% s: interior field deviatoric-stress
  s       = Multiply2D(Ce, e);
% w: interior field vorticity
  de      = e - D_bar;                  
  t2      = Multiply2D(Contract2D(PI, invS), de);
  w       = t2 + W_bar;
% p: interior field pressure
  b1      = R_Multiply2D(LAMBDA, Cm);
  b2      = contract1(R_Multiply2D(b1, invS), de);
  p       = b2 ; %+ P;

  
 wE = Wd2D(a,W_bar,e);
 theta= w-wE;
    
 
  
 Pdev_in= p;
 Sigma_ellipse = Multiply2D(Ce,e);
 Sigma_Inva_ellipse = inva2D(Sigma_ellipse);
 Pdev_steps_alpha1(i+1,1)= Pdev_in/SI;            % normalized with respect to far field stress invariant
 Pdev_steps_alpha2(i+1,1) = Pdev_in/Sigma_Inva_ellipse ; 
 
end

   
