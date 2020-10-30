 % Internal Pressure field for an isotropic ellipse in Linear Anisotropic Material, with Pressure
% calcutions inside the ellipsoid.

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 ;0 -1]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 1000;
%  number of computation steps between output sets 
   mm=10; 
%  two viscosities for anisotropic medium
   m = 25;                              % strength of anisotropy
   Nn =  1;

%  Viscousity of isotropic inclusion(for r=10)
   etaE= 10* Nn;           
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   


% initial state of the ellipse
                                     

   a= [5;2];                                             % initial length of the ellipse axis

   ang = pi/4 ;                                          % initial angle of the ellipse axis
 
   q = [cos(ang), sin(ang); -sin(ang), cos(ang)];
      
 Pdev_steps_alpha1= zeros(steps,1);
 Pdev_steps_alpha2= zeros(steps,1);
 Q_evl= zeros(2,2,steps);
 A_evl= zeros(2,steps);
 strainininvacheck = zeros(1,steps);
 SICheck = zeros(1,steps);
 angcheck = zeros(1,steps);
 Qcheck = zeros(2,2,steps);
 sigmaininvacheck = zeros(1,steps);
 
 Cm    = zeros(2,2,2,2);
 
for i=1:steps  
   
   % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + (sin(2*ang))^2/m); Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + (cos(2*ang))^2/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)=Cm(1,1,1,2);
Cm(2,1,1,1)= Cm(1,2,1,1);
Cm(1,1,2,1)= Cm(1,1,1,2);
Cm(2,1,2,2)= Cm(1,2,2,2);
Cm(2,2,2,1)= Cm(2,2,1,2);
Cm(2,1,1,2)= Cm(1,2,1,2);
Cm(2,1,2,1)= Cm(2,1,1,2);
Cm(1,2,2,1)= Cm(2,1,2,1);
  
  
  

  
% stiffness tensor of the ellipsoid, 
  Ce    = 2*etaE*Jd; 


% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S_el,PI_el) for Interior points  
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Ja,ang,Nn,m,Cm);

% Interior fields (e,w,p)

% describe D,W,C in the clast's coordinate system 
  D_bar  = q * D * q';
  W_bar  = q * W * q';
  

  Ce    = Transform2D(Ce, q);  
  S_bar  = Multiply2D(Cm, D_bar);
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization               
  SICheck(1,i)= SI;
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
 qanew= QL2D(a,q,theta,e,tincr);   
 %  write updated Q to Q_evl        
 Q_evl(:,:,i)=qanew(1:2,1:2);
 %  write updated a to A_evl
 A_evl(:,i)= qanew(:,3);
  
 Pdev_in= p;
 Sigma_ellipse = Multiply2D(Ce,e);
 Sigma_Inva_ellipse = inva2D(Sigma_ellipse);
 Pdev_steps_alpha1(i,1)= Pdev_in/SI;            % normalized with respect to far field stress invariant
 Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipse ; 
 a=qanew(:,3);
 q=qanew(1:2,1:2); 
 Qcheck(:,:,i)= q;
 angcheck(1,i) = ang;
 ang= atan(q(1,2)/q(1,1));
end
x=1:steps;
figure
plot(x,Pdev_steps_alpha1,'g')
hold on
xlabel('deformation steps');
ylabel('alpha1');
   
