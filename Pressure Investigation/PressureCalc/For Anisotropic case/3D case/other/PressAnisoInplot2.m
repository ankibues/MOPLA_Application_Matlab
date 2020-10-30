function PressAnisoInplot2(r)
%Input parameters--------------------------------------------------------
%  Matrix flow field (Velocity field for simple shear flow)
   L     = [1 0 0;0 -1 0;0 0 0]; 

%  Power law coefficients for matrix and clast( since Linear case is
%  considered, both coefficients as 1)

   
   
     gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
   
%  Ratio of ellipsoid's viscosity to the matrix viscosity ( r here is the viscosity of inclusion)  
   %r= .1;
   
   load('AnisotropicC_200inclusionsC','C_bar_evl');
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
   
%  generating 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   Ce=2*r*Jd;
% Calculating Far field Stress Invariant for normalization :
   Cm = C_bar_evl(:,:,:,:,200);      % stiffness tensor for the matrix
   
   Sigma = Multiply(Cm,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   strainInva = inva(D);
 % load('shearstrain10(boudinageat25).mat');
 %  Ce=  C_bar_evl(:,:,:,:,2);
% initial state of the ellipsoid
                                     % initial length of the ellipsoid axis
   a= [10;5;1];                   % Point to note: as ellipsoidal long axis is reduced, results vary with respect 
                                     % to analytical 2-D result
  
                                     % initial angles of the ellipsoid axis
   
 % assigning variables for calculations to follow:
 Pdev_steps_alpha1= zeros(180,1);
 Pdev_steps_alpha2= zeros(180,1);
  for i=1:180
    
    ang = [0; 0;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
    q = Q(ang);
    
    d= q*D*q';                                    % Transforming D and W to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    Carray  = C2OneDarray(Cm_a);
    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
    S = Contract(Jd, Contract(T,Cm_a));
    
    lambda = SnpInAniso3DPress(a,Cm_a);
    invS = Inverse((S-Jd),b);
    h= Contract(Cm_a, invS);                        
    p_in = R_Multiply(lambda,h);
    H = Contract(Cm_a,(invS-Jd));           % Eq.12a Jiang 2014
    A = Contract(FourTensorInv(H + Ce),(H + Cm_a));                     
    e= Multiply(A,d);    % e is strain rate in ellipsoid
    u1= e-d;                  
    sigmainside= Multiply(Ce,e);
    sigmaininva= inva(sigmainside);
    Pdev_in= contract1(p_in,u1);       
    Pdev_steps_alpha1(i,1)= Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
    Pdev_steps_alpha2(i,1) = Pdev_in/sigmaininva ; 
    
  end   
    
x=1:180;
%figure
%plot(x,Pdev_steps_alpha2,'r')
xlabel('Angle in degrees');
ylabel('P/\Sigma');
%hold on
plot(x,Pdev_steps_alpha1,'g')
xlabel('Angle in degrees');
ylabel('P/\Sigma');

end
