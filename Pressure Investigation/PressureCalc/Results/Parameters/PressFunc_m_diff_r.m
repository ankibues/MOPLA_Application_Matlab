function PressFunc_m_diff_r(r,m)

% Internal Pressure field for an isotropic ellipse in Linear Anisotropic Material, with Pressure
% calcutions inside the ellipsoid.
 [Alp,ww] = Gauss1(8000);
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 ;0 -1];  
  
   
%  two viscosities for anisotropic medium
   %m =10;
   Nn=1;
   

  
etaE = r* Nn;      
   
%  ------------------------------------------------------------------------   
   D = 0.5 * (L + L');
   W = 0.5 * (L - L'); 
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   
    sqrtm= (m)^0.5;
 %  length of the ellipse axises
   a= [5;1];                        
   R= 5;
 Pdev_steps_alpha1= zeros(180,1);
 Pdev_steps_alpha2= zeros(180,1);
 lambda_check= zeros(2,2,180);
 S_check= zeros(2,2,2,2,180);
  lambdainva= zeros(180,1);
for i=0:180  
   ang = (pi/180)*i;
   
   q = [cos(ang), sin(ang); -sin(ang), cos(ang)];

   
   
% decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D_bar  = q * D * q';
   W_bar  = q * W * q';
    
   
   % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  
Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + ((sin(2*ang))^2)/m);
Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + ((cos(2*ang))^2)/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)= Cm(1,1,1,2);
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
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Ja,ang,Nn,m,Cm,Alp,ww);
 lambda_check(:,:,1+i)= LAMBDA; 
 lambdainva(i+1,1) = inva2D(LAMBDA);
 S_check(:,:,:,:,1+i)= S; 
% Interior fields (e,w,p)
  
  
% describe D,W,C in the clast's coordinate system 

 
  Ce    = Transform2D(Ce, q);  
  S_bar  = Multiply2D(Cm, D_bar);
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization 
   
% e: interior field strain-rate
  invS = FourTensorInv2D(S);
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
% p:  field pressure difference
  b1      = R_Multiply2D(LAMBDA, Cm);
  b2      = contract1(R_Multiply2D(b1, invS), de);
  p       = b2 ; 

 
 theta= t2;
 
 
 Sigma_ellipse = Multiply2D(Ce,e);
 
 Pdev_steps_alpha1(i+1,1)= p/SI;            % normalized with respect to far field stress invariant
 
 
end

x=0:180;
if r<1
plot(x,Pdev_steps_alpha1,':')
else
plot(x,Pdev_steps_alpha1)
end

hold on
xlabel('Angle in degrees');



end





