
% Internal Pressure field for an isotropic ellipse in Linear Anisotropic Material, with Pressure
% calcutions inside the ellipsoid.
 [Alp,ww] = Gauss1(8000);
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 ;0 -1];  
   ep=1;
   v= 0;
   
%  two viscosities for anisotropic medium
   m =10;
   Nn=1;
   
%  Viscousity of isotropic inclusion(for r=10)
r = .10;   
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
 e11= zeros(180,1);
 e12= zeros(180,1);
 w12= zeros(180,1);
  e11_anal= zeros(180,1);
   e12_anal= zeros(180,1);
   w12_anal= zeros(180,1);
    pdiff_anal= zeros(180,1);
 angcheck= zeros(180,1);
 lambdainva= zeros(180,1);
for i=0:180  
   ang = (pi/180)*i;
   angcheck(i+1,1)=ang;
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
  
% analytical results from Fletcher 2009  (as in Jiang 2016)
  
   E11= ep*cos(2*ang) + v*sin(2*ang)/2;
   E11m= ep*cos(2*ang) + (v/m)*sin(2*ang)/2;
   E12= -ep*sin(2*ang) + v*cos(2*ang)/2;
   E12m= -ep*sin(2*ang) + (v/m)*cos(2*ang)/2;
   W12= v/2;
   E11bar= -r*E11 + E11m;
   E12bar= -r*E12 + E12m;
   e11analytic = 2*sqrtm*R*E11bar/ (R^2 + 2*r*sqrtm*R + 1);
   e12analytic = sqrtm*(R^2 + 1)* E12bar/(r*sqrtm*(R^2 +1) + 2*R);
   w12analytic = sqrtm*(R^2 - 1)* E12bar/(r*sqrtm*(R^2 +1) + 2*R);
   pdiff_analytic = 2*Nn*(R^2 - 1)*(E11bar)/(R^2 + 2*r*sqrtm*R + 1);
   E_anal= [E11, E12; E12, -E11];
  
% describe D,W,C in the clast's coordinate system 

 
  Ce    = Transform2D(Ce, q);  
  S_bar  = Multiply2D(Cm, D_bar);
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization 
  S_anal=  Multiply2D(Cm, E_anal);                            
  SI_anal= inva2D(S_anal);                      % Analytical stress invariant for far-field for normalization 
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

 e11(i+1,1)= de(1,1)/SI; 
 e12(i+1,1)= de(1,2)/SI;
 
 theta= t2;
 w12(i+1,1)= theta(1,2)/SI;
 
 %analytical results update 
 e11_anal(i+1,1)= e11analytic/SI_anal;
   e12_anal(i+1,1)= e12analytic/SI_anal;
   w12_anal(i+1,1)= w12analytic/SI_anal;
   pdiff_anal(i+1)= pdiff_analytic/SI_anal;
  

 Sigma_ellipse = Multiply2D(Ce,e);
 Sigma_Inva_ellipse = inva2D(Sigma_ellipse);
 Pdev_steps_alpha1(i+1,1)= p/SI;            % normalized with respect to far field stress invariant
 Pdev_steps_alpha2(i+1,1) = p/Sigma_Inva_ellipse ; 
 
end

x=0:180;
figure
plot(x,Pdev_steps_alpha1,'g')
hold on
plot (x,pdiff_anal,'--')
xlabel('Angle in degrees');
ylabel('alpha1 analytical vs empirical');

figure 
plot(x,w12,'g')
hold on 
plot(x,w12_anal,'b--')
xlabel('Angle in degrees');
ylabel('w12 analytical vs empirical');

figure
plot(x,e11,'g')
hold on
plot(x,e12,'b')
hold on 
plot(x,e11_anal,'b--')
hold on
plot(x,e12_anal,'g--')
hold on 
xlabel('Angle in degrees');
ylabel('e11 and e12 analytical vs empirical');





