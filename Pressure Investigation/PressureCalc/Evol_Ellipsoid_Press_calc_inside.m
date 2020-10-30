function Evol_Ellipsoid_Press_calc_inside(pp,ee,gg)
% Motion of single ellipsoid in Power Law Isotropic Material, with Pressure
% calcutions inside the ellipsoid.
% This shows the pressure variation inside an ellipsoid as it is deformed.


%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [ee gg 0;0 -ee 0;0 0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 500;
   
%  Power law coefficients for matrix and clast
   Nm=3;
   Nc=3; 
   r= 100;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% for stress invariant for normalization :
   Sigma = 2* Multiply(Jd,D);                 % far field stress value( sigma= 2*Jd:E)i.e Viscosity of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   strainInva = inva(D);
   Cm    = 2*Jd;
% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                      % initial length of the ellipsoid axis
   a= [10;5;1];
   R=5;
   ang = [0;0;pi*pp/180] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
% calculations to follow:
 Pdev_steps_alpha1= zeros(steps,1);
 Pdev_steps_alpha2= zeros(steps,1);
 Pdev_steps_analytical= zeros(steps,1);
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  strainininvacheck= zeros(1,steps);
  sigmaininvacheck= zeros(1,steps);
  r_new_check= zeros(steps,1);
  SI_check = zeros(steps,1);
for i=1:steps
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
    h= Contract(2*Jd/Nm, invS);           % Here is the Cinverse, divided by Nm             
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
    
    qanew= QL(a,q,theta,dE,tincr);       
    
    %  write updated Q to Q_evl        
        Q_evl(:,:,i)=qanew(1:3,1:3);
    %  write updated a to A_evl
        A_evl(:,i)= qanew(:,4);
        R_new= A_evl(2,i)/A_evl(3,i);
    Pdev_in= contract1(p_in,u2);
    Sigma_ellipsoid = 2* Multiply(Jd,dE)* r_new;
    Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);
    sigmaininvacheck(1,i) = sigmaininva;   
    strainininvacheck(1,i)= strainininva;
    Pdev_analytical = 2*d(2,2)*(1-r_new)*((R_new^2)-1)/((R_new^2) +2*r_new*R_new+1);
    Pdev_steps_analytical(i,1)= Pdev_analytical/Sigma_Inva;
    Pdev_steps_alpha1(i,1)= Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
    Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipsoid ; 
    a=qanew(:,4);
    q=qanew(1:3,1:3);
end 
%     compute two spherical angles for three axes
    %  [a1_ang, a2_ang, a3_ang] = ConvertQ2Angs(Q_evl);
gamma= tincr*strainInva*ones(1,steps);
strain= cumsum(gamma);
%Equal_area_projection(a1_ang, a2_ang, a3_ang);
%figure
%x=1:numel(strain);
plot(strain,Pdev_steps_alpha1,'g')
hold on
%plot(x, tincr*strainininvacheck, 'r')
xlabel('\gamma');
ylabel('P');
end

%[m1,i1]= max(Pdev_steps);
%[m2,i2]=min(Pdev_steps);
%Q_max= Q_evl(:,:,i1);
%Q_min= Q_evl(:,:,i2);

%[MAXang1, MAXang2, MAXang3] = ConvertQ2Angs2(Q_max);
%[Minang1, Minang2, Minang3] = ConvertQ2Angs2(Q_min);
%MAXang1= rad2deg(MAXang1);
%MAXang2= rad2deg(MAXang2);
%MAXang3= rad2deg(MAXang3);
%Minang1= rad2deg(Minang1);
%Minang2= rad2deg(Minang2);
%Minang3= rad2deg(Minang3);
%MAXangs= [MAXang1, MAXang2, MAXang3];
%Minangs= [Minang1, Minang2, Minang3];