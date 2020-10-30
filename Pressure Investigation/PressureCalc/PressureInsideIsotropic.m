function PressureInsideIsotropic(r)
% Pressure variation in an ellipsoid in Linear Isotropic Material, with Pressure
% calcutions inside the ellipsoid, for different instantatneous angles of orientation.
%

%  Input parameters--------------------------------------------------------
%  Matrix flow field (Velocity field for simple shear flow)
   L     = [1 0 0;0 -1 0;0 0 0]; 

%  Power law coefficients for matrix and clast( since Linear case is
%  considered, both coefficients as 1)

   Nm=1;
   Nc=1;
   
%  Ratio of ellipsoid's viscosity to the matrix viscosity   
   %r= 10;
   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity tensor W    

   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
   
  
%  generating 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 
   
% Calculating Far field Stress Invariant for normalization :
   Sigma = 2* Multiply(Jd,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   strainInva = inva(D);
   Cm    = 2*Jd;                     % stiffness tensor for the matrix
   
% initial state of the ellipsoid
                                     % initial length of the ellipsoid axis
   a= [6;3;2];                   % Point to note: as ellipsoidal long axis is reduced, results vary with respect 
                                     % to analytical 2-D result
 
                                     % initial angles of the ellipsoid axis
   
 % assigning variables for calculations to follow:
 Pdev_steps_alpha1= zeros(180,1);
 Pdev_steps_alpha2= zeros(180,1);
 
for i=1:180
    ang = [0;0;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
    % why is this showing anamolas behavior at this orientation: ang = [pi/2 ;pi/2;(pi/180)*i];
    % if a2 oriented along Z-axis, use ang=  theta, pi/2,0
    q = Q(ang);
    d= q*D*q';                                      % Transforming D and W to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);                    
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    
    
    
   % [S1,p1,PI1] = SnpIn(a,Jd,Js,Ja);% this method uses quadratures. So, since
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
   %**************************
  
    
    invS= Inverse(S,b);
    
   % RR= R_Multiply(p,Inverse((S-Jd),b)); % value of psi we use in paper
    %*********************
    h= Contract(Cm_a/Nm, invS);                        
    p_in = R_Multiply(p,h);
    dE       = Multiply(fdE(Nm,r,S,b,Jd),d);
   
    %E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
   % r_new= E(1,4);
   % r_new_check(i,1)= r_new;                         % new r values being stored to see the variation of viscousity ratio
   % dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    LL = dE + wE;
   
    theta= we-wE;
    sigmainside= 2*dE*r;
 
    strainininva = inva(dE);
    sigmaininva= inva(sigmainside);
    
    
    Pdev_in= contract1(p_in,u2);
    Sigma_ellipsoid = 2* Multiply(Jd,dE)* r;
    Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);

  
    

   
    %Pdev_steps_analytical(i,1)= Pdev_analytical/Sigma_Inva;
   % Pdev_steps_analyticalalpha2(i,1)= Pdev_analytical/Sigma_Inva_ellipsoid;
    Pdev_steps_alpha1(i,1)= Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
    Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipsoid ; 
    
  
  
    
    

end

x=1:180;
%figure
plot(x,Pdev_steps_alpha1,'--')
%xlabel('Angle in degress');
%ylabel('\alpha1');
hold on
%plot(x,Pdev_steps_analytical,'b--')
%%MaxAlpha2= max(Pdev_steps_alpha2);
%{
if r>1
    
    %plot(x,Pdev_steps_alpha2,'r')
    xlabel('Angle in degress');
    ylabel('\alpha2');
    hold on
   % plot(x,Pdev_steps_analyticalalpha2,'b--')
end
%}
end



