function viscosity_ratio_parameter2(a)

% for considering rod shaped element with longest axis parallel to Y-axis.
r= .0001:.001:1;

nn=numel(r);
 % assigning variables for calculations to follow:
  Pdev_steps_alpha1= zeros(nn,1);
  Pdev_steps_alpha2= zeros(nn,1);


for ii=1:nn
    
    
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
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Calculating Far field Stress Invariant for normalization :
   Sigma = 2* Multiply(Jd,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);
   
   Cm    = 2*Jd;                     % stiffness tensor for the matrix
   
% initial state of the ellipsoid
                                     % initial length of the ellipsoid axis
                                % Point to note: as ellipsoidal long axis is reduced, results vary with respect 
                                     % to analytical 2-D result
   R=5;
                                     % initial angles of the ellipsoid axis
   
  epsilonII=.5;                      % initial strain rate invariant( This is useful when considering non-linear case)

 
 
       i=0;                                 % angle in degrees corresponding to Max pressure
       ang = [0; pi/2;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
       q = Q(ang);
    
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
            h= Contract(2*Jd, invS);                        
            p_in = R_Multiply(p,h);
            E= Ed(Nm,Nc,r(ii),S,d,epsilonII,b,Jd);
          
                                  
            dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
            u1= Contract(PI, invS);
            u2= dE-d;
            we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
            wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
           
            
           
    
    
            Pdev_in= contract1(p_in,u2);
            Sigma_ellipsoid = 2* Multiply(Jd,dE)* r(ii);
            Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);
           
                   
            
            Pdev_steps_alpha1(ii,1)= Pdev_in/Sigma_Inva;  % normalized with respect to far field stress invariant
            Pdev_steps_alpha2(ii,1) = Pdev_in/Sigma_Inva_ellipsoid ; 
    
   
  
     
      
         
end
%{
plot(r,Pdev_steps_alpha1)
    xlabel('Viscosity Ratio');
    ylabel('\alpha1');

 %}
  
       
        
        %plot(r,Pdev_steps_alpha2)
        plot(r,Pdev_steps_alpha1)
        xlabel('Viscosity Ratio');
        ylabel('\alpha1');
       
       
    
end

  