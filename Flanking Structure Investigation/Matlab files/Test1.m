% This file is to test the velocity field calculation
% Putting r=1, the velocity fields are similar as expected in homogeneous field.

% Test of velocity field calculation
L = [0 -1 0;0 0 0;0 0 0];
 D = 0.5 * (L + L');
 W = 0.5 * (L - L');  
 [Jd, ~, ~, ~, b] = Jnb();
  a=[1000;3;1];
% 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = -5:.5:5;           % grid vector: y'axis,a2
        zgv     = -5:.5:5;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
   
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        epp    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(epp);
      
   

   vel_f= zeros(2,1,n);
   
    YY= squeeze(Y);
    ZZ= squeeze(Z);
    [nn,~]= size(squeeze(Z));
   
   
    
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
 
   pp=0*pi/180;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;
   Nm=1;
   Nc=1; 
   r= 1;
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q';
    invS= Inverse(S,b);
    h= Contract(2*Jd/Nm, invS);           % Here is the Cinverse, divided by Nm             
    p_in = R_Multiply(p,h);
    E= Ed(Nm,Nc,r,S,d,epsilonII,b,Jd);
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;                                     
    ll= dE + we;
    
    
  %  
    
    parfor ir=1:n
        vel_f(:,:,ir)= Velocity_field_calc4(epp(2,ir),epp(3,ir),a,invS,Jd,d,u2,w,ll);
        disp(ir);
    end
   
   
    VV = zeros(2,1,nn);
    VV(:,:,ind)   = vel_f;
    VV(:,:,~ind)=0;
    V = reshape(VV,2,1,nn,nn);

    
    u= squeeze(V(1,1,:,:));
    v=squeeze(V(2,1,:,:));
    figure
    quiver(YY,ZZ,u,v)
    pbaspect([1 1 1])
    xlim([-4 4])
    ylim([-4 4])
   %}