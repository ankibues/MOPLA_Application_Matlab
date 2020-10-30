% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure

%Here, marker is defined by an outline only. That will be a quicker way to
%see how the shape of kinematic indicator comes up !
 
 Wk= .8;
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

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 50;
   
%  Power law coefficients for matrix and inclusion
   Nm=1;
   Nc=1; 
   r= 5;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a v
% initial state of the ellipsoid                                    
   a= [100;.4;.2];                                                     % initial length of the ellipsoid axis
   pp=0;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
 
 % 3D meshgrid-------- this meshgrid is defined in ellipsoid coordinateorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Stiffness tensor of the matrix   
   Cm    = 2*Jd;
   
 % system.
 
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =  -2:.03:2;           % grid vector: y'axis,a2
        zgv     =  -2:.03:2;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);   


 % calculation  of points of the marker-------------------------------------     
        
    % size of the marker ( Here, I consider an elliptical marker, whose shape can be varied)

     %------------------------------This is defined in ellipsoid coordinate system.-------------   
        % Exterior points of marker to be followed these are the points
        % corresponding to the marker     

      
  A=.6;
  B=.2;
   
C= (A^2 - B^2)^(1/2);
 
eta= 0:.03:2*pi;
epsilo= .2;       % this defines the distance of the marker from the ellipsoidal surface

x= C*cosh(epsilo)*cos(eta);
y= C*sinh(epsilo)*sin(eta);  
theta=0;
 [~,nn]=size(x);
 R= [cos(theta) -sin(theta); sin(theta) cos(theta)];
 k=zeros(2,nn);
 for i=1:nn
     k(:,i)= R*[x(i); y(i)];
 end
 x_r= k(1,:);
 y_r= k(2,:);
        
        x_ex  = zeros(size(x_r));
        y_ex  = x_r;
        z_ex  = y_r;
        epp    = [x_ex; y_ex; z_ex];   % So these are the points to be followed at each time step.
   
        [~,n]=size(epp);
% calculations to follow:
  
  F(steps)= struct('cdata',[],'colormap',[]);
  
  
  
for i=1:steps
   tic 
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    
    
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
    ll= dE + we;
   
    % Here, I use the algorithm to calculate velocity field from the
    % velocity gradient field of the marker points.
     % to check if any value is slightly less than 1(e.g. .999) :A bug in the
   % code
   K=(epp(2,:)./a(2)).^2 + (epp(3,:)./a(3)).^2;
   epp(:,K<1)=[];
    
    xxx= epp(2,:);
    yyy= epp(3,:);
    xxx(isnan(xxx)==1)=[];
    yyy(isnan(yyy)==1)=[];
     [~,n]=size(epp);
     vel_f= zeros(2,1,n);
    
     parfor ir=1:n
        vel_f(:,:,ir)= Velocity_field_calc4(xxx(ir),yyy(ir),a,invS,Jd,d,u2,w,ll); 
     end
    
       
       
    
    
       disp_f= vel_f*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form
       
       kkk= zeros(1,n);
       DISP_F= [kkk; squeeze(disp_f)];
       newpos_plot= zeros(size(DISP_F));
       newpos= epp + DISP_F;  % this is the new position calculated for the marker points.
                   
        parfor rrr=1:n                                % bringing to general coordinate system
           newpos_plot(:,rrr)= q'*newpos(:,rrr);
        end
        
        % For interior points of inclusion
       % 3D meshgrid in internal coordinate -----interior points of
        % inclusion-------------------------------------------------------------------------
        
        xgv1     = 0; %-5:.05:5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv1    = -1:.03:1;           % grid vector: y'axis,a2
        zgv1     = -1:.03:1;  % grid vector: z'axis,a3
        [X1,Y1,Z1] = meshgrid(xgv1,ygv1,zgv1);
        ind1= (X1./a(1)).^2 + (Y1./a(2)).^2 + (Z1./a(3)).^2;
        x_ex1  = X1(ind1<1);
        y_ex1  = Y1(ind1<1);
        z_ex1  = Z1(ind1<1);
        epp    = cat(1,x_ex1',y_ex1',z_ex1'); % These are points corresponding the inclusion
        [~,NN]= size(epp);
        EPP= zeros(size(epp));
        for rrr=1:NN                                % bringing to general coordinate system
           EPP(:,rrr)= q'*epp(:,rrr);    % in exterior coordinate system for further conversion
        end
        
        qanew= QL(a,q,theta,dE,tincr);  
          a=qanew(:,4);
         q=qanew(1:3,1:3); 
        
       
        
        % plotting part for each step
    
       XX= newpos_plot(1,:);
       YY= newpos_plot(2,:);
       
       XXX= EPP(1,:);
       YYY= EPP(2,:);
        f= figure('Visible','off');
        plot(XX,YY,'.','MarkerSize',2)     % this is for the marker !
        hold on
        plot(XXX,YYY,'.','MarkerSize',2)    % this is for the inclusion !
        
        
        xlim([-3 3])
        ylim([-3 3])
        pbaspect([1 1 1])
        F(i)= getframe(f);
        clf;
       %updating for next step  
        parfor rrr=1:n   % bringing from general coordinate system to the new coordinate system of rotated ellipsoid
           epp(:,rrr)= q*newpos_plot(:,rrr);
        end
        
             
       
         
     disp('Step completed');
     disp(i);
     toc
end 
    

v= VideoWriter('test_Wk=.8_r=5_aspect_.4-.2_lambda_.6-.2-.2_ia_0_marker0.avi');
v.FrameRate =2;
open(v)
writeVideo(v,F)
close(v)

 %}       
        
        
        %}