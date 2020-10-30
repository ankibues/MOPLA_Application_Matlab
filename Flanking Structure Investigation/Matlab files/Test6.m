%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 1 0;0 -1 0;0 0 0]; 
   
   
%  Power law coefficients for matrix and clast
   Nm=1;
     % since linear isotropic case is being considered
   r= 100;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm, b] = Jnb(); 
   
% for stress invariant for normalization :
   Sigma = 2* Multiply(Jd,D);                 % far field stress value( sigma= 2*Jd:E)i.e Viscosity of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);


% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                      % initial length of the ellipsoid axis
   a= [6;5;1];
   
   ang = [0; 0;(pi/180)*0];                       % pressure variations for 135 degree angle
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant



    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q';                                      % w is the vorticity of material, in ellipsoid's frame of reference 
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
    h= Contract(2*Jd, invS);                        
    p_in = R_Multiply(p,h);
    

    dE       = Multiply(fdE(Nm,r,S,b,Jd),d);                      
                                      % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    ll= dE + we; 
%function V = Velocity_field_calc5(X,Y,a,invS,Jd,d,u2,w,ll)
% This function calculates velocity field at any coordinate point outside the inclusion
% Input: X,Y : coordinates of the point; a: 3*1 matrix representing 

% This function will take the input as the coordinates of the marker point,and gives back the velocity field at that point.
% Once we get the coordinate, we decide how it should be divided into the
% steps.

% Instead of using the cartesian coordinates, Here, I use Elliptical
% coordinates for efficient integration of velocity gradient field.
% Note: The actual calculation is done in Cartesian coordinates; Elliptic
% coordinates are only used to get the shortest path of integrating
% velocity gradient field. Method is more accurate and efficient.
X=10;
Y=10;
A= a(2:3);
[epsilon,eta]= CartesiantoElliptic(X,Y,A);     % bringing the points to Elliptic coordinates

% choosing the known points for velocity field
% these known points in elliptic coordinates
[epsilon0,~]= CartesiantoElliptic(A(1),0,A);    
eta0= real(eta);
% these known points in cartesian coordinates, so that velocity field can
% be calculated.
[Xo,Yo]= ElliptictoCartesian(epsilon0,eta0,A);

Vo= ll*[0;Xo;Yo];                % reference velocity field at the boundary


%{
Later On...Using the Lambda approach, we can give variable lengths of
 intervals. For closer points, more intervals, while for farther points,
less intervals at farther distances.
%}
lambda = solve_eq(a,[0;X;Y]);

% Lambda_External 
         if lambda<= .6920
             kk= .0001;
         elseif lambda>.64 && lambda<= 1.80
             kk= .001;
         else 
             kk= .01;            
         end
%}
%kk=.0001;

EPS= epsilon0:kk:epsilon;
ETA= eta0*ones(size(EPS)); % these are interval coordinates in elliptic coordinate system

% bringing them to cartesian system

aa= A(1);  % long semi axis
b= A(2);  % short semi axis
c= (aa^2 - b^2)^(1/2);

x= c*cosh(EPS).*cos(ETA);
y= c*sinh(EPS).*sin(ETA);

% choosing the intervals for integration; this will be along eta axis only.

ep= [zeros(size(x)); x;y];

[~,N1]=size(ep);
  
    
    
    
    
  % Exterior Fields  along this path
   %--------------------------------------------------------------------------
     % G
          G1 = Ex_Gtensor(a,ep);
        % LAMBDA_Ex 
        LAMBDA_Ex1   = zeros(3,3,N1);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex1(ii,j,:) = -1/3.*(G1(ii,j,1,1,:)+ G1(ii,j,2,2,:)+ G1(ii,j,3,3,:));
                LAMBDA_Ex1(j,ii,:) = LAMBDA_Ex1(ii,j,:);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G1(1,1,iii,iii,:)+ G1(2,2,iii,iii,:)+ G1(3,3,iii,iii,:));
            LAMBDA_Ex1(iii,iii,:) = squeeze(LAMBDA_Ex1(iii,iii,:))+t;
        end
        % S_Ex & PI_Ex
        S_Ex1   = zeros(3,3,3,3,N1);
        PI_Ex1  = zeros(3,3,3,3,N1);
        
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex1(ii,j,k,l,:)  = squeeze(G1(ii,j,k,l,:))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex1(ii,j,:));
                        S_Ex1(j,ii,k,l,:)  = S_Ex1(ii,j,k,l,:);
                        S_Ex1(j,ii,l,k,:)  = S_Ex1(ii,j,k,l,:);
                        S_Ex1(ii,j,l,k,:)  = S_Ex1(ii,j,k,l,:);
                        %PI_Ex
                        PI_Ex1(ii,j,k,l,:) = 1/2.*(delt(j,k).*LAMBDA_Ex1(ii,l,:) +...
                                            delt(j,l).*LAMBDA_Ex1(ii,k,:) - delt(ii,k)...
                                            .*LAMBDA_Ex1(j,l,:) - delt(ii,l).*...
                                            LAMBDA_Ex1(j,k,:));
                        PI_Ex1(ii,j,l,k,:) = PI_Ex1(ii,j,k,l,:);
                        PI_Ex1(j,ii,k,l,:) = -PI_Ex1(ii,j,k,l,:);
                        PI_Ex1(j,ii,l,k,:) = -PI_Ex1(ii,j,k,l,:);
                    end
                end
            end
        end
        S_Ex1(1,1,1,1,:) = -(S_Ex1(1,1,2,2,:)+S_Ex1(1,1,3,3,:));
        S_Ex1(2,2,2,2,:) = -(S_Ex1(2,2,1,1,:)+S_Ex1(2,2,3,3,:));
        S_Ex1(3,3,3,3,:) = -(S_Ex1(3,3,1,1,:)+S_Ex1(3,3,2,2,:));
  
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex1    = zeros(3,3,N1);
        s_Ex1    = zeros(3,3,N1);
        s_Exinva1= zeros(1,N1);
        w_Ex1    = zeros(3,3,N1);
        L_Ex1    = zeros(3,3,N1);
        for rr=1:N1
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex1(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex1(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex1(:,:,rr) = Multiply(2*Jd, e_Ex1(:,:,rr));
            s_Exinva1(:,rr)= inva(s_Ex1(:,:,rr));
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex1(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex1(:,:,rr) = v4 + w;
        % L_Ex: exterior field of velocity gradient tensor
            L_Ex1(:,:,rr)     = e_Ex1(:,:,rr)+ w_Ex1(:,:,rr)  ;
        end
        

        
%------------------------------------------------------------------------------------------------------------------------
% Removing one term from the series of Ls calculated from all points.
    L_Ex1(:,:,N1)= [];
   
  ep1= [zeros(3,1), ep(:,1:N1-1)]; 
  % This process is done to calculate the delta x and delta y values for
  % each intervals. Note: the intervals are equal in Elliptic form, but not
  % equal in Cartesian form. So, this has to be done.
  delta= ep-ep1;
  delta(:,1)=[];
  deltaX= delta(2,:);
  deltaY= delta(3,:);
  
  V1= Vo(2) + sum(deltaX.*(squeeze(L_Ex1(2,2,:)))',2) + sum(deltaY.*(squeeze(L_Ex1(2,3,:)))',2);
  V2= Vo(3) + sum(deltaX.*(squeeze(L_Ex1(3,2,:)))',2) + sum(deltaY.*(squeeze(L_Ex1(3,3,:)))',2);
  
  
  
    
    V= [V1;V2];




%end
