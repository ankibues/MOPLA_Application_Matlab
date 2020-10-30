% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure

% Note : In this code, mesh is defined with respect to inclusion. To model
% the flanking structures in a way to help understand the bulk flow, the
% mesh should be in global coordinate system !

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 1 0;0 0 0;0 0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 100;
   
%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1; 
   r= .001;
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
   a= [20;8;.1];                                                     % initial length of the ellipsoid axis
   pp= pi/3;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
   
 % calculation  of points of the marker----------------------------------------  
      % 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =  -5:.05:5;           % grid vector: y'axis,a2
        zgv     = -5:.05:5;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
       
    % size of the marker ( Here, I consider an elliptical marker, whose shape can be varied)
    a1= [30;.4;30];
    
        % Exterior points of marker to be followed
        ppp=0; %this is because this is for the marker.
        ind1   = (X./1).^2 + (((Y).*cos(ppp)+ Z.*sin(ppp))./a1(2)).^2 + (((Y).*sin(ppp)-Z.*cos(ppp))/a1(3)).^2;
        ind2   = (X./1).^2 + (((Y).*cos(ppp)+ Z.*sin(ppp))./a(2)).^2 + (((Y).*sin(ppp)-Z.*cos(ppp))/a(3)).^2; % this corresponds to the inclusion
        x_ex  = X(ind1<1 & ind2>1);
        y_ex  = Y(ind1<1 & ind2>1);
        z_ex  = Z(ind1<1 & ind2>1);
        ep    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(ep);
        num   = numel(X);
       
        % total points number
        % exterior points number
        %{
%----------------------------------------------- assuming that the mesh is
in ellipsoid's coordinate system.
        epp = zeros(3,n);
        
        qq= [1, 0 ,0; 0 ,cos(pp), sin(pp); 0 ,-sin(pp) ,cos(pp)]; % this part transforms(rotates) the coordinates for their use in G calculation. 
        for kk=1:n
            epp(:,kk)= qq*ep(:,kk);
        end
        %}
        New_Pos= zeros(3,n);
        New_Pos_step= zeros(3,n,steps);
        %-------------------------------------------------------------------------- 
   
   
   
   
% calculations to follow:
 
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  F(steps)= struct('cdata',[],'colormap',[]);
  
  
for i=1:steps
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    
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
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
  
     %Exterior Fields  
     %--------------------------------------------------------------------------
     % G
          G = Ex_Gtensor(a,ep);
        % LAMBDA_Ex 
        LAMBDA_Ex   = zeros(3,3,n);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j,:) = -1/3.*(G(ii,j,1,1,:)+ G(ii,j,2,2,:)+ G(ii,j,3,3,:));
                LAMBDA_Ex(j,ii,:) = LAMBDA_Ex(ii,j,:);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G(1,1,iii,iii,:)+ G(2,2,iii,iii,:)+ G(3,3,iii,iii,:));
            LAMBDA_Ex(iii,iii,:) = squeeze(LAMBDA_Ex(iii,iii,:))+t;
        end
        % S_Ex & PI_Ex
        S_Ex   = zeros(3,3,3,3,n);
        PI_Ex  = zeros(3,3,3,3,n);
        
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex(ii,j,k,l,:)  = squeeze(G(ii,j,k,l,:))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex(ii,j,:));
                        S_Ex(j,ii,k,l,:)  = S_Ex(ii,j,k,l,:);
                        S_Ex(j,ii,l,k,:)  = S_Ex(ii,j,k,l,:);
                        S_Ex(ii,j,l,k,:)  = S_Ex(ii,j,k,l,:);
                        %PI_Ex
                        PI_Ex(ii,j,k,l,:) = 1/2.*(delt(j,k).*LAMBDA_Ex(ii,l,:) +...
                                            delt(j,l).*LAMBDA_Ex(ii,k,:) - delt(ii,k)...
                                            .*LAMBDA_Ex(j,l,:) - delt(ii,l).*...
                                            LAMBDA_Ex(j,k,:));
                        PI_Ex(ii,j,l,k,:) = PI_Ex(ii,j,k,l,:);
                        PI_Ex(j,ii,k,l,:) = -PI_Ex(ii,j,k,l,:);
                        PI_Ex(j,ii,l,k,:) = -PI_Ex(ii,j,k,l,:);
                    end
                end
            end
        end
        S_Ex(1,1,1,1,:) = -(S_Ex(1,1,2,2,:)+S_Ex(1,1,3,3,:));
        S_Ex(2,2,2,2,:) = -(S_Ex(2,2,1,1,:)+S_Ex(2,2,3,3,:));
        S_Ex(3,3,3,3,:) = -(S_Ex(3,3,1,1,:)+S_Ex(3,3,2,2,:));
  
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex    = zeros(3,3,n);
        s_Ex    = zeros(3,3,n);
        s_Exinva= zeros(1,n);
        w_Ex    = zeros(3,3,n);
        L_Ex    = zeros(3,3,n);
        F_Ex    = zeros(3,3,n);
        for rr=1:n
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex(:,:,rr) = Multiply(2*Jd, e_Ex(:,:,rr));
            s_Exinva(:,rr)= inva(s_Ex(:,:,rr));
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex(:,:,rr) = v4 + w;
        % L_Ex: exterior field of velocity gradient tensor
            L_Ex(:,:,rr)     = e_Ex(:,:,rr)+ w_Ex(:,:,rr)  ; 
        % F_Ex: position gradient tensor field
            F_Ex(:,:,rr)     = expm(L_Ex(:,:,rr)*tincr);
        % New Position vectors for the marker points
            New_Pos(:,rr)= F_Ex(:,:,rr)*ep(:,rr);    %
        
        end
        
    
    qanew= QL(a,q,theta,dE,tincr);       
    
    %  write updated Q to Q_evl        
        Q_evl(:,:,i)=qanew(1:3,1:3);
    %  write updated a to A_evl
        A_evl(:,i)= qanew(:,4);
        R_new= A_evl(2,i)/A_evl(3,i);
    
    %updating for next step  
    a=qanew(:,4);
    q=qanew(1:3,1:3);
    ep= New_Pos;
    New_Pos_step(:,:,i)= New_Pos;
    
    % plotting part for each step
    
       XX= ep(2,:);
       YY= ep(3,:);
        ind   = (X./1).^2 + (((Y).*cos(ppp)+ Z.*sin(ppp))./a(2)).^2 + (((Y).*sin(ppp)-Z.*cos(ppp))/a(3)).^2 < 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        epp    = cat(1,x_ex',y_ex',z_ex');
        XXX= epp(2,:);
        YYY= epp(3,:);
        
        plot(XX,YY,'.','MarkerSize',5)
        %fill(XX,YY,'r')
        hold on
        plot(XXX,YYY,'.','MarkerSize',5)
        %fill(XXX,YYY,'b')
        xlim([-20 20])
        ylim([-20 20])
        F(i)= getframe;
        clf;
        drawnow
end 
    
     
v= VideoWriter('vidtest21.avi');
v.FrameRate =5;
open(v)
writeVideo(v,F)
close(v)

        
        
        
        