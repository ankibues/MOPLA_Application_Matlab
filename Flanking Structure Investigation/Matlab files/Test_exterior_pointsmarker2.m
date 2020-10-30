% Deformation of marker with exterior points around inclusion- modelling a
% flanking structure
tic
% Note : In this code, mesh is defined with respect to inclusion, but coordinates will be converted back to far field flow for plotting. 

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 1 0;0 0 0;0 0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 1;
   
%  Power law coefficients for matrix and inclusion
   Nm=1;
   Nc=1; 
   r= 100000;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
% Stiffness tensor of the matrix   
   Cm    = 2*Jd;
   
% initial state of the ellipsoid                                    
   a= [100;2;1.999];                                                     % initial length of the ellipsoid axis
   pp= pi/6;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
 
 % 3D meshgrid-------- this meshgrid is defined in ellipsoid coordinate
 % system.
 
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =  -10:.05:10;           % grid vector: y'axis,a2
        zgv     = -10:.05:10;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);   

 % calculation  of points of the marker-------------------------------------     
        
    % size of the marker ( Here, I consider an elliptical marker, whose shape can be varied)
    a1= [100;3;3.0000001];
    a2= [100;4;4.0000001];
 %------------------------------This is defined in ellipsoid coordinate system.-------------   
        % Exterior points of marker to be followed these are the points
        % corresponding to the marker     
      
        ind1 = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2;          
        ind2 = (X./a1(1)).^2 + (Y./a1(2)).^2 + (Z./a1(3)).^2;
        ind3 = (X./a2(1)).^2 + (Y./a2(2)).^2 + (Z./a2(3)).^2;
        x_ex  = X(ind1>1 & ind2<1);
        y_ex  = Y(ind1>1 & ind2<1);
        z_ex  = Z(ind1>1 & ind2<1);
        
        x_exxx  = X(ind2>1);
        y_exxx  = Y(ind2>1);
        z_exxx  = Z(ind2>1);
        
        epp    = cat(1,x_ex',y_ex',z_ex');   % So these are the points to be followed at each time step.
        epp2    = cat(1,x_exxx',y_exxx',z_exxx'); % other marker to follow !
   
% calculations to follow:
  
  F(steps)= struct('cdata',[],'colormap',[]);
  
  
  
for i=1:steps
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
    
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
     % NOTE: External field has to be calculated in every time step since
     % the inclusion is deforming and rotating.
   %-----------------------------------------------------------------------------------------------------------------------------
   %-----------------------------------------------------------------------------------------------------------------------------
   %---------------------THIS PROCEDURE IS TO CALCULATE EXTERNAL AND INTERNAL VELOCITY GRADIENT FIELD ----------------------------
     % Exterior points
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        %ind   = (X./1).^2 + (((Y).*cos(pp)+ Z.*sin(pp))./a(2)).^2 + (((Y).*sin(pp)-Z.*cos(pp))/a(3)).^2 > 1;
        x_exx  = X(ind);
        y_exx  = Y(ind);
        z_exx  = Z(ind);
        ep    = cat(1,x_exx',y_exx',z_exx');
        [~,n] = size(ep);
        num   = numel(X);
        
        
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
        parfor rr=1:n
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
        end
        
        LL        = zeros(3,3,num);
        lll= ones(size(LL(:,:,~ind)));   % this is some matrix manipulation to 
        kkk= ll.*lll(:,:,:);
        LL(:,:,~ind)  = kkk;
        LL(:,:,ind)   = L_Ex;
        [nn,~]= size(squeeze(Z));
        LL       = reshape(LL,3,3,nn,nn); 
        YY= squeeze(Y);
        ZZ= squeeze(Z);
   %-----------------------------------------------------------------------------------------------------------------------------
   %-----------------------------------------------------------------------------------------------------------------------------
   %---------------------------------------------------------------------------------------------------------------------------------     
    
    % Now here, use the algorithm to calculate velocity field from the
    % velocity gradient field.
       Vel_field= VelGradtoVelfield(a,YY,ZZ,nn,ll,LL);   % this function gives me the velocity field in the mesh
      
    % for 1st marker   
       EP= [epp(2,:);epp(3,:)];
       vel_f = Velocity_field(EP,YY,ZZ,Vel_field);     % this function gives me the velocity fields at the marker points.
       disp_f= vel_f*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form
       [~,N]= size(EP);
       kkk= zeros(1,N);
       DISP_F= [kkk; squeeze(disp_f)];
       newpos_plot= zeros(size(DISP_F));
       newpos= epp + DISP_F;  % this is the new position calculated for the marker points.
                   
        for rrr=1:N                                % bringing to general coordinate system
           newpos_plot(:,rrr)= q'*newpos(:,rrr);
        end
     % for 2nd marker   
       EP2= [epp2(2,:);epp2(3,:)];
       vel_f2 = Velocity_field(EP2,YY,ZZ,Vel_field);     % this function gives me the velocity fields at the marker points.
       disp_f2= vel_f2*tincr; % this corresponds to 2-D coordinates/ has to converted to 3D coordinate form
       [~,N2]= size(EP2);
       kkk2= zeros(1,N2);
       DISP_F2= [kkk2; squeeze(disp_f2)];
       newpos_plot2= zeros(size(DISP_F2));
       newpos2= epp2 + DISP_F2;  % this is the new position calculated for the marker points.
                   
        for rrr=1:N                                % bringing to general coordinate system
           newpos_plot(:,rrr)= q'*newpos(:,rrr);
        end
            
        for rrr=1:N2                                % bringing to general coordinate system
           newpos_plot2(:,rrr)= q'*newpos2(:,rrr);
        end
        
       % For interior points of inclusion
         inddd   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 < 1;
        x_in  = X(inddd);
        y_in  = Y(inddd);
        z_in  = Z(inddd);
        EPP    = cat(1,x_in',y_in',z_in'); 
        %----------------these points are in inclusion coordinate system/so
        %has to be converted into global coordinate system.
        [~,NN]= size(EPP);
        for rrr=1:NN                                % bringing to general coordinate system
           EPP(:,rrr)= q'*EPP(:,rrr);
        end
       % plotting part for each step
    
       XX= newpos_plot(1,:);
       YY= newpos_plot(2,:);
       XX2= newpos_plot2(1,:);
       YY2= newpos_plot2(2,:);
       XXX= EPP(1,:);
       YYY= EPP(2,:);
        
        plot(XX,YY,'.','MarkerSize',5)     % this is for the marker !
        hold on
        plot(XX2,YY2,'.','MarkerSize',5)     % this is for the second marker !
        hold on
        plot(XXX,YYY,'.','MarkerSize',5)    % this is for the inclusion !
        xlim([-5 5])
        ylim([-5 5])
        F(i)= getframe;
        clf;
      
        
             
        %updating for next step  
        qanew= QL(a,q,theta,dE,tincr);  
         a=qanew(:,4);
         q=qanew(1:3,1:3);
         epp= newpos;
         epp2= newpos2;
        
end 
 %{   
     
v= VideoWriter('vid1-circular-rigid_2markers.avi');
v.FrameRate =5;
open(v)
writeVideo(v,F)
close(v)

 %}       
        
  toc
        