%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 1 0;0 0 0;0 0 0]; 
%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 10;

%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1; 
   r= 100;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
   Cm    = 2*Jd;
   

% initial state of the ellipsoid                                    
   a= [100;2;1.9999];                                                     % initial length of the ellipsoid axis
   pp= pi/6;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;

 % 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = -5:.05:5;           % grid vector: y'axis,a2
        zgv     = -5:.05:5;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
  
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);
         
        % Exterior points
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        
        %ind   = (X./1).^2 + (((Y).*cos(pp)+ Z.*sin(pp))./a(2)).^2 + (((Y).*sin(pp)-Z.*cos(pp))/a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(ep);
        num   = numel(X);  
   
   
   
   
for iji=1:steps    
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
    dE= E(1:3,1:3);                                  % dE is strain rate tensor in ellipsoid
    u1= Contract(PI, invS);
    u2= dE-d;
    we= Multiply(u1,u2)+ w ;                        % this is eq. 12b from Jiang 2016/ we is the vorticity of ellipsoid
    wE = Wd(a,w,dE);                                  % wE: shear spin part of rotation of ellipsoid.
    theta= we-wE;
    ll= dE + we;
 
        
        
         %Exterior Fields  
     %--------------------------------------------------------------------------
     % G
          G = Ex_Gtensor(a,ep);  
           [~,n] = size(ep);
           
          
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
        lll= ones(size(LL(:,:,~ind)));   % this is some matrix manipulation to make same kind of data structure.
        kkk= ll.*lll(:,:,:);
        LL(:,:,~ind)  = kkk;
        LL(:,:,ind)   = L_Ex;
        [nn,~]= size(squeeze(Z));
        LL       = reshape(LL,3,3,nn,nn); 
        YY= squeeze(Y);
        ZZ= squeeze(Z);
        Vel_field= zeros(2,1,nn,nn);
       
        
        % -------getting the velocity field using integration approach---------------------------------
       % first quadrant of mesh
        for i= (nn+1)/2 :(nn-1)      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :(nn-1)  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j+1)/a(2))^2 + (ZZ(i,j+1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j+1)-ZZ(i,j)); 
                Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j+1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + LL(2,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + LL(3,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
            end
            end
            ind= (YY(i+1,(nn+1)/2)/a(2))^2 + (ZZ(i+1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
        end
       
        % second quadrant of mesh
            for i= (nn+1)/2 :-1:2      % this corresponds to movement along Y-direction of mesh
                  for j= (nn+1)/2 :(nn-1)  % this corresponds to movement along X-direction of mesh
                     ind= (YY(i,j+1)/a(2))^2 + (ZZ(i,j+1)/a(3))^2;
                     if ind<1
                        Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j+1)-ZZ(i,j)); 
                        Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j+1)-ZZ(i,j));
                     end
                     if ind>1
                        Vel_field(1,1,i,j+1)= Vel_field(1,1,i,j) + LL(2,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                        Vel_field(2,1,i,j+1)= Vel_field(2,1,i,j) + LL(3,2,i,j+1)*(ZZ(i,j+1)-ZZ(i,j));
                     end
                  end
                  ind= (YY(i-1,(nn+1)/2)/a(2))^2 + (ZZ(i-1,(nn+1)/2)/a(3))^2;
                  if ind<1
                     Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                     Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                  end
                  if ind>1
                     Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                     Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                  end
            
            end
        
% third quadrant of mesh
        for i= (nn+1)/2 :-1:2      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :-1:2  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j-1)/a(2))^2 + (ZZ(i,j-1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j-1)-ZZ(i,j)); 
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j-1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + LL(2,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + LL(3,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
            end
            end
            ind= (YY(i-1,(nn+1)/2)/a(2))^2 + (ZZ(i-1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i-1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i-1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i-1,(nn+1)/2)*(YY(i-1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
       end
     %fourth quadrant of mesh
     
        for i= (nn+1)/2 :(nn-1)      % this corresponds to movement along Y-direction of mesh
            
            for j= (nn+1)/2 :-1:2  % this corresponds to movement along X-direction of mesh
            ind= (YY(i,j-1)/a(2))^2 + (ZZ(i,j-1)/a(3))^2;
            if ind<1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + ll(2,2)*(ZZ(i,j-1)-ZZ(i,j)); 
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + ll(3,2)*(ZZ(i,j-1)-ZZ(i,j));
            end
            if ind>1
                Vel_field(1,1,i,j-1)= Vel_field(1,1,i,j) + LL(2,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
                Vel_field(2,1,i,j-1)= Vel_field(2,1,i,j) + LL(3,2,i,j-1)*(ZZ(i,j-1)-ZZ(i,j));
            end
            end
            ind= (YY(i+1,(nn+1)/2)/a(2))^2 + (ZZ(i+1,(nn+1)/2)/a(3))^2;
            if ind<1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + ll(2,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + ll(3,3)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            if ind>1
                Vel_field(1,1,i+1,(nn+1)/2)= Vel_field(1,1,i,(nn+1)/2) + LL(2,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
                Vel_field(2,1,i+1,(nn+1)/2)= Vel_field(2,1,i,(nn+1)/2) + LL(3,3,i+1,(nn+1)/2)*(YY(i+1,(nn+1)/2)-YY(i,(nn+1)/2));
            end
            
        end
        
        %----------------------------------------------------------------------------------------------------------
        
        
        disp_field= Vel_field*tincr;
        
        % now use this displacement field to get new coordinate points
        % original points are YY and ZZ.
        
        % new position of points given by:
        [~,~,~,NN]= size(disp_field);
        new_pos= zeros(NN);
        
        X_disp= squeeze(disp_field(1,1,:,:));
        Y_disp= squeeze(disp_field(2,1,:,:));
        ZZ1= ZZ + X_disp; % these are the new positions of the points
        YY1= YY + Y_disp;
        
       ind  = (YY1./a(3)).^2 + (ZZ1./a(2)).^2 > 1;
       y_ex1  = Y(ind);
       z_ex1  = Z(ind);
        epp    = cat(1,y_ex1',z_ex1');
        k= zeros(1,numel(epp)/2);
        ep= [k;epp];
        YYY= ep(2,:);
        ZZZ= ep(3,:);
        plot(YYY,ZZZ,'*')
        num   = numel(YY); 
        qanew= QL(a,q,theta,dE,tincr); 
         %updating for next step  
         a=qanew(:,4);
         q=qanew(1:3,1:3);
end
