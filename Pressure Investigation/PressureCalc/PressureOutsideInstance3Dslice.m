% Pressure outside calculation of alpha values and their plots    

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 0;0 -1 0;0 0 0]; 
   
   
%  Power law coefficients for matrix and clast
   Nm=1;
     % since linear isotropic case is being considered
   r= 10;
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
   a= [2.5;1.5;1];
   
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
    sigmainside= 2*dE*r;
    sigmaininva= inva(sigmainside);
    
   

    Pdev_in= contract1(p_in,u2);
    Sigma_ellipsoid = 2* Multiply(Jd,dE)/ r;
    Sigma_Inva_ellipsoid = inva(Sigma_ellipsoid);
    

    % 3D meshgrid in clast's coordinate
        xgv     = -4:.1:4;           % grid vector: x'axis,a1
        ygv     = -3:.1:3;           % grid vector: y'axis,a2
        zgv     = -3:.1:3;   % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);

        % Exterior points
        %ind   = (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        % total points number
        num   = numel(X);
        % exterior points number
        [~,n] = size(ep);

        
        % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,ep);
        % LMABDA_Ex 
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
        w_Ex    = zeros(3,3,n);
        p_Ex    = zeros(1,n);
 
        for rr=1:n
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex(:,:,rr) = Multiply(2*Jd, e_Ex(:,:,rr));
    
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex(:,:,rr) = v4 + w;
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex(:,:,rr), 2*Jd);
            uu2          = R_Multiply(uu1, invS);
            uu3          = contract1(uu2, u2);
            p_Ex(rr)     = uu3 ; 
    
        end
        p1          = Pdev_in /Sigma_Inva;          % normalized with respect to far field stress invariant
        p2          = p_Ex ./ Sigma_Inva;
        %--------------------------------------------------
        
        pressure        = zeros(num,1);
        pressure(~ind)  = p1;
        pressure(ind)   = p2;
        pressure        = reshape(pressure,size(squeeze(Z)));
        
        figure
        isosurface(X,Y,Z,pressure,0)
        isosurface(X,Y,Z,pressure,.4)
        isosurface(X,Y,Z,pressure,-.4)
        %}
        
        