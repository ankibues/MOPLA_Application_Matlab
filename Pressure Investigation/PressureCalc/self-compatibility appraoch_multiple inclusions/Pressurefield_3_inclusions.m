% This script calculates pressure field in and around three heterogeneities close to each other
% using the self compatibility approach by Novak et al (2012)


%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 0;0 -1 0;0 0 0]; 
   
   
%  Power law coefficients for matrix 
   Nm=1;   % linearly isotropic matrix
   
% viscosity ratios of two clasts with respect to the matrix   
   r1= 10;
   
   r2= 10;
   r3=10;
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
  
   a1= [10;2;1];                                                        % initial length of the ellipsoid 1 axis
   a2= [10;2;1]; % initial length of the ellipsoid 2 axis
   a3= [10;2;1];
   
   pp1 = 0;  % these values can be applied for different orientations of both clasts
   pp2 = 0;
   pp3 = 0;
   ang1 = [0; 0; pp1];                                                    % pp1 is the theta angle in Jiang and Bhandari (2018)
   ang2 = [0; 0; pp2];                                                    % pp2 is the theta angle in Jiang and Bhandari (2018)
   ang3 = [0; 0; pp3];
   
   q1 = Q(ang1);
   q2 = Q(ang2);
   q3= Q(ang3);
   
    d1= q1*D*q1';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w1= q1*W*q1';                                      % w is the vorticity of material, in ellipsoid 1's frame of reference 
    
    d2= q2*D*q2';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w2= q2*W*q2';                                      % w is the vorticity of material, in ellipsoid 1's frame of reference 
    
    d3= q3*D*q3';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w3= q3*W*q3';                                      % w is the vorticity of material, in ellipsoid 1's frame of reference 
    
    
    
% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el1,PI_el1] = SnP(a1);
  [S_el2,PI_el2] = SnP(a2);
  [S_el3,PI_el3] = SnP(a3);
  
% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green
% tensor for pressure)   -----------------For both inclusions
% p1_______For 1st inclusion
  p1   = zeros(3,3);
  for j=1:3
      p1(j,j) = -1/3* (S_el1(j,j,1,1)+ S_el1(j,j,2,2)+ S_el1(j,j,3,3)); 
  end
% S  
  S1       = S_el1;
  for k=1:3
      for l=1:3
          S1(k,k,l,l) = p1(k,k)+ S_el1(k,k,l,l);
      end
  end
 
 % p2_______For 2nd inclusion
  p2   = zeros(3,3);
  for j=1:3
      p2(j,j) = -1/3* (S_el2(j,j,1,1)+ S_el2(j,j,2,2)+ S_el2(j,j,3,3)); 
  end
% S  
  S2       = S_el2;
  for k=1:3
      for l=1:3
          S2(k,k,l,l) = p2(k,k)+ S_el2(k,k,l,l);
      end
  end
  
  % p3_______For 3rd inclusion
  p3   = zeros(3,3);
  for j=1:3
      p3(j,j) = -1/3* (S_el3(j,j,1,1)+ S_el3(j,j,2,2)+ S_el3(j,j,3,3)); 
  end
% S  
  S3       = S_el3;
  for k=1:3
      for l=1:3
          S3(k,k,l,l) = p3(k,k)+ S_el3(k,k,l,l);
      end
  end
  
  
  
  
  
  
  % -----Calculation of internal fields-----
  
  % for 1st inclusion 
    invS1= Inverse(S1,b);
    h1= Contract(2*Jd/Nm, invS1);                        
    p_in1 = R_Multiply(p1,h1);
    e1       = Multiply(fdE(Nm,r1,S1,b,Jd),d1);   % dE1 is strain rate tensor in ellipsoid 1               
    u1= e1-d1;                          
    Pdev_in1= contract1(p_in1,u1);                  % Pdev_in1 _pressure field in ellipsoid 1
    
  % for 2nd inclusion 
    invS2= Inverse(S2,b);
    h2= Contract(2*Jd/Nm, invS2);                        
    p_in2 = R_Multiply(p2,h2);
    e2       = Multiply(fdE(Nm,r2,S2,b,Jd),d2);   % dE2 is strain rate tensor in ellipsoid 2               
    u2= e2-d2;                          
    Pdev_in2= contract1(p_in2,u2);                  % Pdev_in2 _pressure field in ellipsoid 2
  
   % for 3rd inclusion 
    invS3= Inverse(S3,b);
    h3= Contract(2*Jd/Nm, invS3);                        
    p_in3 = R_Multiply(p3,h3);
    e3       = Multiply(fdE(Nm,r3,S3,b,Jd),d3);   % dE2 is strain rate tensor in ellipsoid 2               
    u3= e3-d3;                          
    Pdev_in3= contract1(p_in3,u3);                  % Pdev_in2 _pressure field in ellipsoid 2
%}    
  %----------Self-compatibility algorithm
  
   K=6;
   EP1= [0;K;0];
   EP2= [0;-K;0];
   
   
   ll=0;
   dE1_o= 0;
   dE2_o= 0;
   dE3_o= 0;
   
delta_all= zeros(1,1000);
dE_check= zeros(3,3,1000);
  for i=1:1000
   
   % each iterative calculation: 
   %internal field
   dE_1 = Multiply(fdE(Nm,r1,S1,b,Jd),d1);
   dE_2 = Multiply(fdE(Nm,r2,S2,b,Jd),d2);
   dE_3= Multiply(fdE(Nm,r3,S3,b,Jd),d3);
   
   u1= dE_1-d1;
   u2= dE_2-d2;
   u3= dE_3-d3;
   
   %external field
   Ex1_2= Extfield(a2,EP1,invS2,u2,d2);  % at 1 because of 2
   Ex1_2= q1*(q2'*Ex1_2*q2)*q1';  % brought to their corresponding coordinate system
   
   Ex1_3= Extfield(a3,2*EP2,invS3,u3,d3);
   Ex1_3= q1*(q3'*Ex1_3*q3)*q1';  % brought to their corresponding coordinate system
   
   
   Ex2_1= Extfield(a1,EP2,invS1,u1,d1);  
   Ex2_1= q2*(q1'*Ex2_1*q1)*q2';  % brought to their corresponding coordinate system
   
   Ex2_3= Extfield(a3,EP2,invS3,u3,d3);  
   Ex2_3= q2*(q3'*Ex2_3*q3)*q2';  % brought to their corresponding coordinate system
   
   
   Ex3_1= Extfield(a1,2*EP2,invS1,u1,d1);
   Ex3_1= q3*(q1'*Ex3_1*q1)*q3';
   
   Ex3_2= Extfield(a2,EP2,invS2,u2,d2);
   Ex3_2= q3*(q2'*Ex3_2*q2)*q3';
   
   
   d1= Ex1_2 + Ex1_3;    
   d2= Ex2_1 + Ex2_3;
   d3= Ex3_1 + Ex3_2;      
   
   dE_check(:,:,i) = dE_1;
   
   %delta_E = abs(inva(dE_1)+ inva(dE_2)) - abs(inva(dE1_o)+ inva(dE2_o));
   %kkk = delta_E;
   
   delta_E1 = abs(inva(dE_1)-inva(dE1_o)); 
   delta_E2 = abs(inva(dE_2)-inva(dE2_o));
   delta_E3 = abs(inva(dE_3)-inva(dE3_o));
   kkk =  max(delta_E1,delta_E2);
   
   delta_all(1,i) = kkk; 
   
   %if ll~=0 && kkk <.01
    % break
   %end 
   ll=ll+1;
   dE1_o= dE_1;
   dE2_o= dE_2;
  end
 xx= 1:1000;
 plot(xx,delta_all);
  
    %}

  d1= q1*D*q1';
  d2= q2*D*q2';  
   u1= e1-d1;
   u2= e2-d2;
 
  %-------------------------------------------------------------------------------
%{    
  % -----Calculation of external fields-----  

    % 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =  -8:.1:8;           % grid vector: y'axis,a2
        zgv     = -8:.1:8;  % grid vector: z'axis,a3
       [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
       
    % Exterior points
    
    % for 1st inclusion
        %ind   = (X./a(1)).^2 + ((Y-2)./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind1   = (X./1).^2 + (((Y-K).*cos(pp1)+ Z.*sin(pp1))./a1(2)).^2 + (((Y-K).*sin(pp1)-Z.*cos(pp1))/a1(3)).^2 > 1;
        x_ex1  = X(ind1);
        y_ex1  = Y(ind1);
        z_ex1  = Z(ind1);
        ep1    = cat(1,x_ex1',y_ex1',z_ex1');
        ep1(2,:) = ep1(2,:) - K;   % this is done to shift the ellipsoid from origin to a different point along y-axis
        [~,n1] = size(ep1);
        num   = numel(X);  
        % total points number
        
        % exterior points number
        
        %-----------------------------------------------
        epp1 = zeros(3,n1);
        qq1= [1, 0 ,0; 0 ,cos(pp1), sin(pp1); 0 ,-sin(pp1) ,cos(pp1)]; % this part transforms(rotates) the coordinates for their use in G calculation. 
        for kk=1:n1
            epp1(:,kk)= qq1*ep1(:,kk);
        end
        
        %-----------------------------------------------

   % for 2nd inclusion
        ind2   = (X./1).^2 + (((Y+K).*cos(pp2)+ Z.*sin(pp2))./a2(2)).^2 + (((Y+K).*sin(pp2)-Z.*cos(pp2))/a2(3)).^2 > 1;
        x_ex2  = X(ind2);
        y_ex2  = Y(ind2);
        z_ex2  = Z(ind2);
        ep2    = cat(1,x_ex2',y_ex2',z_ex2');
        ep2(2,:) = ep2(2,:) + K;   % this is done to shift the ellipsoid from origin to a different point along y-axis
        [~,n2] = size(ep2);
        
        % total points number
        % exterior points number
        
        %-----------------------------------------------
        epp2 = zeros(3,n2);
        qq2= [1, 0 ,0; 0 ,cos(pp2), sin(pp2); 0 ,-sin(pp2) ,cos(pp2)]; % this part transforms(rotates) the coordinates for their use in G calculation. 
        for kk=1:n2
            epp2(:,kk)= qq2*ep2(:,kk);
        end
      
        %-----------------------------------------------
      
          
  % Exterior Fields for 1st inclusion
        %--------------------------------------------------------------------------
        % G
          G1 = Ex_Gtensor(a1,epp1);
        % LMABDA_Ex 
        LAMBDA_Ex1   = zeros(3,3,n1);
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
        S_Ex1   = zeros(3,3,3,3,n1);
        PI_Ex1  = zeros(3,3,3,3,n1);
        
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
        e_Ex1    = zeros(3,3,n1);
        s_Ex1    = zeros(3,3,n1);
        s_Exinva1= zeros(1,n1);
        p_Ex1    = zeros(1,n1);
 
        for rr=1:n1
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex1(:,:,:,:,rr), invS1);
            v2          = Multiply(v1, u1);
            e_Ex1(:,:,rr) = v2 + d1;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex1(:,:,rr) = Multiply(2*Jd, e_Ex1(:,:,rr));
            s_Exinva1(:,rr)= inva(s_Ex1(:,:,rr));
        
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex1(:,:,rr), 2*Jd/Nm);
            uu2          = R_Multiply(uu1, invS1);
            uu3          = contract1(uu2, u1);
            p_Ex1(rr)     = uu3 ; 
    
        end
        pp1          = Pdev_in1 /Sigma_Inva;          % normalized with respect to far field stress invariant
       pp2          = p_Ex1 ./ Sigma_Inva;           
        pressure1        = zeros(num,1);
        %stress= zeros(num,1);
        
        pressure1(~ind1)  = pp1;
        pressure1(ind1)   = pp2;
        pressure1        = reshape(pressure1,size(squeeze(Z)));
      
        
        
 % Exterior Fields for 2nd inclusion
        %--------------------------------------------------------------------------
        % G
          G2 = Ex_Gtensor(a2,epp2);
        % LMABDA_Ex 
        LAMBDA_Ex2   = zeros(3,3,n2);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex2(ii,j,:) = -1/3.*(G2(ii,j,1,1,:)+ G2(ii,j,2,2,:)+ G2(ii,j,3,3,:));
                LAMBDA_Ex2(j,ii,:) = LAMBDA_Ex2(ii,j,:);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G2(1,1,iii,iii,:)+ G2(2,2,iii,iii,:)+ G2(3,3,iii,iii,:));
            LAMBDA_Ex2(iii,iii,:) = squeeze(LAMBDA_Ex2(iii,iii,:))+t;
        end
        % S_Ex & PI_Ex
        S_Ex2   = zeros(3,3,3,3,n2);
        PI_Ex2  = zeros(3,3,3,3,n2);
        
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex2(ii,j,k,l,:)  = squeeze(G2(ii,j,k,l,:))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex2(ii,j,:));
                        S_Ex2(j,ii,k,l,:)  = S_Ex2(ii,j,k,l,:);
                        S_Ex2(j,ii,l,k,:)  = S_Ex2(ii,j,k,l,:);
                        S_Ex2(ii,j,l,k,:)  = S_Ex2(ii,j,k,l,:);
                        %PI_Ex
                        PI_Ex2(ii,j,k,l,:) = 1/2.*(delt(j,k).*LAMBDA_Ex2(ii,l,:) +...
                                            delt(j,l).*LAMBDA_Ex2(ii,k,:) - delt(ii,k)...
                                            .*LAMBDA_Ex2(j,l,:) - delt(ii,l).*...
                                            LAMBDA_Ex2(j,k,:));
                        PI_Ex2(ii,j,l,k,:) = PI_Ex2(ii,j,k,l,:);
                        PI_Ex2(j,ii,k,l,:) = -PI_Ex2(ii,j,k,l,:);
                        PI_Ex2(j,ii,l,k,:) = -PI_Ex2(ii,j,k,l,:);
                    end
                end
            end
        end
        S_Ex2(1,1,1,1,:) = -(S_Ex2(1,1,2,2,:)+S_Ex2(1,1,3,3,:));
        S_Ex2(2,2,2,2,:) = -(S_Ex2(2,2,1,1,:)+S_Ex2(2,2,3,3,:));
        S_Ex2(3,3,3,3,:) = -(S_Ex2(3,3,1,1,:)+S_Ex2(3,3,2,2,:));
  
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex2    = zeros(3,3,n2);
        s_Ex2    = zeros(3,3,n2);
        s_Exinva2= zeros(1,n2);
        p_Ex2    = zeros(1,n2);
 
        for rr=1:n2
        % e_Ex: exterior field strain-rate in clast's coordinate
            vv1          = Contract(S_Ex2(:,:,:,:,rr), invS2);
            v2          = Multiply(vv1, u2);
            e_Ex2(:,:,rr) = v2 + d2;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex2(:,:,rr) = Multiply(2*Jd, e_Ex2(:,:,rr));
            s_Exinva2(:,rr)= inva(s_Ex2(:,:,rr));
        
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex2(:,:,rr), 2*Jd/Nm);
            uu2          = R_Multiply(uu1, invS2);
            uu3          = contract1(uu2, u2);
            p_Ex2(rr)     = uu3 ; 
    
        end
        ppp1          = Pdev_in2 /Sigma_Inva;          % normalized with respect to far field stress invariant
       ppp2          = p_Ex2 ./ Sigma_Inva;           
      
        pressure2        = zeros(num,1);
       
        pressure2(~ind2)  = ppp1;
        pressure2(ind2)   = ppp2;
        pressure2        = reshape(pressure2,size(squeeze(Z)));       
        
        
  %-------------------------------------------------------------------------------      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
       pressure= pressure1+pressure2; 
        figure
        
        contourf(squeeze(Y),squeeze(Z),pressure)
       
        axis equal
      
        figure
        
        contourf(squeeze(Y),squeeze(Z),pressure1,'ShowText','on')
       
        axis equal
        
        figure
        
        contourf(squeeze(Y),squeeze(Z),pressure2,'ShowText','on')
       
        axis equal
    
        
        
        %}
        
        
        
        