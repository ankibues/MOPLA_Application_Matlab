% Exterior stress calculation in Elastic inclusion problem
  
  D= [0 .5 0;.5 0 0;0 0 0];

 
%  generate 4th-order identity tensors   
  [Jd, Js, Ja, Jm] = FourIdentity();
   
% elastic stiffnes: % Here, K -> bulk modulus taken as 1 G->shear modulus taken as 1, nu->
% poisson ratio     nu= (3K-2G/6K+2G)
%--------------------------------------------------- major parameters
%determining strength of the zone
r1=1;
r2= .00001;
 nu= .25; 
 Km= 1;
 % ------------------------------------------------
 mu_m= (1-2*nu)*3*Km/2/(nu+1);
% bulk modulus and shear modulus for inclusion, in terms of bulk modulus
% and shear modulus of medium
K= r1*Km;
mu= r2*mu_m;


   Cm = 3*Km*Jm + 2*mu_m*Jd;  
   Ce= 3*K*Jm + 2*mu*Jd;
   SIGMA= Multiply(Cm,D);
   SIGMA_Inva = inva(Multiply(Jd,SIGMA));

  
 
% initial state of the ellipsoid
   %a= sort(rand(3,1),'descend');                                      % initial length of the ellipsoid axis
   a= [1000;100;1];
   pp =0;
   ang = [0; 0; pp];                                                    % pp is the theta2 angle in paper           
   q = Q(ang);
   



    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
                                  
   % Eshelby Tensors S for Interior points using elliptical integrals  
    [S,~] = SnP_poisson(a,nu);
   
    invS= FourTensorInv_elastic(S);
   
    H= Contract(Cm, invS-Js); 
    A= Contract(FourTensorInv_elastic(H+Ce),H+Cm);
    
    dE=Multiply(A,d);  % dE is strain tensor in ellipsoid
    
    %{
    % case of a void
    dE= Multiply(FourTensorInv(Jd-S),d);
    %}
    u2= dE-d;
    sigmainside= Multiply(Ce,dE);
    sigmaininva= inva(Multiply(Jd,sigmainside));
      %}
   
    

    % 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = 95:.1:105;           % grid vector: y'axis,a-3
        zgv     = -15:.1:15;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);
         
        % Exterior points
        %ind   = (X./a(1)).^2 + ((Y-2)./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind   = (X./1).^2 + (((Y).*cos(pp)+ Z.*sin(pp))./a(2)).^2 + (((Y).*sin(pp)-Z.*cos(pp))/a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(ep);
        num   = numel(X);
        

        % total points number
        
        % exterior points number
        
        %-----------------------------------------------
        epp = zeros(3,n);
        qq= [1, 0 ,0; 0 ,cos(pp), sin(pp); 0 ,-sin(pp) ,cos(pp)]; % this part transforms(rotates) the coordinates for their use in G calculation. 
        for kk=1:n
            epp(:,kk)= qq*ep(:,kk);
        end
        %}
        %-----------------------------------------------

        % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Gtensor_El(a,epp,nu);
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
      
 
        for rr=1:n
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex(:,:,rr) = Multiply(Cm, e_Ex(:,:,rr));
            s_Exinva(:,rr)= inva(Multiply(Jd,s_Ex(:,:,rr)));
        
    
        end
     
        stress= zeros(num,1);
        Ss1= sigmaininva./SIGMA_Inva;
        Ss2= s_Exinva./SIGMA_Inva;
     
        stress(~ind)  = Ss1;
        stress(ind)   = Ss2;
        stress        = reshape(stress,size(squeeze(Z)));
        

            
               
        figure
        contourf(squeeze(Y),squeeze(Z),stress)
       
      
        
    
        
        
        
        