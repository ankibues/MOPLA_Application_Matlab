
% Motion of single ellipsoid in Power Law Isotropic Material, with finite
% strain calculated at an external point, that changes its position with
% every time step



%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 0;0 -1 0;0 0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.01;                         % step length for computation
%  total steps of the computation
   steps = 100;
   
%  Power law coefficients for matrix and clast
   Nm=1;
   Nc=1; 
   r=100;
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
   %a= sort(rand(3,1),'descend');                                      % initial length of the ellipsoid axis
   a= [10;9.9;1];
 
   ang = [0;0;0] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
   epsilonII=0.5;                                                      % initial strain rate invariant
% calculations to follow:
  Q_evl= zeros(3,3,steps);
  A_evl= zeros(3,steps);
  f=eye(3);
  Xo= 12;         % The exterior point to be followed
for i=1:steps
    d= q*D*q';                                      % Transforming D and W with respect to Ellipsoidal axis system
    w= q*W*q'; 
    Cm_a    = Transform(Cm, q);
    S_bar  = Multiply(Cm_a, d);
    SI     = inva(S_bar);
    P      = SI;
      
% ****************
%Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
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
%****************    
    
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
    sigmainside= 2*dE*r_new;
    strainininva = inva(dE);
    sigmaininva= inva(sigmainside);
    qanew= QL(a,q,theta,dE,tincr);       
    
    %  write updated Q to Q_evl        
        Q_evl(:,:,i)=qanew(1:3,1:3);
    %  write updated a to A_evl
        A_evl(:,i)= qanew(:,4);
        
        
        
    %*****************************************************
    
          
        % Exterior points (Instead of making a grid, this time I consider a
        % point independently !)
       
        EP= zeros(3,1);
        EP(2,:)= Xo;
        
        % exterior points number
        [~,n] = size(EP);
        
        % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,EP);
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
        e_Exinva= zeros(1,n);
        w_Ex    = zeros(3,3,n);
        p_Ex    = zeros(1,n);
 
        for rr=1:n
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex(:,:,rr) = v2 + d;
            e_Exinva(:,rr)= inva(e_Ex(:,:,rr));
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex(:,:,rr) = v4 + w;
            %{
        don't need these here !
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex(:,:,rr) = Multiply(2*Jd, e_Ex(:,:,rr));
            s_Exinva(:,rr)= inva(s_Ex(:,:,rr));
        
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex(:,:,rr), 2*Jd/Nm);
            uu2          = R_Multiply(uu1, invS);
            uu3          = contract1(uu2, u2);
            p_Ex(rr)     = uu3 ; 
         %}
        end
        l_Ex= e_Ex+ w_Ex;
        
        f =(eye(3) + l_Ex*tincr)*f; 
        
        [Ain,~]=finite_strain_ellipsoid(f);
        % Flinn diagram evolution
        x = log(Ain(2)./Ain(3));
        y = log(Ain(1)./Ain(2));
        
        plot(x,y,'.r',0:0.5:1.5,0:0.5:1.5,'-k')
        hold on
        xlabel('ln(a2/a3)')
        ylabel('ln(a1/a2)')
        axis square
        hold on
        
        
        
        
        Xo= Xo + l_Ex(1,1)*Xo*tincr;
        a=qanew(:,4);
        q=qanew(1:3,1:3);
end 

 

