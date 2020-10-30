% Test of Exterior Pressure calculation Analytical vs Numerical 2D
% Anisotropic case !
 L     = [1 0 0;0 -1 0;0 0 0]; 
  D = 0.5 * (L + L');
   ang = [0; 0; 0];
    q = Q(ang);
    d= q*D*q'; 
a=[1000;5;1];
 [Jd, Js, Ja, Jm, b] = Jnb(); 
 Cm= 2*Jd;
 Nm=1;
 r=.1;
 %{
 Cm1    = zeros(3,3,3,3);
m= 1;
Nn=1;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end
 %}
%-------------------------------------------------------------------------
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
 
  invS= Inverse(S,b);
    h= Contract(2*Jd, invS);                        
    p_in = R_Multiply(p,h);
    dE       = Multiply(fdE(Nm,r,S,b,Jd),d);                  
                                                                     % dE is strain rate tensor in ellipsoid
    
    u2= dE-d;
 %----------------------------------------------------------------------
 % some exterior points
 
 ep= [0 ,0; 4,5; 4,4];
 % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,ep);
        % LMABDA_Ex 
        LAMBDA_Ex   = zeros(3,3,2);
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
         p_Ex    = zeros(1,2);
 
        for rr=1:2
       
        % p_Ex: exterior field pressure
            uu1          = R_Multiply(LAMBDA_Ex(:,:,rr), Cm);
            uu2          = R_Multiply(uu1, invS);
            uu3          = contract1(uu2, u2);
            p_Ex(rr)     = uu3 ; 
    
        end
        
        %--------------FOr anisotropic formulation.............
        

        
     L     = [1 0 ;0 -1]; 


%  two viscosities for anisotropic medium
   m = 1;                              % strength of anisotropy
   Nn =  1;
   tic
   [Alp,ww] = Gauss1(20000);  
   toc
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   

% initial state of the ellipse
                                     

 a = [5;1];                                             % initial length of the ellipse axis
 

 ang1 = (pi/180)* 0;                                          % angle of the ellipse axis with external coordinate
 
 %  Viscousity of isotropic inclusion(for r=100)
   etaE= r* Nn; 
 q = [cos(ang1), sin(ang1); -sin(ang1), cos(ang1)];
 
   
 % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  
Cm(1,1,1,1)= Nn*((cos(2*ang1))^2 + (sin(2*ang1))^2/m); 
Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang1)*cos(2*ang1);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang1))^2 + (cos(2*ang1))^2/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)= Cm(1,1,1,2);
Cm(2,1,1,1)= Cm(1,2,1,1);
Cm(1,1,2,1)= Cm(1,1,1,2);
Cm(2,1,2,2)= Cm(1,2,2,2);
Cm(2,2,2,1)= Cm(2,2,1,2);
Cm(2,1,1,2)= Cm(1,2,1,2);
Cm(2,1,2,1)= Cm(2,1,1,2);
Cm(1,2,2,1)= Cm(2,1,2,1);
  
% stiffness tensor of the ellipsoid, 
  Ce    = 2*etaE*Jd; 


% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Ja,ang1,Nn,m,Cm,Alp,ww);

% Interior fields (e,w,p)

% describe D,W,C in the clast's coordinate system 
  D_bar  = q * D * q';  % here, pure shear is with respect to anisotropy !(i.e. the general coordinate system)
  W_bar  = q * W * q';   % in  case, when we consider pure shear w.r.t. inclusion, D_bar= D;
   
  Ce    = Transform2D(Ce, q);  
  S_bar  = Multiply2D(Cm, D_bar);
  SI     = inva2D(S_bar);                        % stress invariant for far-field for normalization 
  P      = SI;                               

% e: interior field strain-rate
  invS= FourTensorInv2D(S);
  H       = Contract2D(Cm,(invS - Jd));   %Hill's constraint tensor
  t1      = FourTensorInv2D(H + Ce);
  A       = Contract2D(t1,(H + Cm));      %strain-rate partitioning tensor
  e       = Multiply2D(A, D_bar);
  
% s: interior field deviatoric-stress
  s       = Multiply2D(Ce, e);
% w: interior field vorticity
  de      = e - D_bar;                  
  t2      = Multiply2D(Contract2D(PI, invS), de);
  w       = t2 + W_bar;
% p: interior field pressure
  b1      = R_Multiply2D(LAMBDA, Cm);
  b2      = contract1(R_Multiply2D(b1, invS), de);
  p       = b2 ; 

  

  
 Pdev_in= p;
 

       % Exterior Fields  
        %--------------------------------------------------------------------------
        
        % S_Ex & PI_Ex
        S_Ex = zeros(2,2,2,2,2);
        PI_Ex= zeros(2,2,2,2,2);
        LAMBDA_Ex = zeros(2,2,2);
        
        
       
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex    = zeros(2,2,2);
        s_Ex    = zeros(2,2,2);
        w_Ex    = zeros(2,2,2);
        p_Ex2    = zeros(1,2);
 
        epp= [ep(2,:); ep(3,:)];
        for rr=1:2
         % Exterior Eshelby tensor calculation
         
            [S_Ex(:,:,:,:,rr),LAMBDA_Ex(:,:,rr),PI_Ex(:,:,:,:,rr)] =  SnpExtAniso(a,Jd,Ja,ang1,Nn,m,epp(:,rr),Cm,Alp,ww);
            
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract2D(S_Ex(:,:,:,:,rr), invS);
            v2          = Multiply2D(v1, de);
            e_Ex(:,:,rr) = v2 + D_bar;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex(:,:,rr) = Multiply2D(Cm, e_Ex(:,:,rr));
    
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract2D(PI_Ex(:,:,:,:,rr), invS);
            v4          = Multiply2D(v3, de);                    % u2 difference:  dE-d
            w_Ex(:,:,rr) = v4 + W_bar;
        % p_Ex: exterior field pressure
            uu1          = R_Multiply2D(LAMBDA_Ex(:,:,rr), Cm);
            uu2          = R_Multiply2D(uu1, invS);
            uu3          = contract1(uu2, de);
            p_Ex2(rr)     = uu3 ;      
        end

  

        