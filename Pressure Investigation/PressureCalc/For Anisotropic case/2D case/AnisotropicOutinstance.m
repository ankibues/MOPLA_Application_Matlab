% Internal and External Pressure field for an isotropic ellipse in Linear Anisotropic Material

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 ;0 -1]; 


%  two viscosities for anisotropic medium
   m = 10;                              % strength of anisotropy
   Nn =  1;
   tic
   [Alp,ww] = Gauss1(5000);  
   toc
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   

% initial state of the ellipse
                                     

 a= [5;1];                                             % initial length of the ellipse axis
 
 r = 0;
 ang = (pi/180)* 30;                                          % angle of the ellipse axis with external coordinate
 
 %  Viscousity of isotropic inclusion(for r=100)
   etaE= r* Nn; 
 q = [cos(ang), sin(ang); -sin(ang), cos(ang)];
 
   
 % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  
Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + (sin(2*ang))^2/m); 
Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + (cos(2*ang))^2/m);
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
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Ja,ang,Nn,m,Cm,Alp,ww);

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
 
% Pdev_steps_alpha1(i,1)= Pdev_in/SI;            % normalized with respect to far field stress invariant
% Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipse ; 
 
 %if i==10 % || i==20 || i==30 || i==40 || i==50 || i==60 || i==70 || i==80 || i==90 || i==100 || i==110 || i==120
        % 2D meshgrid in clast's coordinate
              
        ygv     = -8:.1:8;           % grid vector: y'axis,a2
        zgv     = -8:.1:8;    % grid vector: z'axis,a3
        [Y,Z] = meshgrid(ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);

        % Exterior points
       % ind   =  (Y./a(1)).^2 + (Z./a(2)).^2 > 1;
       ind   = ((Y.*cos(ang)+ Z.*sin(ang))./a(1)).^2 + ((Y.*sin(ang)-Z.*cos(ang))/a(2)).^2 > 1;
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,y_ex',z_ex');
        % total points number
        num   = numel(Y);
        % exterior points number
        [~,n] = size(ep);
 %-----------------------------------------------
        
        epp = zeros(2,n);
       % this part transforms(rotates) the coordinates for their use in calculation. 
        for kk=1:n
            epp(:,kk)= q*ep(:,kk);
        end

        % Exterior Fields  
        %--------------------------------------------------------------------------
        
        % S_Ex & PI_Ex
        S_Ex = zeros(2,2,2,2,n);
        PI_Ex= zeros(2,2,2,2,n);
        LAMBDA_Ex = zeros(2,2,n);
        
        
       
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex    = zeros(2,2,n);
        s_Ex    = zeros(2,2,n);
        w_Ex    = zeros(2,2,n);
        p_Ex    = zeros(1,n);
 
        parfor rr=1:n
         % Exterior Eshelby tensor calculation
         
            [S_Ex(:,:,:,:,rr),LAMBDA_Ex(:,:,rr),PI_Ex(:,:,:,:,rr)] =  SnpExtAniso(a,Jd,Ja,ang,Nn,m,epp(:,rr),Cm,Alp,ww);
            
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
            p_Ex(rr)     = uu3 ; 
    
        end
        p1          = Pdev_in /SI;
        p2          = p_Ex ./ SI;
        %--------------------------------------------------
        pressure        = zeros(num,1);
        pressure(~ind)  = p1;
        pressure(ind)   = p2;
        pressure        = reshape(pressure,size(squeeze(Z)));
        figure
        %subplot(4,3,i/10)
        contourf(squeeze(Y),squeeze(Z),pressure)
        %hold on
        axis equal
        
        save('pressuredata',pressure);
        
   



   
   
