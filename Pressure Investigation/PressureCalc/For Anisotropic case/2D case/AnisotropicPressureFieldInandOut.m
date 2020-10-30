% Internal and External Pressure field for an isotropic ellipse in Linear Anisotropic Material

%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [0 2 ;0 0]; 

%  time increment of each step during the computation; see Jiang (2007) for the choice of tincr.   
   tincr = 0.05;                         % step length for computation
%  total steps of the computation
   steps = 10;
%  number of computation steps between output sets 
   mm=10; 
%  Viscousity of isotropic inclusion
   r= 100;
%  two viscosities for anisotropic medium
   m = 25;                              % strength of anisotropy
   Nm =  1;
   Nn = m*Nm;
            
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   [Jd, Js, Ja, Jm] = FourIdentity2D();
   


% initial state of the ellipse
                                     

   a= [5;2];                                             % initial length of the ellipse axis
   R=5;
   ang = 0 ;                                          % initial angle of the ellipse axis
 
   q = [cos(ang), sin(ang); -sin(ang), cos(ang)];
   epsilonII=0.5;   
   
 Pdev_steps_alpha1= zeros(steps,1);
 Pdev_steps_alpha2= zeros(steps,1);
 Q_evl= zeros(2,2,steps);
 A_evl= zeros(2,steps);
 strainininvacheck= zeros(1,steps);
 sigmaininvacheck= zeros(1,steps);
for i=1:steps  
   
   % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  Cm    = zeros(2,2,2,2);
  
Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + (sin(2*ang))^2/m); Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + (cos(2*ang))^2/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)=Cm(1,1,1,2);
Cm(2,1,1,1)= Cm(1,2,1,1);
Cm(1,1,2,1)= Cm(1,1,1,2);
Cm(2,1,2,2)= Cm(1,2,2,2);
Cm(2,2,2,1)= Cm(2,2,1,2);
Cm(2,1,1,2)= Cm(1,2,1,2);
Cm(2,1,2,1)= Cm(2,1,1,2);
Cm(1,2,2,1)= Cm(2,1,2,1);
  
% stiffness tensor of the ellipsoid, 
  Ce    = 2*r*Jd; % lets confirm this later, what should be done here


% Interior Fields 
%--------------------------------------------------------------------------
% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
   [S,LAMBDA,PI] = SnpInAniso(a,Jd,Ja,ang,Nn,m,Cm);

% Interior fields (e,w,p)

% describe D,W,C in the clast's coordinate system 
  D_bar  = q * D * q';
  W_bar  = q * W * q';
   
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
  p       = b2 ; %+ P;

  
 wE = Wd2D(a,W_bar,e);
 theta= w-wE;
 qanew= QL2D(a,q,theta,e,tincr);   
 %  write updated Q to Q_evl        
 Q_evl(:,:,i)=qanew(1:2,1:2);
 %  write updated a to A_evl
 A_evl(:,i)= qanew(:,3);
  
 Pdev_in= p;
 Sigma_ellipse = 2* Multiply2D(Jd,e)* r;
 Sigma_Inva_ellipse = inva2D(Sigma_ellipse);
 Pdev_steps_alpha1(i,1)= Pdev_in/SI;            % normalized with respect to far field stress invariant
 Pdev_steps_alpha2(i,1) = Pdev_in/Sigma_Inva_ellipse ; 
 
 if i==10 % || i==20 || i==30 || i==40 || i==50 || i==60 || i==70 || i==80 || i==90 || i==100 || i==110 || i==120
        % 3D meshgrid in clast's coordinate
        xgv     = 0;           % grid vector: x'axis,a1
        ygv     = -5:0.5:5;           % grid vector: y'axis,a2
        zgv     = -5:0.5:5;    % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);

        % Exterior points
        ind   = (X./1).^2 + (Y./a(1)).^2 + (Z./a(2)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,y_ex',z_ex');
        % total points number
        num   = numel(X);
        % exterior points number
        [~,n] = size(ep);


        % Exterior Fields  
        %--------------------------------------------------------------------------
        
        % S_Ex & PI_Ex
        [S_Ex,LAMBDA_Ex,PI_Ex] = SnpExtAniso(a,Jd,Ja,ang,Nn,m,ep,Cm);
        
       
        % Exterior fields (e_Ex,s_Ex,w_Ex,p_Ex)
        e_Ex    = zeros(2,2,n);
        s_Ex    = zeros(2,2,n);
        w_Ex    = zeros(2,2,n);
        p_Ex    = zeros(1,n);
 
        for rr=1:n
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
        
        %subplot(4,3,i/10)
        contourf(squeeze(Y),squeeze(Z),pressure)
        %hold on
        axis equal
        
        
 end
end

 a=qanew(:,3);
 q=qanew(1:2,1:2); 
 ang= atan(q(1,2)/q(1,1));


   
   
