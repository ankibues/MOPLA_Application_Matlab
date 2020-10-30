% Internal and External Pressure field for an isotropic ellipse in Linear Anisotropic Material
% Getting 2-D cross sections using 3D formulations
%  Input parameters--------------------------------------------------------
%  Matrix flow field
   L     = [1 0 0 ;0 -1 0; 0 0 0]; 

   gp                = 20;
    [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
  
      
   %  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   W = 0.5 * (L - L');  
%  generate 4th-order identity tensors   
   
   [Jd, Js, Ja, Jm, b] = Jnb(); 

% initial state of the ellipse
                                     

 a= [3;1.5;1];                                             % initial length of the ellipse axis
 
 r = .100;
 

Cm    = zeros(3,3,3,3);
m= 10;
Nn=1;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end
  %}
   % load('AnisotropicC_200inclusionsB','C_bar_evl');
 % stiffness tensor of the matrix(in ellipse coordinate system),(Jiang 2016, Eq-52)
  %Cm =  C_bar_evl(:,:,:,:,200);      % stiffness tensor for the matrix
 
  eta= 1;%3.3; % Average HEM viscosity
% stiffness tensor of the ellipsoid, 
  Ce    = 2*eta*r*Jd; 
  
   Sigma = Multiply(Cm,D);        % far field stress value( sigma= 2*eta*Jd:E)i.e Viscosity(eta)of matrix taken as 1.  
   Sigma_Inva = inva(Sigma);


% Interior Fields 
%--------------------------------------------------------------------------
     i= 0;
     ang = [0; 0;(pi/180)*i];                        % See Jiang 2007a for choice of angles.
    q = Q(ang);
    d= q*D*q';                                    % Transforming D and W to Ellipsoidal axis system
   
    Cm_a    = Transform(Cm, q);
    Carray  = C2OneDarray(Cm_a);
    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
    S = Contract(Jd, Contract(T,Cm_a));
    
    lambda = SnpInAniso3DPress(a,Cm_a);
    invS = Inverse((S),b);
    h= Contract(Cm_a, invS);                        
    p_in = R_Multiply(lambda,h);
    H = Contract(Cm_a,(invS-Js));           % Eq.12a Jiang 2014
   
    A = Contract(FourTensorInv(H + Ce),(H + Cm_a));                     
    e= Multiply(A,d);    % e is strain rate in ellipsoid
    u1= e-d;                  
    sigmainside= Multiply(Ce,e);
    sigmaininva= inva(sigmainside);
       
    Pdev_in= contract1(p_in,u1);
    
 

 % 3D meshgrid in clast's coordinate ( For Y-Z cross-section, X has 0 values)
        xgv     = 0;           % grid vector: x'axis,a1
        ygv     = 0:.1:3;           % grid vector: y'axis,a2
        zgv     = 0:.1:2;   % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
      %{   
   % 3D meshgrid in clast's coordinate ( For X-Y cross-section, Z has 0 values)
        xgv     = 0;           % grid vector: x'axis,a1
        ygv     = 0:.1:3;           % grid vector: y'axis,a2
        zgv     = 0:.1:2;   % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);      
    %}    
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);

        % Exterior points
        %ind   = (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        % total points number
        num   = numel(Y);
        % exterior points number
        [~,n] = size(ep);
      
 lambda = solve_eq(a, ep);  % function for calculating lamda value for each point

 %
        % Exterior Fields  
 %--------------------------------------------------------------------------
% Gaussian point calculation

[Alpp1,Bett1,www1] = GaussGGQ(280);
theta1 = reshape(Alpp1,1,[]);
phi1   = reshape(Bett1,1,[]);
Wout1 = reshape(www1,1,[]);
[psi1,Win1] = Gauss1(25);

[Alpp2,Bett2,www2] = GaussGGQ(200);
theta2 = reshape(Alpp2,1,[]);
phi2   = reshape(Bett2,1,[]);
Wout2 = reshape(www2,1,[]);
[psi2,Win2] = Gauss1(25);

[Alpp3,Bett3,www3] = GaussGGQ(100);
theta3 = reshape(Alpp3,1,[]);
phi3   = reshape(Bett3,1,[]);
Wout3 = reshape(www3,1,[]);
[psi3,Win3] = Gauss1(25);

[Alpp4,Bett4,www4] = GaussGGQ(70);
theta4 = reshape(Alpp4,1,[]);
phi4   = reshape(Bett4,1,[]);
Wout4 = reshape(www4,1,[]);
[psi4,Win4] = Gauss1(25);

[Alpp5,Bett5,www5] = GaussGGQ(60);
theta5 = reshape(Alpp5,1,[]);
phi5   = reshape(Bett5,1,[]);
Wout5 = reshape(www5,1,[]);
[psi5,Win5] = Gauss1(25);

[Alpp6,Bett6,www6] = GaussGGQ(50);
theta6 = reshape(Alpp6,1,[]);
phi6   = reshape(Bett6,1,[]);
Wout6 = reshape(www6,1,[]);
[psi6,Win6] = Gauss1(25);
%-----------------------------------------------------------------

%{

  
        % Exterior fields 
     
        T_Ex= zeros(3,3,3,3,n);
        sigma_out= zeros(1,n);
        e_Ex= zeros(3,3,n);
        parfor rr=1:n

        % Lambda_External 
         if lambda(1,rr)<= .41
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta1,phi1,Wout1,psi1,Win1);
         elseif lambda(1,rr)>.41 && lambda(1,rr)<= .84
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta2,phi2,Wout2,psi2,Win2);
         elseif lambda(1,rr)> .84 && lambda(1,rr)<= 1.76
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta3,phi3,Wout3,psi3,Win3);
         elseif lambda(1,rr)> 1.76 && lambda(1,rr)<= 2.25
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta4,phi4,Wout4,psi4,Win4);
         elseif lambda(1,rr)> 2.25 && lambda(1,rr)<= 2.76
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta5,phi5,Wout5,psi5,Win5);    
         else 
             T_Ex(:,:,:,:,rr) = T_ExtAniso3D_vec(a,ep(:,rr),Cm_a,theta6,phi6,Wout6,psi6,Win6);
         end
         
        
        % S_Ex: exterior Eshelby tensor     strain rate e_Ex
            S_Ex = Contract(Contract(Jd,T_Ex(:,:,:,:,rr)),Cm_a);
            hh=  Contract(S_Ex, invS);
            e_Ex(:,:,rr)= Multiply(hh,u1)+ d;
           
            sigma= Multiply(Cm_a,e_Ex(:,:,rr));
            sigma_out(:,rr)= inva(sigma);
        end
        s1          = sigmaininva /Sigma_Inva;          % normalized with respect to far field stress invariant
        s2          = sigma_out./  Sigma_Inva;
      
        
        stress        = zeros(num,1);
        stress(~ind)  = s1;
        stress(ind)   = s2;
        stress        = reshape(stress,size(squeeze(Y)));
        
        figure
        
        contourf(squeeze(Y),squeeze(Z),stress)
        save('cross_section_press_high_res.mat','stress');
%}