function Plot_ep_X_stress(a)

%Plotting pressure for exterior points along X-axis 

% shape and orientation of inclusion
%a=[10;2;1]; % a1, a2 and a3 are the semi-axes with a1>=a2>=a3     
ang = [0; 0; 0];  
% these spherical angles correspond to the case when a1 lie along Z-axis, a2 along X-axis and a3 along Y-axis 
q = Q(ang);

%exterior points
x=0:.5:8;
%y=0:.5:10;
Xpts=[x;zeros(1,length(x));zeros(1,length(x))];    % points along X-axis
%Ypts=[zeros(1,length(x));y;zeros(1,length(x))];    % points along Y-axis

% note, these points are defined in coordinate system of the flow, and
% should be converted to inclusion's coordinate system for calculation of
% exterior fields.
Xpts_inc= zeros(size(Xpts));
%Ypts_inc= zeros(size(Ypts));
for i=1:length(Xpts)
    Xpts_inc(:,i)= q*Xpts(:,i);
    %Ypts_inc(:,i)= q*Ypts(:,i);
end


 ind   = (Xpts_inc(1,:)./a(1)).^2 + (Xpts_inc(2,:)./a(2)).^2 + (Xpts_inc(3,:)./a(3)).^2 > 1;
 EP=Xpts_inc(:,ind);
 
 
 [~,N]= size(EP);
 
%  Matrix flow field
   L     = [1 0 0;0 -1 0;0 0 0]; 
   
   
%  Power law coefficients for matrix and clast
   Nm = 1;
% since linear isotropic case is being considered
   r= 10;
%  ------------------------------------------------------------------------   
%  decomposition of the bulk flow L into a strain rate tensor D and a vorticity 
%  tensor W        
   D = 0.5 * (L + L');
   
%  generate 4th-order identity tensors   
   [Jd, ~, ~, ~, b] = Jnb(); 
   
   C= 2*Jd ; % Viscous stiffness tensor of the matrix
   
   d= q*D*q';              % Transforming D and W with respect to Ellipsoidal axis system
   
   



 %----------------------------------------------------------------------
 % some exterior point[ Important note: these points are defined in
 % inclusion's coordinate system ]
 S_anal= zeros(1,N);
 S_num= zeros(1,N);
 
 parfor kkk=1:N
     ep= EP(:,kkk);
        
     % Solutions using quasi-analytical solutions after Jiang (2016)

     % Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
       [S_el,~] = SnP(a);

     % Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green tensor for pressure)
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
        dE       = Multiply(fdE(Nm,r,S,b,Jd),d);    %  strain rate tensor in inclusion                  
        u2= dE-d; 
     
        
 
    % Exterior Fields  
        %--------------------------------------------------------------------------
        % G
          G = Ex_Gtensor(a,ep);
        % LMABDA_Ex 
        LAMBDA_Ex   = zeros(3,3);
        for ii=1:3
            for j=1:3
                LAMBDA_Ex(ii,j) = -1/3.*(G(ii,j,1,1)+ G(ii,j,2,2)+ G(ii,j,3,3));
                LAMBDA_Ex(j,ii) = LAMBDA_Ex(ii,j);
            end
        end
        for iii=1:3
            t                = 1/3.*squeeze(G(1,1,iii,iii)+ G(2,2,iii,iii)+ G(3,3,iii,iii));
            LAMBDA_Ex(iii,iii) = squeeze(LAMBDA_Ex(iii,iii))+t;
        end
        % S_Ex & PI_Ex
        S_Ex   = zeros(3,3,3,3);
        PI_Ex  = zeros(3,3,3,3);
        delt   = eye(3);
        for ii=1:3
            for j=ii:3
                for k=1:3
                    for l=k:3
                        %S_Ex
                        S_Ex(ii,j,k,l)  = squeeze(G(ii,j,k,l))+ delt(k,l).*...
                                            squeeze(LAMBDA_Ex(ii,j));
                        S_Ex(j,ii,k,l)  = S_Ex(ii,j,k,l);
                        S_Ex(j,ii,l,k)  = S_Ex(ii,j,k,l);
                        S_Ex(ii,j,l,k)  = S_Ex(ii,j,k,l);
                        %PI_Ex
                        PI_Ex(ii,j,k,l) = 1/2.*(delt(j,k).*LAMBDA_Ex(ii,l) +...
                                            delt(j,l).*LAMBDA_Ex(ii,k) - delt(ii,k)...
                                            .*LAMBDA_Ex(j,l) - delt(ii,l).*...
                                            LAMBDA_Ex(j,k));
                        PI_Ex(ii,j,l,k) = PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,k,l) = -PI_Ex(ii,j,k,l);
                        PI_Ex(j,ii,l,k) = -PI_Ex(ii,j,k,l);
                    end
                end
            end
        end
        S_Ex(1,1,1,1) = -(S_Ex(1,1,2,2)+S_Ex(1,1,3,3));
        S_Ex(2,2,2,2) = -(S_Ex(2,2,1,1)+S_Ex(2,2,3,3));
        S_Ex(3,3,3,3) = -(S_Ex(3,3,1,1)+S_Ex(3,3,2,2));
  
       

      
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex, invS);
            v2          = Multiply(v1, u2);
            e_Ex = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            sigma_out1 = Multiply(C, e_Ex);     % deviatoric stress at an exterior point around the inclusion using isotropic solutions
            
            S_anal(1,kkk)= inva(sigma_out1);
  
   %--------------------USING ANISOTROPIC SOLUTIONS-----------------------
   % Gaussian points   
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
  
         
% stiffness tensor of the inclusion, 
   Ce    = 2*r*Jd; 
  
   


    Cm_a    = Transform(C, q);
    Carray  = C2OneDarray(Cm_a);
    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
    S = Contract(Jd, Contract(T,Cm_a));
    
    invS = Inverse((S),b);               
    H = Contract(Cm_a,(invS-Jd));           % Eq.12a Jiang 2014
   
    A = Contract(FourTensorInv(H + Ce),(H + Cm_a));                     
    e= Multiply(A,d);    % e is strain rate in ellipsoid
    u1= e-d;                  
    
            
 
%
        % Exterior Fields  
 %--------------------------------------------------------------------------
% Gaussian point calculation

 %--------------Calculation using code for the general Anisotropic medium.............
        
        
 %--------------------------------------------------------------------------
% Gaussian point calculation
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

    [Alpp3,Bett3,www3] = GaussGGQ(180);
    theta3 = reshape(Alpp3,1,[]);
    phi3   = reshape(Bett3,1,[]);
    Wout3 = reshape(www3,1,[]);
    [psi3,Win3] = Gauss1(25);

    [Alpp4,Bett4,www4] = GaussGGQ(180);
    theta4 = reshape(Alpp4,1,[]);
    phi4   = reshape(Bett4,1,[]);
    Wout4 = reshape(www4,1,[]);
    [psi4,Win4] = Gauss1(25);

    [Alpp5,Bett5,www5] = GaussGGQ(150);  
    theta5 = reshape(Alpp5,1,[]);
    phi5   = reshape(Bett5,1,[]);
    Wout5 = reshape(www5,1,[]);
    [psi5,Win5] = Gauss1(25);

    [Alpp6,Bett6,www6] = GaussGGQ(100);
    theta6 = reshape(Alpp6,1,[]);
    phi6   = reshape(Bett6,1,[]);
    Wout6 = reshape(www6,1,[]);
    [psi6,Win6] = Gauss1(25);

           % Calculating (T_Ex)
  
           
  
 
          tic
           
        % --------------------------------------------------------------
             lambda = solve_eq(a, ep);
             if lambda<= .41
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta1,phi1,Wout1,psi1,Win1);
             elseif lambda>.41 && lambda<= .84
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta2,phi2,Wout2,psi2,Win2);
             elseif lambda> .84 && lambda<= 1.76
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta3,phi3,Wout3,psi3,Win3);
             elseif lambda> 1.76 && lambda<= 2.25
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta4,phi4,Wout4,psi4,Win4); 
             elseif lambda> 2.25 && lambda<= 2.76
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta5,phi5,Wout5,psi5,Win5);    
             else 
                T_Ex = T_ExtAniso3D_vec(a,ep,Cm_a,theta6,phi6,Wout6,psi6,Win6);
             end

            % This approach allows us to choose number of Gaussian points 
            % and grid nodes based on exterior point's position from the 
            % inclusion-medium boundary-------------------------------------

% Exterior fields 
          
           
        % S_Ex: exterior Eshelby tensor     strain rate e_Ex
            S_Ex2 = Contract(Contract(Jd,T_Ex),Cm_a);
            hh=  Contract(S_Ex2, invS);
            e_Ex2= Multiply(hh,u1)+ d;
           
            sigma_out2= Multiply(Cm_a,e_Ex2); % deviatoric stress at an exterior point around the inclusion using anisotropic solutions
            S_num(1,kkk)= inva(sigma_out2);
           
          toc
      
   
 end
 
 XX= EP(2,:);
 
 scatter(XX, S_anal,'*')
 hold on
 scatter(XX, S_num,'s')
 
 
end