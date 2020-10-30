%function V = Velocity_field_calc2(X,Y,a,invS,Jd,d,u2,w,ll)
% this file is just to test the velocity field calculation at individual
% points.

X= 1.5;
Y= 1.5;
tic
% This function will take the input as the coordinates of the marker point,and gives back the velocity field at that point.
% Once we get the coordinate, we decide how it should be divided into the
% steps.
%lambda = solve_eq(a,[0;X;Y]);
%{
% Lambda_External 
         if lambda<= .64
             kk= .0001;
         elseif lambda>.64 && lambda<= 1.36
             kk= .001;
         elseif lambda> 1.36 && lambda<= 2.16
             kk= .01;
         else 
             kk= .01;
         end
%}
kk=.0001;

% (X>0 && X<a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y<0 && Y>-a(3)) || (X>0 && X<a(2) && Y<0 && Y>-a(3))
 %    kk=.0001;
%end
% For choosing the known point for the velocity field Xo and Yo.
if (Y>=0) && (Y>=a(3))
    Xo= 0;
    Yo= a(3);
elseif (Y<0) && (Y<=-a(3))
    Xo= 0;
    Yo= -a(3);
elseif X>=0
    Xo= a(2);
    Yo= 0;
elseif X<0
    Xo= -a(2);
    Yo= 0;
end

% choosing the intervals for integration

if (X<0 && X<Xo) || (X>0 && X<Xo)
    x= X:kk:Xo;
elseif (X<0 && X>Xo) || (X>0 && X>Xo) 
    x= Xo:kk:X;
else
    x=0;
end

if Y<0
    y= Y:kk:Yo;
elseif Y>0 
    y= Yo:kk:Y;
else
    y=0;
end

[~,N2]=size(y); 
[~,N1]=size(x);
                                        
  kkk1= zeros(1,N1);
  kk1= Yo*ones(1,N1);
  KK1= Y*ones(1,N1);
  kk2= Xo*ones(1,N2);
  kkk2= zeros(1,N2);
  
  % these are the points where external velocity gradient tensor needs to be calculated
  if (X>0 && X<a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y<0 && Y>-a(3)) || (X>0 && X<a(2) && Y<0 && Y>-a(3))
      ep1    = cat(1,kkk1,x,KK1);        % moving along X axis  
      ep2    =  cat(1,kkk2,kk2,y);       % moving along Y axis
      
  else
      ep1    = cat(1,kkk1,x,kk1);        % moving along X axis  -----These are applicable for almost all cases except some regions.
      ep2    =  cat(1,kkk2,kk2,y);       % moving along Y axis
  end
  
  
    Vo= ll*[0;Xo;Yo];                % reference velocity field at the boundary
    
    
    
  % Exterior Fields  along X-axis
   %--------------------------------------------------------------------------
     % G
          G1 = Ex_Gtensor(a,ep1);
        % LAMBDA_Ex 
        LAMBDA_Ex1   = zeros(3,3,N1);
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
        S_Ex1   = zeros(3,3,3,3,N1);
        PI_Ex1  = zeros(3,3,3,3,N1);
        
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
        e_Ex1    = zeros(3,3,N1);
        s_Ex1    = zeros(3,3,N1);
        s_Exinva1= zeros(1,N1);
        w_Ex1    = zeros(3,3,N1);
        L_Ex1    = zeros(3,3,N1);
        parfor rr=1:N1
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex1(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex1(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex1(:,:,rr) = Multiply(2*Jd, e_Ex1(:,:,rr));
            s_Exinva1(:,rr)= inva(s_Ex1(:,:,rr));
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex1(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex1(:,:,rr) = v4 + w;
        % L_Ex: exterior field of velocity gradient tensor
            L_Ex1(:,:,rr)     = e_Ex1(:,:,rr)+ w_Ex1(:,:,rr)  ;
        end
        
%------------------------------------------------------------------------------------------------------------------------

 % Exterior Fields  along Y-axis
   %--------------------------------------------------------------------------
     % G
          G2 = Ex_Gtensor(a,ep2);
        % LAMBDA_Ex 
        LAMBDA_Ex2   = zeros(3,3,N2);
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
        S_Ex2   = zeros(3,3,3,3,N2);
        PI_Ex2  = zeros(3,3,3,3,N2);
        
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
        e_Ex2    = zeros(3,3,N2);
        s_Ex2    = zeros(3,3,N2);
        s_Exinva2= zeros(1,N2);
        w_Ex2    = zeros(3,3,N2);
        L_Ex2    = zeros(3,3,N2);
        parfor rr=1:N2
        % e_Ex: exterior field strain-rate in clast's coordinate
            v1          = Contract(S_Ex2(:,:,:,:,rr), invS);
            v2          = Multiply(v1, u2);
            e_Ex2(:,:,rr) = v2 + d;
        % s_Ex: exterior field stress in clast's coordinate
            s_Ex2(:,:,rr) = Multiply(2*Jd, e_Ex2(:,:,rr));
            s_Exinva2(:,rr)= inva(s_Ex2(:,:,rr));
        % w_Ex: exterior field vorticity in clast's coordinate
            v3          = Contract(PI_Ex2(:,:,:,:,rr), invS);
            v4          = Multiply(v3, u2);                    % u2 difference:  dE-d
            w_Ex2(:,:,rr) = v4 + w;
        % L_Ex: exterior field of velocity gradient tensor
            L_Ex2(:,:,rr)     = e_Ex2(:,:,rr)+ w_Ex2(:,:,rr)  ;
        end
        
%------------------------------------------------------------------------------------------------------------------------
% Removing one term from the series of Ls calculated from all points.
    if (X<0 && X>Xo) || (X>0 && X>Xo)
        L_Ex1(:,:,N1)= [];
    else
        L_Ex1(:,:,1)= [];
    end
    
    if (Y>0)
        L_Ex2(:,:,N2)= [];
    else
        L_Ex2(:,:,1)= [];  
    end
    
% assigning negative signs for fields where the movement is towards a
% negative direction.

  if Y<Yo
      L_Ex2= -L_Ex2;
  end
  
  if (X<0 && X<Xo) || (X>0 && X<Xo)
      L_Ex1= -L_Ex1;
  end
  
  
  %---------------Performing the integration.
  if (X>0 && X<a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y>0 && Y<a(3)) || (X<0 && X>-a(2) && Y<0 && Y>-a(3)) || (X>0 && X<a(2) && Y<0 && Y>-a(3))
      %Integration along Y-axis first in this case
    V1 = Vo(2,1) + (kk)* sum(L_Ex2(2,3,:));
    V2 = Vo(3,1) + (kk)* sum(L_Ex2(3,3,:)); 
    %Integration along X-axis
    V11 = V1 + (kk)*sum(L_Ex1(2,2,:));
    V22 = V2 + (kk)*sum(L_Ex1(3,2,:));
  else   
%Integration along X-axis
    V1 = Vo(2,1) + (kk)*sum(L_Ex1(2,2,:));
    V2 = Vo(3,1)  + (kk)*sum(L_Ex1(3,2,:));

%Integration along Y-axis
    V11 = V1 + (kk)* sum(L_Ex2(2,3,:));
    V22 = V2 + (kk)* sum(L_Ex2(3,3,:)); 
  end
  
    
    V= [V11;V22];
toc

%}

%end
