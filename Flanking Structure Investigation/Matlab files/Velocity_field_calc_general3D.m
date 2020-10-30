function V = Velocity_field_calc_general3D(X,Y,Z,a,invS,Jd,d,u2,w,ll)
% This function calculates velocity field at any coordinate point outside the inclusion
% Input: X,Y,Z : coordinates of the point; a: 3*1 matrix representing 

% This function will take the input as the coordinates of the marker point,and gives back the velocity field at that point.
% Once we get the coordinate, we decide how it should be divided into the
% steps.

% Instead of using the cartesian coordinates, Here, I use Ellipsoidal
% coordinates for efficient integration of velocity gradient field.
% Note: The actual calculation is done in Cartesian coordinates;
% Ellipsoidal coordinates are only used to get the shortest path of integrating
% velocity gradient field(along the height axis of the ellipsoidal coordinate system)
% This method is more accurate and efficient.



[lat,long,h]= CartesiantoEllipsoidal(X,Y,Z,a);     % bringing the points to Ellipsoidal coordinates

% choosing the known points for velocity field
% these known points in elliptic coordinates
%[~,~,ho]= CartesiantoEllipsoidal(a(1),0,0,a);    this is not needed since
%ho , that lies on the Ellipsoidal surface, has h=0.
ho=0;
% these known points in cartesian coordinates, so that velocity field can
% be calculated.
[Xo,Yo,Zo]= EllipsoidaltoCartesian(lat,long,ho,a);

Vo= ll*[Xo;Yo;Zo];                % reference velocity field at the boundary


%{
Later On...Using the Lambda approach, we can give variable lengths of
 intervals. For closer points, more intervals, while for farther points,
less intervals at farther distances.
%}
lambda = solve_eq(a,[X;Y;Z]);
%
% Lambda_External 
         if lambda<= .6920
             kk= .0001;
         elseif lambda>.64 && lambda<= 1.80
             kk= .001;
         else 
             kk= .01;            
         end
%}
%kk=.0001;
%kk=.001;
H= ho:kk:h;
LAT= lat*ones(size(H)); % these are interval coordinates in ellipsoidal coordinate system
LONG= long*ones(size(H)); 

% bringing them to cartesian system
ex= (a(1)^2- a(3)^2)/a(1)^2;
ee= (a(1)^2- a(2)^2)/a(1)^2;
N= a(1)./sqrt(1-(ex*(sin(LAT)).^2)-(ee*(cos(LAT)).^2).*(sin(LONG)).^2);

x= (N+H).*cos(LAT).*cos(LONG);

y=(N.*(1-ee)+H).*cos(LAT).*sin(LONG);

z=(N*(1-ex)+H).*sin(LAT);

% choosing the intervals for integration; this will be along h axis only.

ep= [x;y;z];

[~,N1]=size(ep);
  
    
    
    
    
  % Exterior Fields  along this path
   %--------------------------------------------------------------------------
     % G
          G1 = Ex_Gtensor(a,ep);
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
        for rr=1:N1
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
% Removing one term from the series of Ls calculated from all points.
    L_Ex1(:,:,N1)= [];
   
  ep1= [zeros(3,1), ep(:,1:N1-1)]; 
  % This process is done to calculate the delta x and delta y values for
  % each intervals. Note: the intervals are equal in Elliptic form, but not
  % equal in Cartesian form. So, this has to be done.
  delta= ep-ep1;
  delta(:,1)=[];
  deltaX= delta(1,:);
  deltaY= delta(2,:);
  deltaZ= delta(3,:);
  
 
  V1= Vo(1) + sum(deltaX.*(squeeze(L_Ex1(1,1,:)))',2) + sum(deltaY.*(squeeze(L_Ex1(1,2,:)))',2) + sum(deltaZ.*(squeeze(L_Ex1(1,3,:)))',2);
  V2= Vo(2) + sum(deltaX.*(squeeze(L_Ex1(2,1,:)))',2) + sum(deltaY.*(squeeze(L_Ex1(2,2,:)))',2) + sum(deltaZ.*(squeeze(L_Ex1(2,3,:)))',2);
  V3= Vo(3) + sum(deltaX.*(squeeze(L_Ex1(3,1,:)))',2) + sum(deltaY.*(squeeze(L_Ex1(3,2,:)))',2) + sum(deltaZ.*(squeeze(L_Ex1(3,3,:)))',2) ;
  
  
    
    V= [V1;V2;V3];




end
