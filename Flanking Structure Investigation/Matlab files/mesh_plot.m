
 a= [100;5;.2];   
 pp= pi/4;
   ang = [0;0;pp] ;                                          % initial angles of the ellipsoid axis
   q = Q(ang);
% 3D meshgrid in external coordinate
        xgv     = -5:.01:5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = -5:.01:5;           % grid vector: y'axis,a2
        zgv     = 0; %-5:.01:5;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
        x_ex  = X(Y==2);
        y_ex  = Y(Y==2);
        z_ex  = Z(Y==2);
        ep    = cat(1,x_ex',y_ex',z_ex');
        
        
        x_ex_1  = X(Y==0);
        y_ex_1  = Y(Y==0);
        z_ex_1  = Z(Y==0);
        ep_1    = cat(1,x_ex_1',y_ex_1',z_ex_1');
        
        x_ex_2  = X(Y==-2);
        y_ex_2  = Y(Y==-2);
        z_ex_2  = Z(Y==-2);
        ep_2    = cat(1,x_ex_2',y_ex_2',z_ex_2');
        
        % 3D meshgrid in internal coordinate
        xgv1     = 0; %-5:.05:5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv1    = -5:.01:5;           % grid vector: y'axis,a2
        zgv1     = -5:.01:5;  % grid vector: z'axis,a3
        [X1,Y1,Z1] = meshgrid(xgv1,ygv1,zgv1);
        ind1= (X1./a(1)).^2 + (Y1./a(2)).^2 + (Z1./a(3)).^2;
        x_ex1  = X1(ind1<1);
        y_ex1  = Y1(ind1<1);
        z_ex1  = Z1(ind1<1);
        epp    = cat(1,x_ex1',y_ex1',z_ex1'); % These are points corresponding the inclusion
        [~,NN]= size(epp);
        
        for rrr=1:NN                                % bringing to general coordinate system
           epp(:,rrr)= q'*epp(:,rrr);
        end
        
           
        
       
        X_marker1= ep(1,:);
        Y_marker1= ep(2,:);
        
        X_marker2= ep_1(1,:);
        Y_marker2= ep_1(2,:);
        
        X_marker3= ep_2(1,:);
        Y_marker3= ep_2(2,:);
        
        XX= epp(1,:);
        YY= epp(2,:);
        
        plot(X_marker1,Y_marker1,'.')
        hold on
        plot(X_marker2,Y_marker2,'.')
        hold on
        plot(X_marker3,Y_marker3,'.')
        hold on
        plot(XX,YY,'.')
  %}

% input should be a x,y coordinates 
% what the code does: check the point in the mesh closest to the
% coordinate, and use that point's velocity field for the input coordinates
% output should be the velocity field at the given coordinates.


%{
%-----------------------------------------------------------------------------------------------------------------
% first lets try with simple case of just single coordinates.

x= 1.232;
y= 2.342;

load('test_data.mat','YY','ZZ','Vel_field'); % loads mesh and velocity field for the mesh.
dist= ((x-ZZ).^2 + (y-YY).^2).^(0.5);
minimum = min(min(dist));
[x,y]=find(dist==minimum);
Vel= Vel_field(:,:,x,y);

%----------------------------------------------------------------------------------------------------------------
%}

%{
x= 1.232;
y= 2.342;
n=251*251;
ep= [x;y];
ep= repmat(ep,1,n);

load('test_data.mat','YY','ZZ','Vel_field'); % loads mesh and velocity field for the mesh.
[~,N]= size(YY);

vel_f = zeros(2,1,n);
tic
for i=1:n
    dist= ((ep(1,i)-ZZ).^2 + (ep(2,i)-YY).^2).^(0.5);
    minimum = min(min(dist));
    [xx,yy]=find(dist==minimum);
    vel_f(:,:,i)= Vel_field(:,:,xx,yy);
end
toc
%}
