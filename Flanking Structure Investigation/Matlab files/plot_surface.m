% plot surface
% 3D meshgrid in external coordinate
        xgv     = -1.5:.25:1.5;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = zeros(size(xgv));    % grid vector: y'axis,a2 ----since the points are in X-Z plane
        zgv     = -1.5:.25:1.5;  %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
        
        [Xx,Yy,Zz]= meshgrid(xgv,ygv,zgv);
        [~,~,nnn]= size(Xx);
        k1= Xx(:);
        k2= Yy(:);
        k3= Zz(:);
        
        %ep    = [k1';k2';k3'];   % basically, if I have the points, I can generate the plot.

        [~,N]= size(ep);

sf= fit([ep(1,:)', ep(2,:)'], ep(3,:)','linearinterp');
%set(sf, 'edgecolor','none')
figure
plot(sf)
shading interp
camlight
hold on
%
% code to plot a randomly oriented ellipsoid
% a=[1;.5;.1];
% ang=[0;0;3*pi/4];  % one can define this and plot ellipsoidds of different orientations
% q = Q(ang);


       [X1,Y1,Z1]= ellipsoid(0,0,0,a(1),a(2),a(3));
        XX1=X1(:);
        YY1=Y1(:);
        ZZ1=Z1(:);
        pts= [XX1';YY1';ZZ1'];    % these are the points on the surface of ellipsoid in ellipsoids coordinate system
        [~,NN]= size(pts);
        R= [1 0 0; 0 0 1;0 -1 0];
        for rrr=1:NN                                % bringing pts to general coordinate system
            pts(:,rrr)= q'*pts(:,rrr);
            pts(:,rrr)= R*pts(:,rrr);
        end 
        
        XX1= pts(1,:);
        YY1= pts(2,:);
        ZZ1= pts(3,:);
        
        % getting points in the form of meshgrid
        X_mesh= reshape(XX1,[21,21]);
        Y_mesh= reshape(YY1,[21,21]);
        Z_mesh= reshape(ZZ1,[21,21]);
        % this can be plotted using the surf command
        surf(X_mesh,Y_mesh,Z_mesh)
        shading interp
        %}