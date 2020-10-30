% code to plot a randomly oriented ellipsoid

 ang=[0;0;3*pi/4];  % one can define this and plot ellipsoidds of different orientations
 q = Q(ang);


       [X1,Y1,Z1]= ellipsoid(0,0,0,a(1),a(2),a(3));
        XX1=X1(:);
        YY1=Y1(:);
        ZZ1=Z1(:);
        pts= [XX1';YY1';ZZ1'];    % these are the points on the surface of ellipsoid in ellipsoids coordinate system
        [~,NN]= size(pts);
        
        for rrr=1:NN                                % bringing pts to general coordinate system
            pts(:,rrr)= q'*pts(:,rrr);
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
        