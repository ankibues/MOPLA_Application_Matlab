% test of coordinate system conversion

a=[2;.8;.1];
%
 % 3D meshgrid in clast's coordinate
        xgv     = -2:.1:2 ;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     = -2:.1:2;           % grid vector: y'axis,a2
        zgv     = -2:.1:2;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);
         
        % Exterior points
        %ind   = (X./a(1)).^2 + ((Y-2)./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(ep);
        
        EP= zeros(size(ep));
        
        for i= 1:n
            [lat,long,h]= CartesiantoEllipsoidal(ep(1,i),ep(2,i),ep(3,i),a);
            [X,Y,Z]= EllipsoidaltoCartesian(lat,long,h,a);
            EP(:,i)= [X;Y;Z];
        end
        
        K1= round(EP,3);
        K2= round(ep,3);
        isequal(K1,K2);
        
 %}
 %{
 ep=[5;0;0];
 tic
[lat,long,h]= CartesiantoEllipsoidal(ep(1),ep(2),ep(3),a);
 toc
 
 tic
[X,Y,Z]= EllipsoidaltoCartesian(lat,long,h,a);
toc
 %}