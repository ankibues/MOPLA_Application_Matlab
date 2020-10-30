a= [100;10;1];
 pp=0;
% 3D meshgrid in clast's coordinate
        xgv     = 0;%-4:.05:4;          %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
        ygv     =  -8:.5:8;           % grid vector: y'axis,a2
        zgv     = -8:.5:8;  % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);
         
        % Exterior points
       %ind   = (X./a(1)).^2 + ((Y-2)./a(2)).^2 + (Z./a(3)).^2 > 1;
        indd   = (X./1).^2 + (((Y).*cos(pp)+ Z.*sin(pp))./a(2)).^2 + (((Y).*sin(pp)-Z.*cos(pp))/a(3)).^2 < 1;
        x_ex  = X(indd);
        y_ex  = Y(indd);
        z_ex  = Z(indd);
        ep    = cat(1,x_ex',y_ex',z_ex');
        [~,n] = size(ep);
        num   = numel(X);
        XX= ep(2,:);
        YY= ep(3,:);
        plot (XX,YY,'*')
        xlim([-8 8])
        ylim([-8 8])