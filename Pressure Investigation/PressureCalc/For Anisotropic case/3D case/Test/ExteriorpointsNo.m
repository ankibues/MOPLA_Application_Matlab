% no. of points needed for exterior calculation  
a= [3;1.5;1];
% 3D meshgrid in clast's coordinate
        
        xgv     = 0:.1:4;           % grid vector: x'axis,a1
        ygv     = 0:.1:3;           % grid vector: y'axis,a2
        zgv     = 0:.1:2;   % grid vector: z'axis,a3
        [X,Y,Z] = meshgrid(xgv,ygv,zgv);
 
     %   [a1ang, a2ang, a3ang] = ConvertQ2Angs2(q);

        % Exterior points
        %ind   = (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        ind   = (X./a(1)).^2 + (Y./a(2)).^2 + (Z./a(3)).^2 > 1;
        x_ex  = X(ind);
        y_ex  = Y(ind);
        z_ex  = Z(ind);
        ep    = cat(1,x_ex',y_ex',z_ex');
        % total points number
        num   = numel(X);
        % exterior points number
        [~,n] = size(ep);
 lambda = solve_eq(a, ep);