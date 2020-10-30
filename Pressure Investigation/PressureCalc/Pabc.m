function [xcoord, ycoord, zcoord, P] = Pabc(y1, y2, y3, intcr, a1, a2, a3, x, h, u2, Pin)
% Pressure 3D
% 30/03/2015
% According to Healy2009

    l  = numel(-y1:intcr:y1);
    m  = numel(-y2:intcr:y2);
    n  = numel(-y3:intcr:y3);
    
    
    xcoord = zeros(l,m,n);
    ycoord = zeros(l,m,n);
    zcoord = zeros(l,m,n);
    
    P      = zeros(l,m,n);
    
    xi = 0;
    
    
    
    for x1 = -y1:intcr:y1
        xi = xi + 1;
        yi = 0;
      
        for x2 = -y2:intcr:y2
            yi = yi + 1;
            zi = 0;
            for x3 = -y3:intcr:y3
                zi = zi + 1;
                xcoord(xi,yi,zi) = x1;
                ycoord(xi,yi,zi) = x2;
                zcoord(xi,yi,zi) = x3;
            
                rho = x1^2/(a1^2) + x2^2/(a2^2) + x3^2/(a3^2);
             if rho <= 1
                 P(xi,yi,zi) = Pin;
             elseif rho > 1 && rho <= 1.5
                 P(xi,yi,zi)  = ExP(x,[x1;x2;x3],h,u2,'apgq');
       
             else
                 P(xi,yi,zi)  = ExP(x,[x1;x2;x3],h,u2,'glesh');
    
             end
             
            end
        end
    end
end