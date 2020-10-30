function P = Pac(y1, y3, a1, a3, x, h, u2, Pin)
% Pressure on ac plane
% 29/03/2015

    m  = numel(y1);
    n  = numel(y3);
    mn = m*n;
    P  = zeros(mn,3);
    k  = 1;
    
    for i = 1:n
        for j =1:m
            P(k,2) = j;
            P(k,3) = i;
            rho    = y1(j)^2/(a1^2)+ y3(i)^2/(a3^2);
             if rho <= 1
                 P(k,1) = Pin;
             elseif rho > 1 && rho <= 1.5
                 P(k,1)  = ExP(x,[y1(j);0;y3(i)],h,u2,'apgq');
       
             else
                 P(k,1)  = ExP(x,[y1(j);0;y3(i)],h,u2,'glesh');
    
             end
             k = k+1;
        end
    end
end