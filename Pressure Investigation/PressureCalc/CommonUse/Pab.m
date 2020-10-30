function P = Pab(y1, y2, a1, a2, x, h, u2, Pin, Alp, Bet, ww)
% Pressure on ab plane
% 29/03/2015

    m  = numel(y1);
    n  = numel(y2);
    mn = m*n;
    P  = zeros(mn,3);
    k  = 1;
    
    for i = 1:n
        for j =1:m
            P(k,2) = j;
            P(k,3) = i;
            rho    = y1(j)^2/(a1^2)+ y2(i)^2/(a2^2);
             if rho <= 1
                 P(k,1) = Pin;
             elseif rho > 1 && rho <= 1.01
                 P(k,1)  = ExP(x,[y1(j);y2(i);0],h,u2,'apgq',Alp,Bet,ww);
       
             else
                 P(k,1)  = ExP(x,[y1(j);y2(i);0],h,u2,'ggq',Alp,Bet,ww);
    
             end
             k = k+1;
        end
    end
end