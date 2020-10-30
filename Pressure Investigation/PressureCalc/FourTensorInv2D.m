function m = FourTensorInv2D(X)
% FourTensorInv2D.m
% Inverse of a 4th-order symmetric Tensor for 2D 
%--------------------------------------------------------------------------
    
sqrt2 = 2^0.5;
  
b(:,:,1) = 1/(sqrt2) * [1, 0; 0, -1];
b(:,:,2) = 1/sqrt2 * [0, 1; 1, 0];
b(:,:,3) = 1/sqrt2 * [1, 0; 0, 1];
    
    
    h = 2;
    M = zeros(h,h);
    for lambda = 1:h
        for xi = 1:h
            M (lambda,xi)= contract1(Multiply2D(X,b(:,:,xi)),b(:,:,lambda));  
        end
    end
    
    m = ConvertBK(eye(2)/M,h,b);
end

function M = ConvertBK(x,h,b)
    M = zeros(2,2,2,2);
    bb = zeros(h,h);
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                     for m=1:h
                         for n=1:h
                             bb(m,n)= b(i,j,m)*b(k,l,n); 
                         end
                     end
                   M(i,j,k,l)=sum(sum(x.*bb));
                end
            end
        end
    end
end
