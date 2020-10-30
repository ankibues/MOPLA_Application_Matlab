function m=Inverse(X,b) 
    % The inverse operation of a 4th-order tensor. Choose 5 or 6 for
    % incompressible or compressible materials.
    m=ConvertBK(Convert(X,5,b)^-1,5,b);
end

function M = Convert(X,h,b)
    M = zeros(h,h);
    for lambda = 1:h
        for xi = 1:h
            M (lambda,xi)= contract1(Multiply(X,reshape(b(xi,:,:),3,3)),reshape(b(lambda,:,:),3,3));  
        end
    end
end

function M = ConvertBK(x,h,b)
    M = zeros(3,3,3,3);
    bb = zeros(h,h);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                     for m=1:h
                         for n=1:h
                             bb(m,n)= b(m,i,j)*b(n,k,l); 
                         end
                     end
                   M(i,j,k,l)=sum(sum(x.*bb));
                end
            end
        end
    end
end
