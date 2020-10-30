 function m = Contract(X,Y)
    % Double-index contraction between 2 4th-order tensors.
    % Input: X and Y are two 4th-order tensors.
    % Output:The Contract(X,Y)is a 4th-order tensor.
    m = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    x= reshape(X(i,j,:,:),3,3);
                    y= reshape(Y(:,:,k,l),3,3);
                    m(i,j,k,l)=sum(sum(x.*y));
                end
            end
        end
    end
end

