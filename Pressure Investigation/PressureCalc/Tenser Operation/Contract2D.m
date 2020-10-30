 function m = Contract2D(X,Y)
    % Double-index contraction between 2 4th-order tensors.
    % Input: X and Y are two 4th-order tensors.
    % Output:The Contract(X,Y)is a 4th-order tensor.
    m = zeros(2,2,2,2);
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    x= reshape(X(i,j,:,:),2,2);
                    y= reshape(Y(:,:,k,l),2,2);
                    m(i,j,k,l)=sum(sum(x.*y));
                end
            end
        end
    end
end

