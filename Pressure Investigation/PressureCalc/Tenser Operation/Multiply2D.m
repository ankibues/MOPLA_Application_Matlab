function m = Multiply2D(X,y)
    % Double-index contraction of a 4th-order tensor and 2nd-order tensor.
    % Input: X is a 4th order tensor and y is a second order tensor.
    % Output: The Multiply(X,Y)is a second order tensor. 
    m = zeros(2,2);
    for i=1:2
        for j=1:2
            x = reshape(X(i,j,:,:),2,2);
            m(i,j) = sum(sum(x.*y));
        end
    end
end