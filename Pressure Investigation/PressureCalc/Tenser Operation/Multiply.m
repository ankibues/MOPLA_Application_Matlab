function m = Multiply(X,y)
    % Double-index contraction of a 4th-order tensor and 2nd-order tensor.
    % Input: X is a 4th order tensor and y is a second order tensor.
    % Output: The Multiply(X,Y)is a second order tensor. 
    m = zeros(3,3);
    for i=1:3
        for j=1:3
            x = reshape(X(i,j,:,:),3,3);
            m(i,j) = sum(sum(x.*y));
        end
    end
end