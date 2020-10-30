function m = R_Multiply(x,Y)
    % Double-index contraction of a 2nd-order tensor and 4th-order tensor.
    % Input : x is a second order tensor and Y is a 4th order tensor.
    % Output: The R-Multiply(X,Y)is a second order tensor. 
    m = zeros(3,3);
    for i=1:3
        for j=1:3
            y = reshape(Y(:,:,i,j),3,3);
            m(i,j) = sum(sum(x.*y));
        end
    end
end