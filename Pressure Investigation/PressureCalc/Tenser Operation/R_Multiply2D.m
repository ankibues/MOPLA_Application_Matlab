function m = R_Multiply2D(x,Y)
    % Double-index contraction of a 2nd-order tensor and 4th-order tensor.
    % Input : x is a second order tensor and Y is a 4th order tensor.
    % Output: The R-Multiply(X,Y)is a second order tensor. 
    m = zeros(2,2);
    for i=1:2
        for j=1:2
            y = reshape(Y(:,:,i,j),2,2);
            m(i,j) = sum(sum(x.*y));
        end
    end
end