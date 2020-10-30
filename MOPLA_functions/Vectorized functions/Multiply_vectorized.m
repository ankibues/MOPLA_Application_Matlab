function m = Multiply_vectorized(X,y)
% Multiply_vectorized.m
% Double-index contraction of 4th-order tensor to a 2nd-order tensor for N cases all
% at once.
% again, this works when all 4th order tensors are to be multiplied by a
% single 2nd order tensor
% Input: X, a 4th-order tensor, 3*3*3*3*N matrix;
%        y, a 2nd-order tensor, 3*3 matrix;

% Output: m, a 2nd-order tensor, 3*3*N matrix. 
%--------------------------------------------------------------------------
[~,~,~,~,N] = size(X);    
m = zeros(3,3,N);
    for i=1:3
        for j =1:3
             x = reshape(X(i,j,:,:,:),3,3,N);
             m(i,j,:) = sum(sum(x.*repmat(y,1,1,N)));
        end
    end
end