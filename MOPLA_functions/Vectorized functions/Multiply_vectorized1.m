function m = Multiply_vectorized1(X,y)
% Multiply_vectorized.m
% Double-index contraction of 4th-order tensor to a 2nd-order tensor for N cases all
% at once.
% *********, this is for the case when all 4th order tensors are to be
% multiplied by N different 2nd order tensors
% Input: X, a 4th-order tensor, 3*3*3*3*N matrix;
%        y, a 2nd-order tensor, 3*3*N matrix;

% Output: m, a 2nd-order tensor, 3*3*N matrix. 
%--------------------------------------------------------------------------
[~,~,~,~,N] = size(X);    
m = zeros(3,3,N);
    for i=1:3
        for j =1:3
             x = reshape(X(i,j,:,:,:),3,3,N);
             m(i,j,:) = sum(sum(x.*y));
        end
    end
end