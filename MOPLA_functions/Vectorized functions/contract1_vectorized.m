function m = contract1_vectorized(x,y)
% contract1.m
% Double-index contraction between two 2nd-order tensors.
% x is 3*3*N matrix      y is 3*3 matrix
%--------------------------------------------------------------------------
    [~,~,N]= size(x);
    m = sum(sum(x.*repmat(y,1,1,N)));
end