function m = contract1(x,y)
    % Double-index contraction between 2 2nd-order tensors. (this can be
    % used for both 2D and 3D)
    % (This is the contract(x,y) function in mathcad!)
    % Input:x and y are second order tensors.
    % Output is a scalar.
    m = sum(sum(x.*y));
end