function m = Contract_vectorized2(X,Y)
% Contract.m
% Double-index contraction between two 4th-order tensors for 'N' inclusions.

% Input:  X, a 4th order tensor 
% Y : 4th-order tensors, for 'N' number of inclusions;
% Output: m, a 4th-order tensor;
%--------------------------------------------------------------------------
    [~,~,~,~,N] = size(Y);
    m = zeros(3,3,3,3,N);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    x = repmat(reshape(X(i,j,:,:),3,3),1,1,N);
                    y = squeeze(Y(:,:,k,l,:));
                    m(i,j,k,l,:) = sum(sum(x.*y));
                end
            end
        end
    end
end
