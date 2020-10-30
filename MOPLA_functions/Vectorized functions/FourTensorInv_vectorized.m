function m = FourTensorInv_vectorized(X)
% FourTensorInv.m
% Vectorized form to calculate Fourth order tensor inverse for N tensors at once
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% still need to work on this !!
    [~,~,~,~,N]= size(X);
  sqrt2 = 2^0.5;
  sqrt3 = 3^0.5;
  
    b(:,:,1) = 1/(sqrt2*sqrt3) * [-1, 0, 0; 0, -1, 0; 0, 0, 2];
    b(:,:,2) = 1/sqrt2 * [-1, 0, 0; 0, 1, 0; 0, 0, 0];
    b(:,:,3) = 1/sqrt2 * [0, 0, 0; 0, 0, 1; 0, 1, 0];
    b(:,:,4) = 1/sqrt2 * [0, 0, 1; 0, 0, 0; 1, 0, 0];
    b(:,:,5) = 1/sqrt2 * [0, 1, 0; 1, 0, 0; 0, 0, 0];
    b(:,:,6) = 1/sqrt3 * [1, 0, 0; 0, 1, 0; 0, 0, 1]; 
    
    h = 5;
    M = zeros(h,h,N); %% This M is for N tensors simultaneously !
    for lambda = 1:h
        for xi = 1:h
            M(lambda,xi,:)= contract1_vectorized(Multiply_vectorized(X,b(:,:,xi)),b(:,:,lambda));  
        end
    end
    
    m = ConvertBK(multinv(M),h,b); %% work on this part!
end

function M = ConvertBK(x,h,b)
    M = zeros(3,3,3,3);
    bb = zeros(h,h);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                     for m=1:h
                         for n=1:h
                             bb(m,n)= b(i,j,m)*b(k,l,n); 
                         end
                     end
                   M(i,j,k,l)=sum(sum(x.*bb));
                end
            end
        end
    end
end
