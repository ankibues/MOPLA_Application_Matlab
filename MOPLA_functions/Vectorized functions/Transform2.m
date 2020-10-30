function ss=Transform2(X,Q,N)
% Transform2.m    
% 4th-order tensor transformation between coordinate systems.( Vectorized form) 
%
% Input:  X, 3*3*3*3 matrix;
%         Q, 3*3*N matrix;
%         N, number of elements 
% Output: ss, 3*3*3*3*N matrix.
%--------------------------------------------------------------------------
    ss = zeros(3,3,3,3,N);
    qq = zeros(3,3,3,3,N);
    
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                    for m=1:3
                        for n=1:3
                            for s=1:3
                                for t=1:3
                                    qq(m,n,s,t,:)= Q(i,m,:).*Q(j,n,:).*Q(k,s,:).*Q(l,t,:);
                                end
                            end
                        end
                    end
                    kkk= reshape(qq,[],1) .* repmat(reshape(X,[],1),N,1);
                
                    ss(i,j,k,l,:)= sum(reshape(kkk,[81,N]))';
                end
            end
        end
    end
end
