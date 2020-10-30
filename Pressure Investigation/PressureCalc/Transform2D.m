function ss=Transform2D(X,Q)
% Transform.m    
% 4th-order tensor transformation between coordinate systems. 
%
% Input:  X, 2*2*2*2 matrix;
%         Q, 2*2 matrix;
%         N, number of elements 
% Output: ss, 2*2*2*2 matrix.
%--------------------------------------------------------------------------
    ss = zeros(2,2,2,2);
    qq = zeros(2,2,2,2);
    
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    for m=1:2
                        for n=1:2
                            for s=1:2
                                for t=1:2
                                    qq(m,n,s,t)= Q(i,m)*Q(j,n)*Q(k,s)*Q(l,t);
                                end
                            end
                        end
                    end
                    ss(i,j,k,l)=sum(reshape(qq,[],1).*reshape(X,[],1));
                end
            end
        end
    end
end