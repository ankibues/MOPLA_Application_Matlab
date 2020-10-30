function k = MultimixMvec(X,y)
    %Mixing operation as in Mathcad
    % Input: X is a 4th order tensor (3*3*3*3) matrix and y is a second order tensor(3*3*n) matrix.
    % n is for number of gaussian points
    
    % Output: The Multimix(X,y)should be a second order tensor. Here, it would be a 3*3*n matrix
    [~,~,N]= size(y);
    k = zeros(3,3,N);
   
    for i=1:3
        for j=1:3
            x = repmat(reshape(X(i,:,j,:),3,3),1,1,N);
            k(i,j,:) = sum(sum(x.*y));      
        end
    end
   end
   

