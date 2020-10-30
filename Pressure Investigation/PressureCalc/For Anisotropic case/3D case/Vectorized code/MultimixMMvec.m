function k = MultimixMMvec(X,y)
    % Performs a mixing operation 
    % Input: X is a 4th order stiffness tensor (3*3*3*3)of the matrix and y is a (3*3*N1*N2) matrix.
    % N1, N2 are number of gaussian points corresponding to theta/phi and
    % psi respectively.
    
    % Output:   3*3*N1*N2 matrix
     [~,~,N1,N2]= size(y);
    k = zeros(3,3,N1,N2);
   
    for i=1:3
        for j=1:3
            x = repmat(reshape(X(i,:,j,:),3,3),1,1,N1,N2);
            k(i,j,:,:) = sum(sum(x.*y));      
        end
    end
   end
   

