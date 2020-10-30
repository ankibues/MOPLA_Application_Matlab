function k = MultimixM(X,y)
    %Mixing operation 
    % Input: X is a 4th order tensor (3*3*3*3) matrix and y is a second order tensor(3*3*n) matrix.
    % n is for number of gaussian points
    
    % Output:  a 3*3*n matrix
    [~,~,N]= size(y);
    k = zeros(3,3,N);
   for kk=1:N
    for i=1:3
        for j=1:3
            x = reshape(X(i,:,j,:),3,3);
            k(i,j,kk) = sum(sum(x.*y(:,:,kk)));
                
        end
    end
   end
   
end
