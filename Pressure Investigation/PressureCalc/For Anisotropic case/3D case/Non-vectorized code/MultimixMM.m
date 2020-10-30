function k = MultimixMM(X,y)
    %Mixing operation as in Mathcad
    % Input: X is a 4th order tensor (3*3*3*3) matrix and y is(3*3*N1*N2) matrix.
    % N1 and N2 is for number of gaussian points
    
    % Output:  it would be a 3*3*N1*N2 matrix
    [~,~,N1,N2]= size(y);
    k = zeros(3,3,N1,N2);
   for kk=1:N1
       for jj=1:N2
           for i=1:3
               for j=1:3
                   x = reshape(X(i,:,j,:),3,3);
                   k(i,j,kk,jj) = sum(sum(x.*y(:,:,kk,jj)));
               end
           end
       end
   end
   
end