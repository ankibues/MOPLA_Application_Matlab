function kk= check2fun2(X)  % this makes xi as 3*N1*N2 to 3*3*N1*N2
[~,N1,N2]= size(X);
kk= zeros(3,3,N1,N2);
 for i=1:N1
    for j=1:N2
        kk(:,:,i,j)= squeeze(X(:,i,j)).*squeeze(X(:,i,j))';
    end
 end
 
end