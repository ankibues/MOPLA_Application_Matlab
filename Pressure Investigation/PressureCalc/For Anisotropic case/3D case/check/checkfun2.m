function kk= checkfun2(X)
[~,n]= size(X);
kk= zeros(3,3,n);

 for i=1:n
    kk(:,:,i)= X(:,i).*X(:,i)';
 end



end