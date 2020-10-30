function k= fun2check(X)
[~,n]= size(X);
k= zeros(3,3,n);

 for i=1:n
    k(:,:,i)= X(:,i).*X(:,i)';
 end



end