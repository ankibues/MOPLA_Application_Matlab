function k= AKm(C,xi)
% C is a 3*3*3*3 matrix, xi will be a 3*n matrix
% Output k will be a 4*4*N
[~,n]= size(xi);
k=zeros(4,4,n);
a= MultimixM(C, fun2(xi)); % 3*3*n

for i=1:n
    p= [a(:,:,i) xi(:,i); xi(:,i)'  0];
    k(:,:,i)= inv(p);  % this can be speeded up !
end

end

function kk= fun2(X)
[~,n]= size(X);
kk= zeros(3,3,n);

 for i=1:n
    kk(:,:,i)= X(:,i).*X(:,i)';
 end



end
