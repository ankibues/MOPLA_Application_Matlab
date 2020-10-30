function k= AKm_vec(C,xi)
% C is a 3*3*3*3 matrix, xi will be a 3*N1*N2 matrix
% Output k will be a 4*4*N1*N2
[~,N1,N2]= size(xi);
k=zeros(4,4,N1,N2);
a= MultimixMMvec(C, fun2(xi)); % 3*3*N1*N2
for i=1:N1
    for j=1:N2
        p = [ a(:,:,i,j) squeeze(xi(:,i,j)); squeeze(xi(:,i,j))'  0];
        k(:,:,i,j)= inv(p);  % this can be speeded up !
    end
end

end

function kk= fun2(X)  % this makes xi as 3*N1*N2 to 3*3*N1*N2
[~,N1,N2]= size(X);
kk= zeros(3,3,N1,N2);
 for i=1:N1
    for j=1:N2
        kk(:,:,i,j)= squeeze(X(:,i,j)).*squeeze(X(:,i,j))';
    end
 end
 
end