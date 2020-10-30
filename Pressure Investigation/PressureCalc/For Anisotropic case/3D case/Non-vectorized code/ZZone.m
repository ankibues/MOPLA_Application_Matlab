function [modx,y,z,zT]= ZZone(x,psi)

% x is 3*1 matrix representing coordinates 
% z will be 3*n matrix. zT is 3*3*n matrix


k= sqrt((x(1)).^2 + (x(2)).^2);

modx= sqrt(x(1).^2 + x(2).^2 + x(3).^2); % 1*1 matrix



n= [x(2); -(x(1)); 0]./k;  % 3*1



m= [(x(1)).*(x(3)); (x(2)).*(x(3)); -((x(1))^2 + (x(2))^2)]/k/modx; %3*1


y= (x)./modx;  % 3*1


z= repmat(cos(psi'),3,1).*n + repmat(sin(psi'),3,1).*m; %3*n

zT= fun1(z,y) + fun2(y,z); % 3*3*n matrix

end


function k1= fun1(x,y)
% x  3*n matrix
% y  3*1 matrix
[~,n]= size(x);
k1= zeros(3,3,n);
    for i=1:n
        k1(:,:,i)= x(:,i).*(y');
    end
end
%used for testing

function k= fun2(x,y)
% x  3*1 matrix
% y  3*n matrix
[~,n]= size(y);
k= zeros(3,3,n);
 for i=1:n
    k(:,:,i)= x.*y(:,i)';
 end
end   %used for testing
%}

function k= fun3(x,y)
% x and y are both 3*n matrix
[~,n]= size(y);
k= zeros(3,3,n);
 for i=1:n
    k(:,:,i)= x(:,i).*y(:,i)';
 end
end