function [y,z,zT]= ZZ_vec(x,psi)
% x is 3*N1 matrix representing coordinates 

% z will be 3*N1*N2 matrix. zT is 3*3*N1*N2 matrix
[~,N1]= size(x);   
[~,N2]=size(psi);
k= sqrt(x(1,:).^2 + x(2,:).^2); %1*N1 matrix
%k= sqrt((x(1)).^2 + (x(2)).^2);
%k= repmat(k,3,1);
%modx= sqrt((x(1)).^2 + (x(2)).^2+(x(3)).^2);
modx= sqrt(x(1,:).^2 + x(2,:).^2 + x(3,:).^2);    %1*N1
%modx= repmat(modx,3,1);

%n= [x(2); -(x(1)); 0]./k;
n= [x(2,:); -x(1,:); zeros(1,N1)]./repmat(k,3,1); %3*N1


%m= [(x(1)).*(x(3)); (x(2)).*(x(3)); -((x(1))^2 + (x(2))^2)]/k/modx; 

m= [x(1,:).*x(3,:); x(2,:).*x(3,:); -(x(1,:).^2 + x(2,:).^2)]./repmat(k,3,1)./repmat(modx,3,1); %3*N1


y= (x)./repmat(modx,3,1);    %3*N1
nn= repmat(n,1,1,N2);
mm= repmat(m,1,1,N2);
% cos(psi)--> 1*N2
psi= reshape(psi,1,1,N2);

z= repmat(cos(psi),3,N1,1).*nn + repmat(sin(psi),3,N1,1).*mm; %3*N1*N2

zT= fun1(z,y) + fun2(y,z); % 3*3*N1*N2 matrix

end

function k= fun1(x,y)
% x 3*N1*N2
% y 3*N1

[~,~,N2]= size(x);
[~,N1]=size(y);
k= zeros(3,3,N1,N2);
for j=1:N1
    for i=1:N2
    k(:,:,j,i)= squeeze(x(:,j,i)).*y(:,j)';
    end
end
end

function k= fun2(x,y)
% x 3*N1
% y 3*N1*N2

[~,N1]= size(x);
[~,~,N2]=size(y);
k= zeros(3,3,N1,N2);
for j=1:N1
    for i=1:N2
    k(:,:,j,i)= x(:,j).*squeeze(y(:,j,i))';
    end
end
end