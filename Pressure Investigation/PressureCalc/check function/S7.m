function s7 = S7(lambda, a, q, q1, omega)
%s7: 3*3*3*3*N
[~,n]  = size(lambda);
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
u= round(u,4);
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
theta  = sum(omega);
v      = zeros(3,3,n);
s7     = zeros(3,3,3,3,n);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                v(k,l,:)=(a(1)*a(2)*a(3))./delta./q.^2.*...
                          omega(i,:).*omega(j,:).*omega(k,:).*omega(l,:).*...
                          (2 - 2.*lambda.*(omega(i,:)+omega(j,:)+omega(k,:)...
                          +omega(l,:))- lambda.*theta + 4.*lambda.*q1./q);
            end
        end
        s7(i,j,:,:,:) = v;
    end
end
end