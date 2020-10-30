function J = J1(lambda, a)
%J: 3*N
[~,n]  = size(lambda);
u1     = [lambda+a(1).^2; lambda+a(2).^2; lambda+a(3).^2];
u2     = [a(1)^2-a(2)^2; a(1)^2-a(3)^2; a(2)^2-a(3)^2];
u3     = repmat(u2,1,n);
u      = cat(1,u1,u3);
u= round(u,4);
theta  = asin((u(5,:)./u(1,:)).^0.5);  
k      = (u(4,:)./u(5,:)).^0.5;
f      = lellipf(theta, k, 1.0e-10);

e      = lellipe(theta, k, 1.0e-10);
aa     = (f-e).*2./u(4,:)./u(5,:).^0.5;
b      = 2.* u(5,:).^0.5./u(4,:)./u(6,:).*e - 2./u(4,:)./u(5,:).^0.5.*f - ...
         2.* u(3,:).^0.5./u(6,:)./(u(1,:).*u(2,:)).^0.5;
c      = 2./u(6,:).*(u(2,:)./u(1,:)./u(3,:)).^0.5 - 2./u(6,:)./u(5,:).^0.5.*e;
t      = cat(1,aa,b,c);
J      = 2*pi*a(1)*a(2)*a(3).*t;

end