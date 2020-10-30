function s5 = S5(lambda, a, q, omega)
%s5: 3*N 
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
u= round(u,4);
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
t      = (a(1)*a(2)*a(3))./delta./q;
t1     = repmat(t,3,1);
t2     = (repmat(lambda,3,1).*omega-0.5);
s5     = t1.*t2;
end