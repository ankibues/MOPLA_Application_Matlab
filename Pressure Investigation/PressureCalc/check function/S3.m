function s3 = S3(lambda, a, q, omega)
%s3: 3*N
u      = [a(1)^2+lambda; a(2)^2+lambda; a(3)^2+lambda];
u =round(u,4);
delta  = (u(1,:).*u(2,:).*u(3,:)).^0.5;
t      = (a(1)*a(2)*a(3)).*lambda./delta./q; %t:1*N
tt     = repmat(t,3,1); 
s3     = tt.*omega;
end