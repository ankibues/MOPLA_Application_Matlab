function kk= PHICALC(C,z,x)
% C is 3*3*3*3 matrix
% z is 3*n matrix
% x is 3*1 matrix
[~,nn]= size(z);
kk= zeros(3,3,nn);
s=0;
for i=1:3
        for j=1:3
            for m=1:3
                for n=1:3
                    s = s+ C(i,m,j,n)*(z(m,:)*x(n) + z(n,:)*x(m));
                end
            end
            kk(i,j,:)=s;
            s=0;
        end
end

end
