function k= H(x,C)
k= zeros(3,1);
k1= sqrt(x(1)^2 + x(2)^2);
modx= sqrt(x(1)^2 + x(2)^2+ x(3)^2);

n= [x(2); -x(1); 0]*1/k1;
m= [x(1)*x(3); x(2)*x(3); -(x(1)^2 + x(2)^2)]*1/modx/k1;

y= x/modx;

z= @(psi) cos(psi).*n + sin(psi).*m;
zT= @(psi) (z(psi)).*(y') + y.*(z(psi)');
B= @(psi) AK(C,z(psi));
phi= @(psi) Multimix(C,zT(psi));
 for i=1:3
     Eta= @(psi) eta(B(psi),phi(psi),i,y);
     k(i)= integral(@(psi)Eta(psi),0,pi,'ArrayValued',true)*1/4/(pi^2)/modx^2;
 end
 end

function k= AK(C,xi)
a= Multimix(C, xi*xi');
A= [a xi;xi' 0];
k= inv(A);
end

function sum= eta(B,phi,i,y)
sum= 0;
for m= 1:3
    for n= 1:3
        sum= sum + B(i,n)*phi(m,n)*B(4,m) + B(4,4)*B(i,m)*y(m)+ B(4,i)*B(4,m)*y(m);
    end
end


end

function k = Multimix(X,y)
    %Mixing operation as in Mathcad
    % Input: X is a 4th order tensor and y is a second order tensor.
    % Output: The Multimix(X,y)is a second order tensor. 
    k = zeros(3,3);
    for i=1:3
        for j=1:3
            for m=1:3
                for n=1:3
                    k(i,j) = X(i,m,j,n)*y(m,n);
                end
            end
        end
    end
end