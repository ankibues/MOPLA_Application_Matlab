function s2 = S2(lambda, a)
%s2: 3*3*N
[~,n]  = size(lambda);
J      = J1(lambda, a);
JJ     = J2(lambda, a);
s2     = zeros(3,3,n);
for i=1:3
    for j=1:3
      s2(i,j,:) = (1/4/pi).*((J(i,:)-J(j,:))./2+a(i)^2.*(squeeze(JJ(i,j,:)))');
    end
end
end