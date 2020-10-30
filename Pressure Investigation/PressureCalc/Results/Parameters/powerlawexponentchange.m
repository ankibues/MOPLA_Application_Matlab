% varying power law exponent
% plotting figure for paper
Nm = 1:7;
K = zeros(1,7);
J = zeros(1,7);
a1= [10;5;1];
a2= [10;9.9999;1];
%a3= [10;.9999;1];
a4= [100;99.9999;1];
%a5= [100;.9999;1];

A= [a1 a2 a4];

 for j= 1:3;
    for i=1:7
        K(i)= PasfuncN(A(:,j),i,100);
        J(i)= PasfuncN(A(:,j),i,.01);
    end
    
    plot(Nm,K)
    xlabel('Power law stress exponent');
    hold on
    plot(Nm,J)
    xlabel('Power law stress exponent');
    hold on
    K = zeros(1,7);
    J = zeros(1,7);
 end






