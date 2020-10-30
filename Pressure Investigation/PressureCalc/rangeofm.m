m= 1:100;
sqrtm= (m).^0.5;
range= zeros(100,1);
for i=1:100
    range(i)= ((1/m(i))-100)*(48/(26+sqrtm(i)*1000));
end
