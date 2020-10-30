
%m= 1:100;

b=0;
%PressvsAnisostrength(m);
n= 2000000;
a=zeros(1,n);
steps = 200;
for i=1:steps
    for j=1:100
        for r= 1:n
            a(1,r)= a(1,r)+ 1;
            for k=1:100
                b= b+ 200;
            end
        end
    end
end

    
    