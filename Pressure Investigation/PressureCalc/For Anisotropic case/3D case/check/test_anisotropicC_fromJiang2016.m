ang=0; 
Nn=2; 
 m=10;
 [Jd, Js, Ja, Jm, b] = Jnb(); 
%{
 [Jd, Js, Ja, Jm] = FourIdentity2D();
 Cm    = zeros(2,2,2,2);
Cm1    = zeros(2,2,2,2);
for i=1:2
    for j=1:2
        for k=1:2
            for l=1:2
                if i==j
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm1(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end

        

Cm(1,1,1,1)= Nn*((cos(2*ang))^2 + ((sin(2*ang))^2)/m);
Cm(2,2,2,2)=  Cm(1,1,1,1);
Cm(1,1,2,2)= -Cm(1,1,1,1);
Cm(1,1,1,2)= Nn*((1/m)-1)*sin(2*ang)*cos(2*ang);
Cm(1,2,2,2)= -Cm(1,1,1,2);
Cm(1,2,1,2)= Nn*((sin(2*ang))^2 + ((cos(2*ang))^2)/m);
Cm(2,2,1,2)= Cm(1,2,2,2);
Cm(2,2,1,1)= Cm(1,1,2,2);
Cm(1,2,1,1)= Cm(1,1,1,2);
Cm(2,1,1,1)= Cm(1,2,1,1);
Cm(1,1,2,1)= Cm(1,1,1,2);
Cm(2,1,2,2)= Cm(1,2,2,2);
Cm(2,2,2,1)= Cm(2,2,1,2);
Cm(2,1,1,2)= Cm(1,2,1,2);
Cm(2,1,2,1)= Cm(2,1,1,2);
Cm(1,2,2,1)= Cm(2,1,2,1);
%}
 Cm    = zeros(3,3,3,3);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end