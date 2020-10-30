function k= Hm(x,C,psi,c)
% This function calculates the Green function Hi for pressure.
% x is the coordinates, (3*1)
% C is stiffness tensor (3*3*3*3)
%output is a 3*1 vector corresponding to Hi.
k= zeros(3,1);
%[psi,c] = Gauss1(30);
[n,~]= size(c);

modx= sqrt(x(1).^2 + x(2).^2 + x(3).^2);


[~,y,z,zT]= ZZone(x,psi);
B= AKm(C,z);   % 4*4*n
PHI= MultimixM(C,zT);

s=zeros(n,1);
for i=1:3
    
    % this part is for calculation of Fi ----------------
    for m=1:3
        for nn =1:3
            s = s + squeeze(B(i,nn,:).*PHI(m,nn,:).*B(4,m,:)) ; 
        end
        s= s + squeeze((B(4,4,:).*B(i,m,:)+ B(4,i,:).*B(4,m,:)).*y(m,:));
    end
    %------------------------
    
    
    kk= -1/4/(pi^2)/modx^2;
    k(i,1)= kk*sum(c.*s); % note : initially i multiply pi/2
    s=zeros(n,1);
end
end


