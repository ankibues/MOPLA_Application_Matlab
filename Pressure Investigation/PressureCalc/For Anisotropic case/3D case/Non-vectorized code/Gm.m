function k= Gm(x,C,psi,c)
k= zeros(3,3,3,1);
% x is the coordinates, (3*1)
% C is stiffness tensor (3*3*3*3)

[n,~]= size(c);

modx= sqrt(x(1).^2 + x(2).^2 + x(3).^2);


[~,y,z,zT]= ZZone(x,psi);
B= AKm(C,z);   % 4*4*n
PHI= MultimixM(C,zT);

s=zeros(n,1);
for i=1:3
    for j=1:3
        for l=1:3
    
  % this part is for calculation of Mij  ----------------        
    for m=1:3
        for nn =1:3
            s = s + squeeze(B(i,m,:).*B(j,nn,:).*PHI(m,nn,:)) ; 
        end
        s= s + squeeze((B(i,m,:).*B(4,j,:)+ B(j,m,:).*B(4,i,:)).*y(m,:));
    end
  %----------------------------------    
  % here, s = Mij, now the integrand I becomes:
  
    I= (-y(l,:).*squeeze(B(i,j,:))) + (z(l,:)'.*s);
    kk= 1/4/(pi^2)/modx^2;
    k(i,j,l,1)= kk*sum(c.*I);
    s=zeros(n,1);
    
        end
    end
    
end
end


