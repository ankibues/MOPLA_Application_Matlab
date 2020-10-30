%function k= Hm(x,C)
%{
x= [2;3;5];
[Jd, Js, Ja, Jm, b] = Jnb(); 
Cm= 2*Jd;
k= zeros(3,1);
% x is the coordinates, (3*1)
% C is stiffness tensor (3*3*3*3)
[psi,c] = md_gauss(20,1);

  % to change the limits of integration
psi=psi';
[n,~]= size(c);

modx= sqrt(x(1).^2 + x(2).^2 + x(3).^2);
%modx= sqrt((x(1)).^2 + (x(2)).^2 +(x(3)).^2);

[~,y,z,zT]= ZZone(x,psi);
B= AKm(Cm,z);   % 4*4*n
PHI= MultimixM(Cm,zT);

s=zeros(n,1);
for i=1:3
    
    for m=1:3
        for nn =1:3
            s = s + squeeze(B(i,nn,:).*PHI(m,nn,:).*B(4,m,:)) ; 
        end
        s= s + squeeze((B(4,4,:).*B(i,m,:)+ B(4,i,:).*B(4,m,:)).*y(m,:));
    end
    kk= -1/4/(pi^2)/modx^2;
    k(i,1)= kk*pi*sum(c.*s)/2;
    s=zeros(n,1);
end
%}
a= [3;1.5;1];
x= [0;3;0];
[Jd, Js, Ja, Jm, b] = Jnb(); 
Cm    = zeros(3,3,3,3);
m= 1;
Nn=1;
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
 ang = [0; 0;(pi/180)*0];                       
    q = Q(ang);
 Cm_a    = Transform(Cm, q);
[Alp1,Bet1,ww1] = GaussGGQ(20);
theta = reshape(Alp1,1,[]);
phi   = reshape(Bet1,1,[]);
ww2 = reshape(ww1,1,[]);
 
[psi1,wwww1] = Gauss1(20);

tic
p2= LambdaExtAniso3DPress(a,x,Cm_a,theta,phi,ww2,psi1,wwww1);
toc


