%fr = funIn(sub,a,x,alp,bet,gamma,Cm)
a=[5;3;1];
sub=[1,2];
   gp = 20;
  [point, w] = md_gauss(gp,3);
   alp= point(:,3)';
   bet= point(:,1)';
   gamma= point(:,2)';
theta= pi*alp + pi;   % changing the limits   1*n
phi= (pi/2)*bet + pi/2;
psi= (pi/2)*gamma+ pi/2;
xp= [a(1).*cos(theta).*sin(phi); a(2).*sin(theta).*sin(phi); a(3).*cos(phi)];  %3*n
xi= [ cos(theta).*sin(phi); sin(theta).*sin(phi);cos(phi)];
[modx,y,z,zT]= ZZ(x,xp,psi);         % zT is 3*3*n, z is 3*n, y is 3*n, modx is 3*n, n here is number of gaussian points
B= AKm(C,z);   % B is 4*4*n
PHI= MultimixM(C,zT); %PHI is 3*3*n
[~,~,N]= size(PHI);

s=zeros(1,N);

    for m=1:3
        for nn =1:3
            s = s + (squeeze(B(sub(1),nn,:).*PHI(m,nn,:).*B(4,m,:)))' ; 
        end
        s= s + (squeeze((B(4,4,:).*B(sub(1),m,:)))' + squeeze(B(4,sub(1),:).*B(4,m,:))').*y(m,:);
    end
kkk= a(1)*a(2)*a(3)/4/(pi^2);
fr= kkk.*s.*xi(sub(2),:).*sin(phi)./modx^2./a(sub(2));
 k   = w'.*fr;
  T = -(2*pi-0)*(pi-0)*(pi-0)*sum(k)/8;
