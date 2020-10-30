% testing gaussian quadrature 
%{
% check for 1D gaussian quadrature
[point, c] = md_gauss(20,1);
[Alp,ww] = Gauss1(50);  %% this will work for inner integral 
fun= @(x) sin(x).^2;


I1= integral(fun,0,pi);

kn= point(:,1);
kn= pi/2*kn + pi/2;




fun2= @(x) sin(x).^2;

fun3= @(x) sin(x).^2;
k2= fun3(Alp);

k= fun2(kn);
kk2= pi*sum(c.*k)/2;
kkk2=  sum(ww.*k2);
%}

% check for 2D gaussian quadrature
[point, c] = md_gauss(20,2);
[Alp,Bet,ww] = GaussGGQ(20);
fun= @(x,y) sin(x+y).^2;
theta = reshape(Alp,1,[]);
phi   = reshape(Bet,1,[]);
ww= reshape(ww,1,[]);
I1= integral2(fun,0,pi,0,pi);

kn= point(:,1);
kn= pi/2*kn + pi/2;

kn2= point(:,2);
kn2= pi/2*kn2 + pi/2;


fun2= @(x,y) sin(x+y).^2;

k2= fun2(theta,phi);

k= fun2(kn,kn2); 
kk2= pi*pi*sum(c.*k)/4;
kkk2=  pi*pi*sum(ww.*k2)/4;