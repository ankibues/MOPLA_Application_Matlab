%testing 
tic
a= [5;3;1];
x= [2;3;5];
[Jd, Js, Ja, Jm, b] = Jnb(); 
Cm= 2*Jd;
%{
[psi,c] = md_gauss(20,1);
[p, w]            = Gauss(100);
p=p';
  % to change the limits of integration
psi=psi*pi/2 +pi/2;
psi=psi';
  tic
  k= Hm(x,Cm,psi,c);
  toc
  
  tic
  k1= Hm(x,Cm,p,w);
  toc
  %}
[Alp,Bet,ww] = GaussGGQ(20);
theta1 = reshape(Alp,1,[]);
phi1   = reshape(Bet,1,[]);
Wout1 = reshape(ww,1,[]);
[psi1,Win1] = Gauss1(20);

p1 = LambdaExtAniso3DPress(a,x,Cm,theta1,phi1,Wout1,psi1,Win1);
toc

%{

    gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 
%%%%%%%%%%%%%%%Test of Green tensor calculation for Anisotropic material

% Eshelby Tensors (S_el,PI_el) for Interior points using elliptical integrals  
  [S_el,PI_el] = SnP(a);

% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green
% tensor for pressure)
% p
  pp   = zeros(3,3);
  for j=1:3
      pp(j,j) = -1/3* (S_el(j,j,1,1)+ S_el(j,j,2,2)+ S_el(j,j,3,3)); 
  end
% S  
  S       = S_el;
  for k=1:3
      for l=1:3
          S(k,k,l,l) = pp(k,k)+ S_el(k,k,l,l);
      end
  end
  
  % S  
  S1       = S_el;
  for k=1:3
      for l=1:3
          S1(k,k,l,l) = pp(k,k)+ S_el(k,k,l,l);
      end
  end
%{  
% PI  
lambda = SnpInAniso3DPress(a,Cm);
Carray  = C2OneDarray(Cm);
    T= TGreen(a,  Carray, Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
    S = Contract(Jd, Contract(T,Cm));
    
  %}
 %}
%{
[theta2, phi2, Wout2] = Lebedev(2702) ;
theta2 = reshape(theta2,1,[]);
phi2   = reshape(phi2,1,[]);
Wout2 = reshape(Wout2,1,[]);
tic
p2 = LambdaExtAniso3DPress(a,x,Cm,theta2,phi2,Wout2,psi1,Win1);
toc  
%}