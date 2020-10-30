% test for vectorized code.
a= [5;3;1];
x=[0;5;0];

[Alpp3,Bett3,www3] = GaussGGQ(40);
theta = reshape(Alpp3,1,[]);
phi   = reshape(Bett3,1,[]);
Wout3 = reshape(www3,1,[]);
[psi,ww] = Gauss1(25);

tic

toc
tic

toc