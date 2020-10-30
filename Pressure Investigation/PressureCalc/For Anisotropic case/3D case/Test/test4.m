% Test of MultimixM vs MultimixMvec

a= [5;3;1];
x= [6;4;2];
[Jd, Js, Ja, Jm, b] = Jnb(); 
Cm= 2*Jd;

yy= rand(3,3,400);

tic
k1 = MultimixM(Cm,yy);
toc

tic
k2 = MultimixMvec(Cm,yy);
toc

tic 
k4 = MultimixMM(Cm,yy);
toc

tic 
k3 = MultimixMMvec(Cm,yy);
toc
