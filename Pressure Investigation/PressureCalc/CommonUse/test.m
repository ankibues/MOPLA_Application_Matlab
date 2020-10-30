% Validation of SnP_poisson function when poisson's ratio nu is not equal
% .5
nu=.5;
a=[5;4;1];


[S1,~]=Eshint(nu,a);
[S2,~] = SnP_poisson(a,nu);

S3= [S2(1,1,1,1),  0,      0,      S2(1,1,2,2),  0,      S2(1,1,3,3);
          0,      2*S2(1,2,1,2),  0,      0,      0,      0;
          0,      0,      2*S2(1,3,1,3),  0,      0,      0;
          S2(2,2,1,1),  0,      0,      S2(2,2,2,2),  0,      S2(2,2,3,3);
          0,      0,      0,      0,      2*S2(2,3,2,3),  0;
          S2(3,3,1,1),  0,      0,      S2(3,3,2,2),  0,      S2(3,3,3,3)];
      
    error= S3-S1;
      [S0,~] = SnP(a);