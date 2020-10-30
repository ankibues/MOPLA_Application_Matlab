% test 
clear;
clc;
[Jd, Js, Ja, Jm, b] = Jnb();

L = [1,0,0;0,-1,0;0,0,0];                                                  % bulk flow
D = 0.5*(L+L');
W = 0.5*(L-L');

Nm = 3;                                                                    % stress exponents for the matrix and clast
Nc = 3;

ang = [0;90;90];
ang = degtorad(ang);
q = Q(ang);
x(1:3,1:3) = q;
x(1:3,4) = [1;1/3;1/3];

x(1:3,5) = [1000;0.5;0];
x(1:3,6) = [0;0;0];

tincr = 0.0125;

[h,u2,Pin] = InP(x,D,Nm,Nc,Jd,Js,b);
a1 = x(1,4);
a2 = x(2,4);
a3 = x(3,4);
range = 1;
y1range = range*a1;
y2range = range*a2;
y3range = range*a3;

y1 = -y1range:0.1:y1range;
y2 = -y2range:0.1:y2range;
y3 = -y3range:0.1:y3range;
intcr = 0.1;

gp = 50;
[Alp,Bet,ww] = GaussGGQ(gp);

P_ab = Pab(y1, y2, a1, a2, x, h, u2, Pin, Alp, Bet, ww);
%P_ac = Pac(y1, y3, a1, a3, x, h, u2, Pin);
%P_bc = Pbc(y2, y3, a2, a3, x, h, u2, Pin);

%[x,y,z,P_abc] = Pabc(y1range, y2range, y3range, intcr, a1, a2, a3, x, h, u2, Pin);


TwoDContour(P_ab)


%TwoDContour(P_ac)


%TwoDContour(P_bc)

%ThreeDContour(x,y,z,P_abc)
%[faces,verts,colors] = isosurface(x,y,z,P_abc,Pin,P_abc);
%patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors,...
    %'Facecolor','interp','edgecolor','interp');
%colormap summer
%view(3)
