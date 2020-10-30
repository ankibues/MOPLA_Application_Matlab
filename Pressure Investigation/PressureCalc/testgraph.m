x=-2:1:2;
y=-2:1:3;
[X,Y]=meshgrid(x,y);
Z= X.*exp(-X.^2 - Y.^2);
figure
contourf(X,Y,Z)
