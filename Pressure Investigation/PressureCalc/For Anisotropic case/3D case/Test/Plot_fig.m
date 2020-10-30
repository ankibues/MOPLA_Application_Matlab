%{
a1= [5;2;1];
a2= [5;3;1];
a3=[5;4;1];

figure
Plot_ep_X(a1)
hold on
Plot_ep_X(a2)
hold on
Plot_ep_X(a3)
%}
gp=10:10:150;

for i=1:length(gp)
    plotTvsGP(gp(i))
end