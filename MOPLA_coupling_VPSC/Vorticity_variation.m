% script to plot vorticity variation for randomly generated shapes and
% orientations of RDEs and different viscosities .

Wk=1;% [1];%,.5, .75, .90];
[a, ang]=RandAANG(5,15); 

r= [2 ];
%n=1;
    for j= 1:length(ang)
        for l=1:length(r)
            for k=1:length(a)
                %n=n+1;
                plotVorticity(Wk,r(l),a(:,k),ang(:,k),6)
                hold on
            end
        end 
    end
%saveas(gcf,'Wk=.70test_15RDEs.fig')