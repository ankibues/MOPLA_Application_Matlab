% script to plot vorticity variation for randomly generated shapes and
% orientations of RDEs and different viscosities .

Wk=.75;% [1];%,.5, .75, .90];
[a, ang]=RandAANG(5,50); 

r= [5];% 2 5 10];
n=1;
    
        for k=1:length(a)
            for l=1:length(r)
                n=n+1;
                plotVorticity(Wk,r(l),a(:,k),ang(:,k),6)
                hold on
            end
        end

    
%saveas(gcf,'test50RDEsWk=.7.fig')