% script to plot r variation for randomly generated shapes and
% orientations of RDEs and different viscosities .

Wk=1;% [1];%,.5, .75, .90];
[a, ang]=RandAANG(5,50); 

r= 2;% 3 4];
%n=1;
    for j= 1:length(ang)
        %for l=1:length(r)
            %for k=1:length(a)
                %n=n+1;
                Variationinr(a(:,j), r,ang(:,j));
                hold on
            %end
        %end 
    end
%saveas(gcf,'Wk=.70test_15RDEs.fig')