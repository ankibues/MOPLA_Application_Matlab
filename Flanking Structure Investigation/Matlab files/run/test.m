Wk=linspace(0,1,20);
R= logspace(-2,2); 
steps=50;
ic(length(Wk)*length(R))=  struct('WK',[],'viscR',[],'residual',[]);
n=0;
for ii=1:length(Wk)
     for iii=1:length(R)
         n=n+1;
         ic(n).WK= Wk(ii);
         ic(n).viscR= R(iii);
         ic(n).residual=zeros(1,steps);
     end
end
 