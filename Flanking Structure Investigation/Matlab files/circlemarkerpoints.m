function KKK1 = circlemarkerpoints()
      
        r=0.1;
        eta= 0:.1:2*pi;
        x1=[0 .3 .6 .9 1.2];
        y1=[.25 .50 .75 -.25 -.50 -.75];
       point(length(x1)*length(y1))=  struct('X1',[],'Y1',[],'pts',[]);
        n=0;
        for i=1:length(x1)
            for j=1:length(y1)
                n=n+1;
                xgv1     = x1(i) + r*cos(eta);     %-8:.05:8; %:.1:8;           % grid vector: x'axis,a1
                ygv1     = y1(j) + r*sin(eta);           % grid vector: y'axis,a2
                zgv1     = zeros(size(xgv1)); %-5:.01:5;  % grid vector: z'axis,a3   ----> considering this 0 for now, since only consider 2D here !
                point(n).X1= x1(i);
                point(n).Y1=y1(i);
                point(n).pts=[xgv1;ygv1;zgv1];
            end
        end
        
                
%KKK= zeros(3,length(x1)*length(y1)*length(eta));

KKK1= point(1).pts;
for ii=1:(length(point)-1)
       
    KKK= [KKK1 point(ii+1).pts];
    KKK1=KKK;
end

end
