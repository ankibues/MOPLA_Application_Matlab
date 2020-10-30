function [ANG]= generateANG2()

theta1= 0:45:180;
phi1= 0:45:90;

ANG= zeros(3,numel(theta1)*numel(phi1));
m=1;
for i=1:numel(theta1)
    
    for j=1:numel(phi1)
        
                    ANG(:,m)=[theta1(i); phi1(j); 0];
                    m=m+1;
                
            
    end
    
end

               
end
