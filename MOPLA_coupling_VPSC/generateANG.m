function [ANG]= generateANG()

% generate angles for triaxial ellipsoids

theta1= 0:45:180;
theta1(end)=[];
phi1= 0:45:90;
theta2= 0:45:180;
theta2(5)=[];

ANG= zeros(3,numel(theta1)*numel(phi1));
m=1;
for i=1:numel(theta1)
    
    for j=1:numel(phi1)
        
        for k=1:numel(theta2)
        
                    ANG(:,m)=[theta1(i); phi1(j); theta2(k)];
                    m=m+1;
        end
        
            
    end
    
end


end
