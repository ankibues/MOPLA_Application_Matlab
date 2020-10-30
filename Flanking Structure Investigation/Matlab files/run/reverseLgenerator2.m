function L= reverseLgenerator2(Wk)
% this  function generates  the  reverse flow fields for the reverse
% modelling part
%  assuming the original far field flow is general shearing(sinistral sense of shear plus pure shear).


L= zeros(3,3,length(Wk));
for i=1:length(Wk)
    
    if Wk(i)==0
        gamma=0;
        epsilon=1;
        L(:,:,i)= [-epsilon, gamma ,0; 0 , epsilon, 0; 0,0,0];
    else
        gamma=1;
        epsilon= .5*(((1/(Wk(i)^2))-1)^(.5));
        L(:,:,i)= [-epsilon, gamma ,0; 0 , epsilon, 0; 0,0,0];
    end
end

end
