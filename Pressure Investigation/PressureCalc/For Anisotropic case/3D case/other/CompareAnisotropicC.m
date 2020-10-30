% to compare the self-consistently generated anisotropic C from the
% Jiang(2016) formula .

[Jd, Js, Ja, Jm, b] = Jnb(); 
load('OrthorhombicHEM.mat')
CC= C_bar_evl(:,:,:,:,100);

Cm    = zeros(3,3,3,3);
m= 10;
Nn=2;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                if rem(i+j,2)==0
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l);
                else
                    Cm(i,j,k,l)= 2*Nn*Jd(i,j,k,l)/m;
                end
            end
        end
    end
end