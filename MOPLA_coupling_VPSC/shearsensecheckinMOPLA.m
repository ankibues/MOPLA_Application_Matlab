
load('MOPLAandVPSCstuff.mat','Q_evl','A_evl','C_bar_evl','M_bar_evl','LL');
%for 1st quartz element
LLL     = [0 2 0 ;0 0 0;0 0 0];
WWW = 0.5 * (LLL - LLL');
WWWvec= [2*WWW(3,2); 2*WWW(1,3); 2*WWW(2,1)];
zaxis= [0;0;1];
dotprod= dot(WWWvec,zaxis);
L=zeros(3,3,5,500);
W=zeros(3,3,5,500);
dotprod1=zeros(1,5,500);
Wvec = zeros(3,1,5,500);
for i=1:5
    for j=1:500
        Q = Q_evl(:,:,i,j);
        L(:,:,i,j)   = Q'*LL(:,:,i)*Q;
        W(:,:,i,j) = 0.5 * (L(:,:,i,j) - L(:,:,i,j)');
        Wvec(:,1,i,j)= [2*W(3,2,i,j); 2*W(1,3,i,j); 2*W(2,1,i,j)];
        dotprod1(:,i,j)= dot(Wvec(:,1,i,j),zaxis);
    end
end
