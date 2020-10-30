
aa = [1;1;1/3];
x = [0.5^0.5;0;1/3];

a = zeros(3,3);
J = zeros(3,3,3,3);
J1 = zeros(3,3,3,3);
delta = zeros(3,3,3,3);
delta(:,:,1,1) = diag([1 1 1]);
delta(:,:,2,2) = diag([1 1 1]);
delta(:,:,3,3) = diag([1 1 1]);
for i=1:3
    for j=1:3
        a(i,j)=1;
        J(i,j,:,:)=a;
        a = zeros(3,3); 
    end
end
for i=1:3
    for j=1:3
        J1(:,:,j,i)=J(:,:,i,j);
    end
end
Js = 0.5*(J+J1);
Ja = 0.5*(J-J1);
Jm = 1/3*delta;
Jd = Js-Jm;
[Alp,Bet,ww] = Gauss(1000);
%[s,pi,p] = SPpIn(aa,Jd,Js,Ja); 
%[s,pi,p] = SnPI(aa,x,Jd,Js,Ja);
gp = 2702;
%p = pEx('apgq',aa,x)
p = pEx('glesh',aa,x,Alp,Bet,ww);
%choice:6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,
%  770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810
%p = pEx('ggq',aa,x,gp)