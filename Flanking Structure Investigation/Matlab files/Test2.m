% testing the velocity field function stage 2
%{
load('points.mat')
[~,n]= size(xxx);
vel_f= zeros(2,n);
 for ir=1:n
        vel_f(:,ir)= Velocity_field_calc(xxx(ir),yyy(ir),a,invS,Jd,d,u2,w,ll);
        disp(ir);
 end
%}

a= [100;1;.999];
y= 1:.05:5;
x= zeros(size(y));

[~,n]= size(y);
lambda=zeros(size(y));
for i=1:n
    lambda(i) = solve_eq(a,[0;x(i);y(i)]);
end
figure
plot(y,lambda)