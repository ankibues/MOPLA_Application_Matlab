L     = [0 1;0 0]; 
E = 0.5 * (L + L');
n= 0:.2:4*pi;

alphas= zeros(numel(n),1);
for i= 1:numel(n)
    q= [cos(n(i)), sin(n(i)); -sin(n(i)), cos(n(i))];
    Nm=10;
    e= q*E*q';
    sigma= 2*Nm*E;
    sig= q*sigma*q';
    sigma_inva= inva(sigma);
    r= 15;
    R=5;
    Pdiff= (2*Nm*e(1,1)*(1-r)*((R^2)-1))/((R^2)+2*r*R+1);
    alpha= Pdiff/sigma_inva;
    alphas(i,1)= alpha;
    
end

    