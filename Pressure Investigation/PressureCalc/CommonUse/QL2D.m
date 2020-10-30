function ss=QL2D(a,q,theta,d,t)
% Input: "a" is a 2*1 matrix, represanting 2 semi-axes;
%        "q" is the transformation tensor Q;
%        "theta" is a 2*2 matrix;
%        "d" is a 2*2 matrix, strain rate tensor for ellipsoid;
%        "t" is time.
%   This function is for Rodrigues rotation approximation.
b = Exp2D(-theta*t)*q; % this for change in angle 
aa = a.*exp(diag(d)*t); % this one for change in axis length
s = sortrows([b aa],3);
ss = [s(2,:);s(1,:)];
end