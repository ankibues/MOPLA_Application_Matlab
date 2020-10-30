function ss=QL(a,q,theta,d,t)
% Input: "a" is a 3*1 matrix, represanting 3 semi-axes;
%        "q" is the transformation tensor Q;
%        "theta" is a 3*3 matrix;
%        "d" is a 3*3 matrix, strain rate tensor for ellipsoid;
%        "t" is time.
%   This function is for Rodrigues rotation approximation.
b = Exp(-theta*t)*q; % this for change in angle 
aa = a.*exp(diag(d)*t); % this one for change in axis length
s = sortrows([b aa],4);
ss = [s(3,:);s(2,:);s(1,:)];
end