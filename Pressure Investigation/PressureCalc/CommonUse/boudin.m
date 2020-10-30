function ss=boudin(a,q,theta,d,t)
% Input: "a" is a 3*1 matrix, represanting 3 semi-axes;
%        "q" is the transform tensor Q;
%        "theta" is a 3*3 matrix;
%        "d" is a 3*3 matrix;
%        "t" is time.

s = QL(a,q,theta,d,t);

if s(1,4)/s(3,4)>100 % try to understand why we are doing this ?It is somewhat related to treatment of ellongated ellipsoid
    s(1,4)=0.5*s(1,4);
end

if s(2,4)/s(3,4)>100
    s(2,4)=0.5*s(2,4);
end
s = sortrows(s,4);
ss = [s(3,:);s(2,:);s(1,:)];
end