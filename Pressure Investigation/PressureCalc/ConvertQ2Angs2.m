function [ang1, ang2, ang3] = ConvertQ2Angs2(orien)
% ConvertQ2Angs.m
% Convert tranformation matrix Q to Theta[0,2pi] and Phi[0,pi] (i.e., cover 
% the whole range of 0<=Phi<=pi) of a1, a2, a3.
%
% Input:  qq,    3*3 matrix;

% Output: ang1,   theta, phi of a1, in radian, 2*n matrix;
%         ang2,   theta, phi of a2, in radian, 2*n matrix;
%         ang3,   theta, phi of a3, in radian, 2*n matrix;
%--------------------------------------------------------------------------
  
  aa1  = orien(1,:)';
  aa2  = orien(2,:)';
  aa3  = orien(3,:)';
  
  ang1 = orienMat_whole(aa1);
  ang2 = orienMat_whole(aa2);
  ang3 = orienMat_whole(aa3);
end

function angs = orienMat_whole(u)

% Input:  u, 3*1 matrix, each col stands for a orientation vector

% Output: angs, 2*1 matrix, change to theta and phi

  
  angs = zeros(2,1);
  
  
  if u(1)>0
      angs(1) = atan(u(2)/u(1));
  elseif u(1)==0
      angs(1) = 0.5*pi*sign(u(2));
  elseif u(1)<0
      angs(1)= pi + atan(u(2)/u(1));
  end
  
  angs(2) = acos(u(3));
end