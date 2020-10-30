function wp = Wd(a,m,d)
% Ellipsoid vorticity referred to the frame tracking its semi-axes

% Input: a is three semi-axes, 3*1 matrix
%        m is the W in the ellipsoid coordinate, 3*3 matrix;
%        d is the new strain rate tensor, 3*3 matrix;

% Output: wp updated vorticity, 3*3 matrix;

 
  wp = zeros(3,3);          % wp is simply the shear spin part of rotation of ellipsoid.
  
  for i = 1:3
      for j = 1:3
          if i == j         % find the reason why?
              wp(i,j) = 0;
          elseif a(i) == a(j)
              wp(i,j) = m(i,j);
          else
              wp(i,j) = p(a(i),a(j))*d(i,j);
          end
      end
  end
  
end

function r = p(x,y)
  
r = (x^2+y^2)/(x^2-y^2);

end
