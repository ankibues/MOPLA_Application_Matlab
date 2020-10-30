function Q=Q(ang)

% calculate Q according to three spherical angles;

% Input: ang is the input value indicating the orientation of the inclusion, 
%        which is a set of three spherical angles (1*3 matrix).

% Calculation is based on the book (Jiang,2013,P6) 

% Output: Q is a 3*3 matrix, used to convert the imposed flow from the
%         fixed (C) system to the rotating (C') system.
% for bringing our frame of reference to ellipsoid's frame of reference

    a=[sin(ang(2))*cos(ang(1)),sin(ang(2))*sin(ang(1)),cos(ang(2))];

        if(ang(2)==pi/2)
            
            b=[-sin(ang(1))*sin(ang(3)),cos(ang(1))*sin(ang(3)),cos(ang(3))];
            
        else
            
            b=[cos(atan(tan(ang(2))*cos(ang(1)-ang(3))))*cos(ang(3)),...
                cos(atan(tan(ang(2))*cos(ang(1)-ang(3))))*sin(ang(3)),...
                 -sin(atan(tan(ang(2))*cos(ang(1)-ang(3))))];
             
        end
    
    c=cross(a,b);
    
    Q=[a;b;c];
    
end