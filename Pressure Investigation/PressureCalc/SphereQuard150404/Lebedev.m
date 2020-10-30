function [theta, phi, w] = Lebedev(m)
% m is the the number of pionts generated 24/03/2015

    if m == 6
    [x,y,z,w] = ld0006;
    elseif m == 14
        [x,y,z,w] = ld0014;    
    elseif m == 26
        [x,y,z,w] = ld0026;
    elseif m == 38
        [x,y,z,w] = ld0038;
    elseif m == 50
        [x,y,z,w] = ld0050;
    elseif m == 74
        [x,y,z,w] = ld0074;
    elseif m == 86
        [x,y,z,w] = ld0086;
    elseif m == 110
        [x,y,z,w] = ld0110;
    elseif m == 146
        [x,y,z,w] = ld0146;
    elseif m == 170
        [x,y,z,w] = ld0170;
    elseif m == 194
        [x,y,z,w] = ld0194;
    elseif m == 230
        [x,y,z,w] = ld0230;
    elseif m == 266
        [x,y,z,w] = ld0266;
    elseif m == 302
        [x,y,z,w] = ld0302;
    elseif m == 350
        [x,y,z,w] = ld0350;
    elseif m == 434
        [x,y,z,w] = ld0434;
    elseif m == 590
        [x,y,z,w] = ld0590;
    elseif m == 770
        [x,y,z,w] = ld0770;
    elseif m == 974
        [x,y,z,w] = ld0974;
    elseif m == 1202
        [x,y,z,w] = ld1202;    
    elseif m == 1454
        [x,y,z,w] = ld1454;   
    elseif m == 1730
        [x,y,z,w] = ld1730;   
    elseif m == 2030
        [x,y,z,w] = ld2030;
    elseif m == 2354
        [x,y,z,w] = ld2354;
    elseif m == 2702
        [x,y,z,w] = ld2702;
    elseif m == 3074
        [x,y,z,w] = ld3074;
    elseif m == 3470
        [x,y,z,w] = ld3470;
    elseif m == 3890
        [x,y,z,w] = ld3890;
    elseif m == 4334
        [x,y,z,w] = ld4334;
    elseif m == 4802
        [x,y,z,w] = ld4802;
    elseif m == 5294
        [x,y,z,w] = ld5294;
    elseif m == 5810
        [x,y,z,w] = ld5810;
    else
    fprintf ( 'Hey! Please make a proper choice! ' );
    fprintf ( 'I LOVE SCIENCE XDXDXD ');
    fprintf ( 'The time to work on a problem is after you have solved it.');
    end
    
    [theta, phi] = xyz2thetaphi(x,y,z);
    
    
end