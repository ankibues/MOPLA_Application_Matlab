function [Jd, Js, Ja, Jm, b] = Jnb()
% get Jd,Js,Ja,Jm,b for following calculations 28/03/2015
% J : Fourth order Identity tensor and b for calculating Inverse of Fourth
% order tensor
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

    b(1,:,:)=1/6^0.5*(-eye(3)+3*[0 0 0;0 0 0;0 0 1]);
    b(2,:,:)=1/2^0.5*(-[1 0 0;0 0 0;0 0 0]+[0 0 0;0 1 0;0 0 0]);
    b(3,:,:)=1/2^0.5*([0 0 0;0 0 1;0 0 0]+[0 0 0;0 0 0;0 1 0]);
    b(4,:,:)=1/2^0.5*([0 0 1;0 0 0;0 0 0]+[0 0 0;0 0 0;1 0 0]);
    b(5,:,:)=1/2^0.5*([0 1 0;0 0 0;0 0 0]+[0 0 0;1 0 0;0 0 0]);
    b(6,:,:)=1/3^0.5*eye(3); 
    
end