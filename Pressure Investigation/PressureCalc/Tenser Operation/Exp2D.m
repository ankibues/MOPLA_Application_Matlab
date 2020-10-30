function z= Exp2D(A)
    % Input:A is a 2*2 matrix.

    omega=norm(A)/2^0.5;
    if omega==0
        z=eye(2);
    else
       omega1=A/omega;
       z=eye(2)+omega1*sin(omega)+(1-cos(omega))*omega1^2;
    end
end