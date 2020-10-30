function z= Exp(A)
   % Input:A is a 3*3 matrix.
   % THis function is for Rodrigues rotation approximation.
    omega=norm(A,'fro')/2^0.5;
    if omega==0
        z=eye(3);
    else
       omega1=A/omega;
       z=eye(3)+omega1*sin(omega)+(1-cos(omega))*omega1^2;
    end
end