function k= Hm_vec(x,C,psi,c)
% This function calculates the Green function Hi for pressure.
% Note that since this is a vectorized version, Hi is a 3*N1 matrix rather
% than 3*1 matrix, as in non-vectorized case.

% Input: 

% x is a 3*N1 matrix,here N1 is corresponding to the theta/phi Gaussian
% grid nodes
% C is stiffness tensor , 3*3*3*3 matrix

% output is the vectorized version of Green's function for pressure 3*N1.

[~,N1]= size(x);
k= ones(3,N1);

psi=psi';
[N2,~]= size(c); % N2 corresponding to psi

modx= sqrt(x(1,:).^2 + x(2,:).^2 + x(3,:).^2);    %1*N1
%modx= sqrt((x(1)).^2 + (x(2)).^2 +(x(3)).^2);

[y,z,zT]= ZZ_vec(x,psi);  % z is 3*N1*N2....zT is 3*3*N1*N2....y is 3*N1 
B= AKm_vec(C,z);   % 4*4*N1*N2 
PHI= MultimixMMvec(C,zT); %3*3*N1*N2
c=repmat(c',N1,1);
s=zeros(N1,N2);

% this part is for calculation of Fi --------------------------------------
    for i=1:3
        
     
     
        for m=1:3
            for nn =1:3
               s = s + reshape(B(i,nn,:,:),N1,N2).*reshape(PHI(m,nn,:,:),N1,N2).*reshape(B(4,m,:,:),N1,N2) ; 
            end
               s= s + squeeze(reshape(B(4,4,:,:).*B(i,m,:,:)+ B(4,i,:,:).*B(4,m,:,:),N1,N2).*repmat(y(m,:)',1,N2));
        end
   
        
      kk= -1/4/(pi^2)./modx.^2;
      k(i,:)= kk.*(sum(c.*s,2)');   % 1*N1
      s=zeros(N1,N2);
    end
 %-------------------------------------------------------------------------  
end
