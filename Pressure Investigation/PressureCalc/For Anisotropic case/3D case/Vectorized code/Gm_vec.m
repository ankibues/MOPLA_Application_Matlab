function k= Gm_vec(x,C,psi,c)
% x is a 3*N1 matrix,here N1 is corresponding to the theta/phi gaussian
% points; C is 3*3*3*3
% output is 3*3*3*3*N1.       

[~,N1]= size(x);
k= ones(3,3,3,N1);
% x is the coordinates, (3*1)
% C is stiffness tensor (3*3*3*3)

%psi= (pi/2)*psi + pi/2;
  % to change the limits of integration
psi=psi';
[N2,~]= size(c); % N2 corresponding to psi

modx= sqrt(x(1,:).^2 + x(2,:).^2 + x(3,:).^2);    %1*N1
%modx= sqrt((x(1)).^2 + (x(2)).^2 +(x(3)).^2);

[y,z,zT]= ZZ_vec(x,psi);  % z is 3*N1*N2....zT is 3*3*N1*N2....y is 3*N1 
B= AKm_vec(C,z);   % 4*4*N1*N2 
PHI= MultimixMMvec(C,zT); %3*3*N1*N2
c=repmat(c',N1,1);
s=zeros(N1,N2);

    for i=1:3
       for j=1:3
           for l=1:3
               
       % this part is for calculation of Mij  --------------------------- 
                           
        for m=1:3
            for nn =1:3
               s = s + reshape(B(i,m,:,:),N1,N2).*reshape(B(j,nn,:,:),N1,N2).*reshape(PHI(m,nn,:,:),N1,N2) ; 
            end
               s= s + squeeze(reshape(B(i,m,:,:).*B(4,j,:,:)+ B(j,m,:,:).*B(4,i,:,:),N1,N2).*repmat(y(m,:)',1,N2));
        end
        
         %---------------------------------------------------------------
         % here, s = Mij, now the integrand I becomes:
         % something equivalent of this I= (-y(l,:).*B(i,j,:)) + (z(l,:).*s);
       
      I=  squeeze(repmat(-y(l,:)',1,N2).*reshape(B(i,j,:,:),N1,N2)) + (squeeze(z(l,:,:)).*s);
      
      kk= 1/4/(pi^2)./modx.^2;
      k(i,j,l,:)= kk.*(sum(c.*I,2)');   % 1*N1
      s=zeros(N1,N2);
           end
       end
    end
end


