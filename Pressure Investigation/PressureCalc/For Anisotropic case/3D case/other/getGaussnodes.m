function [m,n]= getGaussnodes(ii,jj,a,x,Cm)
% this code give us appropriate Gaussian nodes as a function of
% lambda for each component of Pressure auxiliary tensor
n=20;
m=30;
diff1=.5;
diff2=.5;
[Alpp,Bett,www] = GaussGGQ(m);
theta = reshape(Alpp,1,[]);
phi   = reshape(Bett,1,[]);
Wout = reshape(www,1,[]);
[psi,Win] = Gauss1(n);

kk0= GGQ3DPressOut([ii,jj], a,x,0,2*pi, 0,pi,theta, phi, Wout,Cm,psi,Win);
k0=kk0;
     while diff1>0.001
         while diff2>0.001
                n=n+5;
                [Alpp,Bett,www] = GaussGGQ(m);
                theta = reshape(Alpp,1,[]);
                phi   = reshape(Bett,1,[]);
                Wout = reshape(www,1,[]);
                [psi,Win] = Gauss1(n);
                k_next= GGQ3DPressOut([ii,jj], a,x,0,2*pi, 0,pi,theta, phi, Wout,Cm,psi,Win);
                diff2=k_next-k0;
                k0=k_next;
         end
         m=m+10;
         [Alpp,Bett,www] = GaussGGQ(m);
         theta = reshape(Alpp,1,[]);
         phi   = reshape(Bett,1,[]);
         Wout = reshape(www,1,[]);
         [psi,Win] = Gauss1(n);
         k_n= GGQ3DPressOut([ii,jj], a,x,0,2*pi, 0,pi,theta, phi, Wout,Cm,psi,Win);
         diff1=k_n-kk0;
         kk0=k_n;
     end
     
         
     end


%end
