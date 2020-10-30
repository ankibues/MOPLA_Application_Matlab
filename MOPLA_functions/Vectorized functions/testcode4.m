%testcode_4_ eshelby tensor calculation for 'N' number of ellipoids. 
% Analytical vs Numerical results
% Old quadrature method vs New vectorized method
N=1;
a1=100;
a = [100;5;1];
a = repmat(a,1,N);
[~, ang] = RandAANG(a1,N);
S= zeros(3,3,3,3,N);
S_anal= zeros(3,3,3,3,N);
PI=zeros(3,3,3,3,N);

[Jd, Js, Ja, Jm] = FourIdentity();
Cm= 2*Jd;
      gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 

 %a = xlsread('a.xlsx','A2:B4');
%ang= xlsread('a.xlsx','A6:B8');

T = zeros(3,3,3,3,N);
   sss = zeros(3,3,3,3,N);
   ss = zeros(3,3,3,3,N);
   c_bar = zeros(3,3,3,3,N); 
   Carray = zeros(N,81);
   q   = zeros(3,3,N);
   m   = zeros(3,3,3,3,N);
   e0  = zeros(3,3,N);
   for i= 1:N
        
       q(:,:,i)= Q(ang(:,i));
      
   end
tic



    

%%% eshelby Tensor calculation (Qu 2016)
   for r=1:N
       c_bar(:,:,:,:,r) = Transform(Cm,q(:,:,r));
       Carray(r,:)  = C2OneDarray(c_bar(:,:,:,:,r));
        
       T(:,:,:,:,r)     = TGreen(a(:,r),  Carray(r,:), Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
      
        s1               = Contract(Jd, Contract(T(:,:,:,:,r),c_bar(:,:,:,:,r)));
       % s                = Transform(s1,q(:,:,r)');                     
        S(:,:,:,:,r)= s1;
   end



toc
%tic
% Eshelby tensor calculation ( Analytical results from Jiang 2016)

 for r=1:N
    [S_el,PI_el] = SnP(a(:,r));
% Eshelby Tensors (S,PI,p) for Interior points ( p, here is the green
% tensor for pressure)
% p
    p   = zeros(3,3);
    for j=1:3
        p(j,j) = -1/3* (S_el(j,j,1,1)+ S_el(j,j,2,2)+ S_el(j,j,3,3)); 
    end
%  S  
    S_anal(:,:,:,:,r)       = S_el;
    for k=1:3
        for l=1:3
            S_anal(k,k,l,l,r) = p(k,k)+ S_el(k,k,l,l);
        end
    end
PI(:,:,:,:,r)   = PI_el;
end
%toc
%}

%tic
% using new vectorized functions
%Cmm= Transform2(Cm,q,N);

%[SS,PPI] = SnpIn_all(a,Cmm,Jd,Ja);
%toc

    
    