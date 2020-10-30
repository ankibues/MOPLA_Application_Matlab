% test Green Tensor with loop vs Vectorized Green tensor function.

   n     = 100;
 [Jd, ~, ~, ~] = FourIdentity();
   C_bar = 2*10*Jd;
   
   tic 
   %  obtain weights and nodes before the loop 
   gp                = 20;
   [p, w]            = Gauss(gp);
   ww                = w * w';
   [Alp1, Bet1, ww1] = Lebedev(86);
   [Alp2, Bet2, ww2] = Lebedev(974);
   [Alp3, Bet3, ww3] = Lebedev(5810);
   [Alp4, Bet4, ww4] = GaussGGLQ(80);
   [Alp5, Bet5, ww5] = GaussGGLQ(200);
   [Alp6, Bet6, ww6] = GaussGGLQ(210); 

%  generate a population of uniformly distributed RDEs with
%  a1:a2:1,a2 from a1 to 1.
   [A, ang] = RandAANG(10,n);
   %A = xlsread('a.xlsx','A2:B4');
   %ang= xlsread('a.xlsx','A7:B9');
   
   T = zeros(3,3,3,3,n);
   sss = zeros(3,3,3,3,n);
   ss = zeros(3,3,3,3,n);
   c_bar = zeros(3,3,3,3,n); 
   Carray = zeros(n,81);
   q   = zeros(3,3,n);
   m   = zeros(3,3,3,3,n);
   e0  = zeros(3,3,n);
   for i= 1:n
        
       q(:,:,i)= Q(ang(:,i));
      
   end
  
   for r=1:n
       c_bar(:,:,:,:,r) = Transform(C_bar,q(:,:,r));
       Carray(r,:)  = C2OneDarray(c_bar(:,:,:,:,r));
        
       T(:,:,:,:,r)     = TGreen(A(:,r),  Carray(r,:), Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
      

   end
   
    toc
    tic
  CC_bar= Transform2(C_bar,q,n);  
  T1 = SnpIn_all(A,CC_bar);
    toc
    
 
   
                           