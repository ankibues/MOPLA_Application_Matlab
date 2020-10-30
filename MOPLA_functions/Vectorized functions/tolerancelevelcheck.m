function tolerancelevelcheck(C_bar,q,r)  

%   S and H_arc for each RDP are calculated starting with
            %   the M_bar and E0 defined at the homogeneous macrosclae
            %   strain-rate state.
            
%  describe C in the RDE's coordinate system 
            c_bar = Transform2(C_bar,q);

            
            
%  compute the 4th-order Green tensor(T)% this has to be changed................

            T(:,:,:,:,r)     = TGreen(a(:,r),  Carray(r,:), Alp1, Bet1, ww1, Alp2, Bet2, ww2, Alp3, Bet3, ww3,...
                               Alp4, Bet4, ww4, Alp5, Bet5, ww5, Alp6, Bet6, ww6, p, ww); 
                           
                           
%  the 4th-order Interaction tensor(H_arc) 


            s1               = Contract(Jd, Contract(T(:,:,:,:,r),c_bar(:,:,:,:,r)));
            s                = Transform(s1,q(:,:,r)');                           
            H_arc(:,:,:,:,r) = Contract(FourTensorInv(FourTensorInv(s)-Jd),M_bar);  %Eq. 12b Jiang 2014
%  record the initial partitioned stress of RDE at the current homogeneous macroscopic state             
            stress           = sigma(:,:,r);
           
            for k=1:100
                %   Inner loop: B and beta(b) for each RDP are calculated starting with
                %   the M_bar and E0 defined at the homogeneous macrosclae
                %   strain-rate state uitil is appraches the balance.

%  the 4th-order stress-partitioning tensor(B)                
                t0              = Ne(1,r)*m(:,:,:,:,r);                            
                t1              = FourTensorInv((t0+H_arc(:,:,:,:,r)));                   
                B               = Contract(t1,(M_bar+H_arc(:,:,:,:,r)));           
                t2              = E0-e0(:,:,r);
%  the second order stress-partitioning tensor(beta)
                beta            = Multiply(t1,t2);   
%  calculate new partitioned stress field inside RDE
                sigma(:,:,r)    = Multiply(B,SIGMA)+beta;                          
%  compare the current partitioned stress and the previous one                 
                alpha           = abs(Inva(sigma(:,:,r)-stress)/Inva(stress));
%  replace the previous partitioned stress with the current one
                stress          = sigma(:,:,r);
%  calculate new partitioned strain rate of RDE               
                e(:,:,r)        = Multiply(m(:,:,:,:,r),sigma(:,:,r));                    
%  calculate new pre-strain rate term of RDE                
                e0(:,:,r)       = (1-Ne(1,r))*e(:,:,r);
%  calculate new strain rate invariant
                e1              = Inva(e(:,:,r));
%  calculate new viscosity of RDE at the new strain rate state               
                t3              = 1/Ne(1,r)-1;
                v1             =(e1/REF(1,r))^t3 *eta(1,r);     
%  calculate new compliance tensor of RDE
                m(:,:,:,:,r)    = 1/(2*v1)*Jd;                                     
               
%  The inner loop for an RDE terminates when the current stress coincides
%  with the previous one within a specific tolerance
               if alpha<0.01 
                   break
               end   
               
            end
%  update the strain-rate invariant for an RDE           
            eI(1,r)  = e1;  
%  update the effective viscosity for an RDE
            vis(1,r) = v1; 
            
            MB = MB+Contract(t0,B);                    %%<M(tan):B>*n%%
            Mb = Mb+Multiply(t0,beta)+e0(:,:,r);       %%<M(tan):beta+e0>*n%%
            BB = BB+B;                                 %%<B>*n$$
            bb = bb+beta;                              %%<beta>*n%%
        end