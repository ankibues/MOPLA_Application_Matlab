function VPSC_run_script(Wk,SI,SS)
% Matlab function to automate VPSC modelling

cd 'D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\';
%WK= [0 .5 .75 .90 1];
tincr=.05;
%SI=[2 4 6 8];

%SS= ['quartz_low_A.sx'; 'quartz_low_B.sx';'quartz_low_C.sx';'quartz_low_D.sx';'quartz_low_E.sx';'quartz_low_F.sx'];


         gamma=1;
         epsilon= .5*(((1/(Wk^2))-1)^(.5));
         if Wk==0
            gamma=0;
            epsilon=1;
            L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
         else
            L= [epsilon, gamma ,0; 0 , -epsilon, 0; 0,0,0];
         end
         
         
         for i=1:1000
              [str,~,~]=Calc_Gamma(L,tincr,i);
              if str>=SI
                 steps=i;
                 break
              end
         end
         
         editdeforminput(steps,L)
         editvpsc7SSinput(SS)
         command= 'D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\vpsc7.win64.exe';
         system(command);
         fname = sprintf('test_3to1_SS=%s_Wk=%2.1f_SI=%2.1f',SS,Wk,SI);
         if Wk==0
            fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\singlescale\Wk=0\' fname];
         elseif Wk==.5
            fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\singlescale\Wk=.50\' fname];
         elseif Wk==.75
            fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\singlescale\Wk=.75\' fname];
         elseif Wk==.90
            fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\singlescale\Wk=.90\' fname];
         elseif Wk==1
            fname1= ['D:\Ankit\Dropbox\VPSC_Modeling2\Example\pureshear\Results\singlescale\Wk=1\' fname];
         end
         
         texfile= 'TEX_PH1.OUT';
         CS = crystalSymmetry('321',[4.9, 4.9, 5.4]);
         data = loadOrientation_generic(texfile,'CS',CS, 'ColumnNames', {'Euler1' 'Euler2' 'Euler3' 'Weights'},'Columns',[1,2,3,4],'Bunge');
         odf = calcODF(data);
         figure;
         plotPDF(odf,Miller(0,0,1,CS),'antipodal');
         setMTEXpref('xAxisDirection','East')
         set(gcf,'Visible','off');
        
         [path,name,~]= fileparts(fname1);
         fnamefig= [path '\' name '.fig'];
         savefig(gcf,fnamefig,'compact')
         filenamenew= [path '\' name, '.jpg'];
         saveas(gcf,filenamenew)
         close(gcf)
         
end