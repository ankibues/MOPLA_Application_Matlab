function editdeforminput(steps,L)
%----------------This code edits VPSC7.in file for VPSC input----------------------------
        % Read txt into cell AAA
         fid = fopen('deform.in','r');
         j = 1;
         tline = fgetl(fid);
         AAA{j} = tline;
         while ischar(tline)
               j = j+1;
               tline = fgetl(fid);
               AAA{j} = tline;
         end
         fclose(fid);
             
         firstlin= sprintf('%d  7   0.05    298.         nsteps  ictrl  eqincr  temp',steps);
         seventhlin= sprintf('%4.2f      %4.2f      %4.2f		udot   |    vel.grad',L(1,1),L(1,2),L(1,3));
         eigthlin= sprintf('%4.2f      %4.2f      %4.2f		       |',L(2,1),L(2,2),L(2,3));
         ninthlin= sprintf('%4.2f      %4.2f      %4.2f		       |',L(3,1),L(3,2),L(3,3));	      	      

         % Change cell AAA
         AAA{1}= firstlin;
         AAA{7}=seventhlin;
         AAA{8}= eigthlin;
         AAA{9}=ninthlin;
         % Write cell AAA into txt
         fid = fopen('deform.in', 'w');
         for jj = 1:numel(AAA)
             if AAA{jj+1} == -1
                fprintf(fid,'%s', AAA{jj});
                break
             else
                fprintf(fid,'%s\r\n', AAA{jj});
             end
         end
         fclose(fid);
end
%}
