function editvpsc7SSinput(fname)
%----------------This code edits slip system input for VPSC7.in file ----------------------------
        % Read txt into cell AAA
         fid = fopen('vpsc7.in','r');
         j = 1;
         tline = fgetl(fid);
         AAA{j} = tline;
         while ischar(tline)
               j = j+1;
               tline = fgetl(fid);
               AAA{j} = tline;
         end
         fclose(fid);

         % Change cell AAA
         AAA{11} = fname;
         % Write cell AAA into txt
         fid = fopen('vpsc7.in', 'w');
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