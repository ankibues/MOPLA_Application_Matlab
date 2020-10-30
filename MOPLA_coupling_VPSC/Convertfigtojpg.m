
% This is the code to convert matlab fig files in the results folder to jpg
% files for quick interpretation

srcfiles= dir('D:\Ankit\Dropbox\VPSC_Modeling\Example\pureshear\Results\IsotropicHEM\Wk=1\r=.5\*.fig');
 %cd 'F:\Dropbox\VPSC_Modeling\Example\pureshear\';
for i=1:numel(srcfiles)
    fig= openfig(srcfiles(i).name);
    [filepath,nam,ext]= fileparts(srcfiles(i).name);
    filenamenew= [filepath, nam, '.jpg'];
    saveas(fig,filenamenew)
end
