% simulation of single marker FS for different initial parameters 
 k=pwd;
 K= extractBefore(k,"Flanking Structure Investigation");
 addpath(genpath(K));

Wk= [0 .50 .75];
r=[0.01 .1 .2 .5 2 5 10 100];

a1=[1000;.5;.1]; 
%a2=[1000;.8;.1];
%a3=[1000;.3;.1];
%a4= [1000;1;.1];
%a5= [1000;2;.1];
%A=[a1];%,a2,a3,a4,a5];
ANG= [pi/4, pi/2, 3*pi/4, pi/3,pi/6,2*pi/3];
n=1; %241 224 done
for i=1:length(r)
    for j=1:length(Wk)
        for k=1:length(ANG)
            %if n>224
                %Single_marker_FS(Wk(j),r(i),a1,ANG(k))
                disp('Simulation number=');
                disp(n);
            %end
            n=n+1;
        end
    end
end

