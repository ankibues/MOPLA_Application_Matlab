WK= [0 .5  1];

SI=[2 6];

%SlipSys= ['quartz_low_A.sx';'quartz_low_B.sx';'quartz_low_C.sx';'quartz_low_D.sx';'quartz_low_E.sx';'quartz_low_F.sx';'quartz_low_G.sx';'quartz_low_H.sx';'quartz_low_I.sx'];


for kk1= 1:numel(WK)
    for kk2=1:numel(SI)
        %for kk3=1:9
            
            Wk=WK(kk1);
            StrI= SI(kk2);
            %SS= SlipSys(kk3,:);
            VPSC_run_script2(Wk,StrI,'quartz_low_F.sx')
        %end
    end
end


