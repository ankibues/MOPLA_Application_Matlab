% testing and optimizing velocity field calculation 

load('test_data1.mat');
tic
vel_f1= Velocity_field_calc4(xxx(30),yyy(30),a,invS,Jd,d,u2,w,ll);
toc
tic
vel_f2= Velocity_field_calc3(xxx(30),yyy(30),a,invS,Jd,d,u2,w,ll);
toc
