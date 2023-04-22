f4m = importdata("tovErr_rk4m_fixed.log");
f5m = importdata("tovErr_rk5m_fixed.log");
f5f = importdata("tovErr_rk5f_fixed.log");
f5dp = importdata("tovErr_rk5dp_fixed.log");
f5ck = importdata("tovErr_rk5ck_fixed.log");

clr=[0.94118 0.50196 0.50196;0.56078 0.73725 0.56078;0.39216 0.58431 0.92941;0.93 0.57 0.13];

hold on;
p4m=plot(f4m(:,2),"k");
p5m=plot(f5m(:,2),"Color",clr(1,:));
p5f=plot(f5f(:,2),"Color",clr(2,:));
p5dp=plot(f5dp(:,2),"Color",clr(3,:));
p5ck=plot(f5ck(:,2),"Color",clr(4,:));

set(gca, "YScale", "log", "Box", "on");
legend([p4m, p5m, p5f, p5dp, p5ck], {'RK-4M','RK-5M','RK-5F','RK-5DP','RK-5CK'});legend('boxoff');