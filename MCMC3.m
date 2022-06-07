%Step 3
% Plots the simulation curves of the model and the confidence intervals
data_long=dateshift(date(1),'start','day',0:length(t)-1);
load('par.mat','M')
load('qujian.mat','FangchaM')
load('dataestp_C.mat','dataestp_C')
load('dataestp_I.mat','dataestp_I')
load('dataestp_RI.mat','dataestp_RI')
load('dataestp_D.mat','dataestp_D')
load('dataestp_A.mat','dataestp_A')
load('dataestp_R0.mat','dataestp_R0')
C1=1.96*std(dataestp_C,0,2);
I1=1.96*std(dataestp_I,0,2);
RI1=1.96*std(dataestp_RI,0,2);
D1=1.96*std(dataestp_D,0,2);
A1=1.96*std(dataestp_A,0,2);
R01=1.96*std(dataestp_R0,0,2);
par1=M;par2=FangchaM;
beta0=par1(1);sigma_beta0=par2(1);
beta1=par1(2);sigma_beta1=par2(2);
tau_beta=par1(3);sigma_tau_beta=par2(3);
delta0=par1(4);sigma_delta0=par2(4);
delta1=par1(5);sigma_delta1=par2(5);
tau_delta=par1(6);sigma_tau_delta=par2(6);
p=par1(7);sigma_p=par2(7);
w=par1(8);sigma_w=par2(8);
gamma0=par1(9);sigma_gamma0=par2(9);
gamma1=par1(10);sigma_gamma1=par2(10);
tau_gamma=par1(11);sigma_tau_gamma=par2(11);

figure(14)
plot(data_long,estp_R0','-k',data_long,estp_R0'+R01,':r',data_long,estp_R0'-R01,'-.g','LineWidth',2)
title('Daily reproduction number')
xlabel('date')
ylabel('R_d(t)')
xlim([data_long(1) data_long(end)])
ylim([0 8])
legend('R0','95% CI upper estimate','95% CI lower estimate');
set(gca,'FontSize',13)
hold on
plot(data_long,1+t-t,'--','LineWidth',2)

figure(15)
plot(data_long,estp_C,'-k',data_long,estp_C+C1,':r',data_long,estp_C-C1,'-.g','LineWidth',3)
hold on
plot(date,C_exp,'*')
title('Cumulative number of symptomatic cases')
legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
xlabel('Date')
ylabel('People')
set(gca,'FontSize',13)
grid on

figure(16)
plot(data_long,estp_D,'-k',data_long,estp_D+D1,':r',data_long,estp_D-D1,'-.g','LineWidth',3)
hold on
plot(date,D_exp,'*')
% ylim([0 1500])
title('Death toll (D(t))')
legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
xlabel('Date')
ylabel('People')
set(gca,'FontSize',13)
grid on

figure(17)
plot(data_long,estp_I,'-k',data_long,estp_I+I1,':r',data_long,estp_I-I1,'-.g','LineWidth',3)
hold on
plot(date,I_exp,'*')
title('Infected (I(t))')
legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
xlabel('Date')
ylabel('People')
set(gca,'FontSize',13)
grid on

figure(18)
plot(data_long,estp_RI,'-k',data_long,estp_RI+RI1,':r',data_long,estp_RI-RI1,'-.g','LineWidth',3)
hold on
plot(date,R_exp,'*')
title('Recovered (R_I(t))')
legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate','Reported data'},'Location','northwest')
xlabel('Date')
ylabel('People')
set(gca,'FontSize',13)
grid on

figure(19)
plot(data_long,estp_A,'-k',data_long,estp_A+A1,':r',data_long,estp_A-A1,'-.g','LineWidth',3)
title('Asymptomatic (A(t))')
legend({'Best fit solutions','95% CI upper estimate','95% CI lower estimate'},'Location','northwest')
xlabel('date')
ylabel('People')
set(gca,'FontSize',13)
grid on