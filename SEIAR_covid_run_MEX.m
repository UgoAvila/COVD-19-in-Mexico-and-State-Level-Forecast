%% Data From The Ministry of Health of Mexico with the number of infections and death by state (We retrieve Mexican data)
     url_Iglobe=('https://raw.githubusercontent.com/UgoAvila/COVD-19-in-Mexico-and-State-Level-Forecast/master/Confirmed_Infected_National.csv');
     urlwrite(url_Iglobe,'Iglobe.csv');
     I_globe=readtable('Confirmed_I.csv','ReadRowNames',true);
     date=datetime(I_globe{1,1:end},'Format','dd/MM/yy');

     url_Dglobe=('https://raw.githubusercontent.com/UgoAvila/COVD-19-in-Mexico-and-State-Level-Forecast/master/Confirmed_Death_National.csv');
     urlwrite(url_Dglobe,'Dglobe.csv');
     D_globe=readtable('Confirmed_D.csv','ReadRowNames',true);

     url_Dglobe=('https://raw.githubusercontent.com/UgoAvila/COVD-19-in-Mexico-and-State-Level-Forecast/master/Confirmed_Recuperated_National.csv');
     urlwrite(url_Dglobe,'Dglobe.csv');
     R_globe=readtable('Confirmed_R.csv','ReadRowNames',true);
 
    
 %%
     C_exp_mx=str2double(I_globe{{'Nacional'},2:end});       %Cumulative number of infected cases
     D_exp_mx=str2double(D_globe{{'Nacional'},2:end});  %Number of deaths
     R_exp_mx=str2double(R_globe{{'Nacional'},2:end});  %Number of recoveries
     date=date(2:end);               %Day 1 corresponds to 13/03/20
     C_exp=C_exp_mx(1:end);
     D_exp=D_exp_mx(1:end);
     R_exp=R_exp_mx(1:end);
     I_exp=C_exp-R_exp-D_exp;
     A_exp=I_exp*8;                   %Approximate number of asymptomatic cases
     
%%
t=linspace(0,140,141);
%x0 contains the best fit parameters from the Mexican outbreak
x0 = [0.4668 0.0100 27.8602 0.0 0.0829 11.2885 0.254 0.1  0.01 0.06 30];
lb = [0.04   0      0       0.0 0      0       0.05  0.08 0    0    0];
ub = [1      1      100     0.1 0.1    100     0.3   0.14 1    1    100];

N=127090000; % The population of Mexico
I0=I_exp(1);
RI0=R_exp(1);
RA0=RI0;
E0=I0;
A0=A_exp(1);

S0 = N-I0-RI0-RA0-E0-A0-D_exp(1);

%% Optimization section
% 1st Optimization with gradient based algorithm
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(@(x)SEIAR_covid_sse_MEX(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,A_exp,N),x0,[],[],[],[],lb,ub,[],options);

% Optimization with patternsearch
x0=x;
options = psoptimset;
options = psoptimset(options,'Display', 'off');
options = psoptimset(options,'PlotFcns', { @psplotbestf });
%options = psoptimset(options, 'MaxIter', 500);
[x,fval,exitflag,output] = patternsearch(@(x)SEIAR_covid_sse_MEX(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,A_exp,N),x0,[],[],[],[],lb,ub,[],options);

% 2nd Optimization with gradient based algorithm
x0=x;
options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @optimplotfval });
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)SEIAR_covid_sse_MEX(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,A_exp,N),x0,[],[],[],[],lb,ub,[],options);

%% Solving the model with the best fit parameters

t=linspace(0,230,231);      %Time span to calculate the predictions
[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_MEX(x,t,S0,I0,RI0,RA0,E0,A0,N);

%% Graphs

%Plot the kinetics of the outbreak
for i=1:length(t)
        delta_t(i)=x(4)*exp( -t(i)/x(6) ) + x(5);
        beta_t(i)=x(1)*exp( -t(i)/x(3) ) + x(2);
end

gamma = ( x(10)./(1+exp(-t+x(11))) + x(9) )';

data_long=dateshift(date(1),'start','day',0:t(end));

figure(12)
plot(data_long,I,':r','LineWidth',2)
hold on
plot(date,I_exp,'r*')
plot(data_long,D,'-k','LineWidth',2)
plot(date,D_exp,'ko')
ylabel('Infected and Dead')
hold on
yyaxis right
plot(data_long,RI,'-.g','LineWidth',2)
plot(date,R_exp,'gd')
plot(data_long,A,'--b','LineWidth',2)
% ylim([0 4e6])
ylabel('Recovered and Asymptomatic')
hold off
legend({'Infected','Infected data','Dead','Dead data','Recovered','Recovered data','Asymptomatic'},'Location','northwest')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figure(13)
plot(data_long,I,':r','LineWidth',3)
hold on
plot(date,I_exp,'r*')
plot(data_long,D,'-k','LineWidth',3)
plot(date,D_exp,'ko')
xlabel('Date')
ylabel('People')
grid on
plot(data_long,RI,'-.g','LineWidth',3)
plot(date,R_exp,'gd')
plot(data_long,A,'--b','LineWidth',3)
ylim([0 .8e6])
legend({'Infected','Infected data','Dead','Dead data','Recovered','Recovered data','Asymptomatic'},'Location','northwest')
% ax = gca;
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'k';

figure(3)
plot(t,beta_t','-b','LineWidth',4)
ylabel('Infection rate (1/day)')
xlabel('Days')
hold on
yyaxis right
plot(t,delta_t',':r','LineWidth',4)
hold on
plot(t,gamma','-.g','LineWidth',4)
ylabel('Death and Recovery rate (1/day)')
hold off
legend('Infection rate','Death rate', 'Recovery rate')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
xlim([0 230])
grid on

figure(4)
p=x(8);
Rdt = beta_t.*p./(delta_t + gamma) + beta_t.*(1-p)./gamma;
plot(data_long,1+t-t,'--','LineWidth',2)
ylim([0 16])
hold on
plot(data_long,Rdt,'-','LineWidth',2)
title('Daily reproduction number')
xlabel('date')
ylabel('R_d(t)')

C = I + RI + D;
figure(5)
plot(data_long,C,'-','LineWidth',2)
hold on
plot(date,C_exp,'*')
title('Cumulative number of symptomatic cases')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(6)
plot(data_long,D,'-','LineWidth',2)
hold on
plot(date,D_exp,'*')
title('Death toll')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(7)
plot(data_long,I,'-k','LineWidth',2)
hold on
plot(date,I_exp,'*')
title('Infected')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(8)
plot(data_long,I,'-','LineWidth',2)
hold on
plot(data_long,A,'--','LineWidth',2)
title('Infected')
legend({'Symptomatic (I(t))','Asymptomatic (A(t))'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(9)
plot(date,RI([1:length(C_exp)]),'-k')
hold on
plot(date,R_exp,'*')
title('Recovered')
legend({'Model solutions','Reported data'},'Location','northwest')
xlabel('days')
ylabel('People')
grid on

figure(10)
plot(date,A([1:length(C_exp)]),'-')
title('Asymptomatic cases')
xlabel('days')
ylabel('People')
grid on