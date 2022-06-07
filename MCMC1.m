% Step 1
% This code computes the random walks via Markov Chain MonteCarlo method
% for obtaining 8000 sets of parameters, using as initial value the best
% fit parameters computed by SEIAR_covid_run_MEX.m
% This same code can be used to compute the parameters for each state,
% replacing SEIAR_covid_solver_MEX with
% SEIAR_covid_solver_(name of the state).

for ini = 1:length(D_exp)   %We compute the initial date at which D_exp and R_exp are positive
    if ( D_exp(ini)>0 )&&( R_exp(ini)>0 )
        break
    end
end

fin=length(C_exp);
Y=[C_exp(1,ini:end),R_exp(1,ini:end),D_exp(1,ini:end),A_exp(1,ini:end)];
% The initial values of parameters (x must be obtained previously)
beta0=x(1);
beta1=x(2);
tau_beta=x(3);
delta0=x(4);
delta1=x(5);
tau_delta=x(6);
w=x(7);
p=x(8); 
gamma0=x(9);
gamma1=x(10);
tau_gamma=x(11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt=linspace(0,113,114)';
xnew1 = [beta0,beta1,tau_beta,delta0,delta1,tau_delta,w,p,gamma0,gamma1,tau_gamma];
[T1,I1,RI1,D1,RA1,E1,A1]=SEIAR_covid_solver_MEX(xnew1,tt,S0,I0,RI0,RA0,E0,A0,N);  %Replace MEX with the name of the state
C1 = I1+RI1+D1;
B1=C1(ini:fin,1); B2=RI1(ini:fin,1); B3=D1(ini:fin,1); B4=A1(ini:fin,1);
n=[B1',B2',B3',B4'];
NewLik=sum(log(n).*Y-n-gammaln(Y+1));
Olik=NewLik;    % Calculate the logarithmic likelihood function for the initial set of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mrange=0.4.*[0.07, 0.005, 2, 0.01, 0.005, 2, 0.005, 0.01 0.1 0.1 2]; % The maximum range of random walks
NN=8000;mm=6000;  % MCMC cycle index and burn-in periods
MHPar=zeros(NN+1,11);
MHPar(1,:)=[beta0,beta1,tau_beta,delta0,delta1,tau_delta,w,p,gamma0,gamma1,tau_gamma];
% MCMC algorithm and Random walks to update parameters
for i=2:NN+1
 NPar=[beta0,beta1,tau_beta,delta0,delta1,tau_delta,w,p,gamma0,gamma1,tau_gamma] ...
     + mrange.*(2*rand(1,11)-1);  % Generate candidate parameters randomly
    Nbeta0=NPar(1);Nbeta1=NPar(2);Ntau_beta=NPar(3);Ndelta0=NPar(4);
    Ndelta1=NPar(5);Ntau_delta=NPar(6);Nw=NPar(7);Np=NPar(8);
    Ngamma0=NPar(9);Ngamma1=NPar(10);Ntau_gamma=NPar(11);
    if(Nbeta0>0&&Nbeta1>0&&Ntau_beta>0&&Ndelta0>0&&Ndelta1>0&& ...
            Ntau_delta>0&&Nw>0&&Np>0&&Ngamma0>0&&Ngamma1>0&&Ntau_gamma>0)  % Accept or reject
        xnew1 = [Nbeta0,Nbeta1,Ntau_beta,Ndelta0,Ndelta1,Ntau_delta,Nw,Np,Ngamma0,Ngamma1,Ntau_gamma];
        [T1,I1,RI1,D1,RA1,E1,A1]=SEIAR_covid_solver_MEX(xnew1,tt,S0,I0,RI0,RA0,E0,A0,N); %Replace MEX with the name of the state
        C1 = I1+RI1+D1;
        B1=C1(ini:fin,1); B2=RI1(ini:fin,1); B3=D1(ini:fin,1); B4=A1(ini:fin,1);
        n2=[B1',B2',B3',B4'];
        Nlik=sum(log(n2).*Y-n2-gammaln(Y+1)); % Calculate log. likelihood function for the new set of parameters
        alfa=min(1,exp(Nlik-Olik));
        if rand<alfa
            beta0=Nbeta0; beta1=Nbeta1; tau_beta=Ntau_beta; delta0=Ndelta0;
            delta1=Ndelta1; tau_delta=Ntau_delta; w=Nw; p=Np;
            gamma0=Ngamma0; gamma1=Ngamma1; tau_gamma=Ntau_gamma;
            Olik=Nlik;
        end
    end
    MHPar(i,:)=[beta0,beta1,tau_beta,delta0,delta1,tau_delta,w,p,gamma0,gamma1,tau_gamma];
end

Yu=MHPar(mm+1:NN,:);
save('Yu.mat','Yu')

% Calculate the estimation values of parameters
beta0=mean(MHPar(mm:end,1))
sigma_beta0=std(MHPar(mm:end,1));
beta0_ci=[beta0 - 1.96*sigma_beta0, beta0 + 1.96*sigma_beta0]

beta1=mean(MHPar(mm:end,2))
sigma_beta1=std(MHPar(mm:end,2));
beta1_ci=[beta1 - 1.96*sigma_beta1, beta1 + 1.96*sigma_beta1]

tau_beta=mean(MHPar(mm:end,3))
sigma_tau_beta=std(MHPar(mm:end,3));
tau_beta_ci=[tau_beta - 1.96*sigma_tau_beta, tau_beta + 1.96*sigma_tau_beta]

delta0=mean(MHPar(mm:end,4))
sigma_delta0=std(MHPar(mm:end,4));
delta0_ci=[delta0 - 1.96*sigma_delta0, delta0 + 1.96*sigma_delta0]

delta1=mean(MHPar(mm:end,5))
sigma_delta1=std(MHPar(mm:end,5));
delta1_ci=[delta1 - 1.96*sigma_delta1, delta1 + 1.96*sigma_delta1]

tau_delta=mean(MHPar(mm:end,6))
sigma_tau_delta=std(MHPar(mm:end,6));
tau_delta_ci=[tau_delta - 1.96*sigma_tau_delta, tau_delta + 1.96*sigma_tau_delta]

w=mean(MHPar(mm:end,7))
sigma_w=std(MHPar(mm:end,7));
w_ci=[w - 1.96*sigma_w, w + 1.96*sigma_w]

p=mean(MHPar(mm:end,8))
sigma_p=std(MHPar(mm:end,8));
p_ci=[p - 1.96*sigma_p, p + 1.96*sigma_p]

gamma0=mean(MHPar(mm:end,9))
sigma_gamma0=std(MHPar(mm:end,9));
gamma0_ci=[gamma0 - 1.96*sigma_gamma0, gamma0 + 1.96*sigma_gamma0]

gamma1=mean(MHPar(mm:end,10))
sigma_gamma1=std(MHPar(mm:end,10));
gamma1_ci=[gamma1 - 1.96*sigma_gamma1, gamma1 + 1.96*sigma_gamma1]

tau_gamma=mean(MHPar(mm:end,11))
sigma_tau_gamma=std(MHPar(mm:end,11));
tau_gamma_ci=[tau_gamma - 1.96*sigma_tau_gamma, tau_gamma + 1.96*sigma_tau_gamma]

M=[mean(MHPar(mm:end,1));mean(MHPar(mm:end,2));mean(MHPar(mm:end,3));
    mean(MHPar(mm:end,4));mean(MHPar(mm:end,5));mean(MHPar(mm:end,6));
    mean(MHPar(mm:end,7));mean(MHPar(mm:end,8));mean(MHPar(mm:end,9));
    mean(MHPar(mm:end,10));mean(MHPar(mm:end,11))]';

FangchaM=[std(MHPar(mm:end,1));std(MHPar(mm:end,2));std(MHPar(mm:end,3));
    std(MHPar(mm:end,4));std(MHPar(mm:end,5));std(MHPar(mm:end,6));
    std(MHPar(mm:end,7));std(MHPar(mm:end,8));std(MHPar(mm:end,9));
    std(MHPar(mm:end,10));std(MHPar(mm:end,11))]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('par.mat','M')
save('qujian.mat','FangchaM')