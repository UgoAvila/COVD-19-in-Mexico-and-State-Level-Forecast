%Step 2
% This code computes the solutions of the model for each of the 2000 sets
% of parameters (after the burn-in period) obtained in Step 1

load('Yu.mat','Yu')
% t=linspace(0,360,361);
dataestp_C=zeros(length(t),NN-mm);
dataestp_I=dataestp_C; dataestp_RI=dataestp_C;
dataestp_D=dataestp_C; dataestp_A=dataestp_C;
dataestp_R0=dataestp_C;
for pp=1:NN-mm
    PARAS=Yu(pp,:);
%     beta0=PARAS(1);
%     beta1=PARAS(2);
%     tau_beta=PARAS(3);
%     delta0=PARAS(4);
%     delta1=PARAS(5);
%     tau_delta=PARAS(6);
%     p=PARAS(7);
%     w=PARAS(8);
%     gamma0=PARAS(9);
%     gamma1=PARAS(10);
%     tau_gamma=PARAS(11);
    [t,I,RI,D,RA,E,A]=SEIAR_covid_solver_MEX(PARAS,t,S0,I0,RI0,RA0,E0,A0,N);
    %Replace MEX with the name of the state
    C = I + RI + D;
    delta_t=zeros(1,length(t));
    beta_t=zeros(1,length(t));
    gamma_t=zeros(1,length(t));
       for i=1:length(t)
        delta_t(i)=PARAS(4)*exp( -t(i)/PARAS(6) ) + PARAS(5);
        beta_t(i)=PARAS(1)*exp( -t(i)/x(3) ) + PARAS(2);
        gamma_t(i) = ( PARAS(10)./(1+exp(-t(i)+PARAS(11))) + PARAS(9) );
       end
    Rdt = beta_t.*p./(delta_t + gamma_t) + beta_t.*(1-p)./gamma_t;
    dataestp_C(:,pp)=C; dataestp_I(:,pp)=I; dataestp_RI(:,pp)=RI;
    dataestp_D(:,pp)=D; dataestp_A(:,pp)=A;
    dataestp_R0(:,pp)=Rdt;
end
estp_C=C;estp_I=I; estp_RI=RI; estp_D=D; estp_A=A; estp_R0=Rdt;
save('estp_C.mat','C')
save('estp_I.mat','I')
save('estp_RI.mat','RI')
save('estp_D.mat','D')
save('estp_A.mat','A')
save('estp_R0.mat','Rdt')
save('dataestp_C.mat','dataestp_C')
save('dataestp_I.mat','dataestp_I')
save('dataestp_RI.mat','dataestp_RI')
save('dataestp_D.mat','dataestp_D')
save('dataestp_A.mat','dataestp_A')
save('dataestp_R0.mat','dataestp_R0')