function sset=SEIAR_covid_sse_MEX(x,t,S0,I0,RI0,RA0,E0,A0,C_exp,R_exp,D_exp,A_exp,N)

[t,I,RI,D,RA,E,A]=SEIAR_covid_solver_MEX(x,t,S0,I0,RI0,RA0,E0,A0,N);

Cmod = I+RI+D;  %Cumulative number of infectives (calculated with the model)
sse_c=0;
sse_r=0;
sse_d=0;
sse_a=0;
for i=1:length(C_exp)
    sse_c=sse_c+(Cmod(i)-C_exp(i))^2;
    sse_d=sse_d+(D(i)-D_exp(i))^2;
    if i>1
        if R_exp(i)>R_exp(i-1)
            sse_r=sse_r+(RI(i)-R_exp(i))^2;
            sse_a=sse_a+(A(i)-A_exp(i))^2;
        end
    end
end

sset = 10*sse_c + sse_r + 10*sse_d + sse_a;

end