function [ val, df ] = augfailTP(T,P,V,f,ft,fp) 
%augfailTP Takes a log linear adjustment to the failure rate, f, with T, P
%   Uses a lever arm between the Maximum variable, the variable failure
%   rate and the actual failure rate to determine the adjusted failure rate
%   for that adjustment in terms of a ratio that is the failure rate is
%   multiplied by, df.

regf=(f*10^6);
regfp=(fp*10^6);
regft=(ft*10^6);

[valT, ind] = max([T-25,25-T]);
maxT=[600,125];
valT= min(valT,maxT(ind));

tauT=log10(regft/regf)/maxT(ind);

df(1)=10^(valT*tauT);

[valP, ind] = max([P-1,1-P]);
maxP=[10^9,1000];
%         valP= max(valP,maxP(ind));


valP = (P*V*(log(P)-(1-1/P)))/101.325; %Joules

valP = min(valP,maxP(ind));



tauP=log10(regfp/regf)/maxP(ind);
df(2)=10^(valP*tauP);

val=regf*df(1)*df(2);

val=val*10^-6;
end