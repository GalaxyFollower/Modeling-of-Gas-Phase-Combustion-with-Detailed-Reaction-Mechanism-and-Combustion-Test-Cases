function F_i = troe_calculation(T,troedata,Pr_i)
alp=troedata.alpha;
T3=troedata.T3;
T2=troedata.T2;
T1=troedata.T1;
F_cent1=(1-alp)*exp(-T/T3)+alp*exp(-T/T1);
if T2==0
    F_cent2=0;
else
    F_cent2=exp(-T2/T);
end
F_cent=F_cent1+F_cent2;

A_troe=log10(Pr_i)-0.67*log10(F_cent)-0.4;
B_troe=0.806-1.1762*log10(F_cent)-0.14*log10(Pr_i);

F_i= exp(((1+(A_troe/B_troe)^2)^-1)* log(F_cent));

end