function F_i=sri_calculation(T,sri,Pr_i)
a=sri.a;
b=sri.b;
c=sri.c;
d=sri.d;
e=sri.e;
X=1/(1+(log10(Pr_i))^2)^-1;
F_i=d*T^e*(a*exp(-b/T)+exp(-T/c))^X;
end

