function [rho,rho_i,X_i] = masstomole(T,P,Y_i,M_i)  

R=8.310;                        % universal gas constant J / mol. K

%molar fraction from mass fraction
M_avg= sum(Y_i./M_i);           %inverse of mixture molecular weight
X_i=(Y_i./M_i)/M_avg  ;          %Mole Fraction of species

%mixture density from mole fraction
rho_i = (X_i.*M_i)*P/(R*T);
rho = sum(rho_i./M_i)/M_avg;



