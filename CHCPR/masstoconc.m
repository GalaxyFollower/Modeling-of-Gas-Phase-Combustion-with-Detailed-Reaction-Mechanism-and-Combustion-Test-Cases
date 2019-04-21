function [X_iconc,X_i] = masstoconc(T,P,Y_i,M_i)   %[X_i,X_avg,...] = molar_fraction(i,c,...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculated molar fraction and concentration used in reaction progression for given species
% mass fractions
% SPECIES: [H2,O2,OH,N2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P=101325;                      % pressure in pascal
R=8.31;                         % universal gas constant J / mol. K
%T=1000                         % temperature K
%Y_i = [0.03225;0.51610;0.0;0.4516];
%Y_i = [0.03809;0.30472;0;0.65719]; 
% initialize stoichoimetric mixture, inert specie->N2

%M_i = [2;32;17;28];         %Molar Weight of species

%% For molar fraction from mass fraction
M_avg= 1/sum(Y_i./M_i);     %Average molecular mass
X_i=(Y_i./M_i)*M_avg;    %Mole Fraction of species

%% For molar concentration from mole fraction
X_iconc=(P.*X_i)/(R*T)*10^(-6);
% X_iconc2=(P*(Y_i./M_i))/(R*sum(Y_i.*T./M_i))*10^-6
%needed in mol/cm^3 because Arhenius constant's unit is consistent with
%concentrations in terms of mol,cm,sec,K.
end