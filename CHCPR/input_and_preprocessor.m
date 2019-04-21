function [reaction,spcs,X_i,X_iconc,a,temp,react_matrix,prod_matrix,...
    rev,third_body,third_body_matrix,low,troe,sri,vr,vp,M_i,Y_i]...
    =input_and_preprocessor(P,R,Tt)
%% PREPROCESSOR AND INPUT PPROGRAM. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMS THE PROCESSING OF CHEMICAL REACTION MECHNISM FILE,
% THERMOCHEMICAL DATA FILE, MASS TABLE FILE.

%%%%%%%%%%%%%%   DEFINITIONS    %%%%%%%%%%%%%%%%%
%% i in parameters denote species and k in parameters denote reactions
%% MAIN PARAMETERS:
% phi                   = equivalence ratio
% elements              = elements in chemical reaction mechanism file
% spcs                  = species in chemical reaction mechanism file
% reaction              = reactions in chemical reaction mechanism file
% atomic mass           = Atomic mass of the 'elements'
% third_body            = third body present or not(if 1 third body present
% third_body_matrix     = enhanced third-body efficiencies matrix for each reaction
% rev                   = reversible (if 1 reversible,if 0 irreversible
% no_of_reactions       = No of reactions
% no_of_species         = No of species in complete reaction
% X_i                   = Species molar fraction
% react_matrix          = Stoichiometric coefficient reactant 
% vr                    = Forward reaction order 
% vp                    = Reverse reaction order
% prod_matrix           = Stochiometric coefficient product
% Tt                    = Temperature at time t
% low              = lower limit Arrhenius paramters presence check(0 or 1)
% sri              = SRI paramters presence check in each reaction 
% troe             = TROE paramters presence check in each reaction 
% Y_i                   = Mass Fraction of species
% M_i                   = Species molecular weight

%% TEMPORARY PARAMETERS
% unbalanced_react_matrix/unbalanced_prod_matrix = coefficients extracted
% from reaction mechanism file for reactants and products respectively
% global reactants      = stoichiometric coeffients for reactants of global reaction 
% reaction_order_matrix = Arbitary reaction order from chemical reaction
% mechanism file for forward reaction
% reverse_reaction_order_matrix= Arbitary reaction order from chemical
% reaction mechanism file for reverse reaction

% fuel_matrix            = fraction of species of fuel (in order of given
%speecies in chemical reaction mechanism file)
% oxidizer_matrix        = fraction of species of oxidizer
% oxidizer_divide        = fraction of O2 or main oxdizer in oxidizer mixture
% a                      = curve fitting coefficients for NASA polynomials
% temp                   = temperature range for validity of coefficients
% elemental_composition  = composition of elements in each species

% clear;
% clc;
%  P=101325;
%  R=8.314;
%  Tt=1500;
%% if INPUT is PHI phi_on=1 else 0 
%if phi_on=1. Supply global reaction in defined format(required) 
%No need to give initial mass fraction. It is calculated from
%global reaction in the chemical reaction mechanism file.
%%%%%%%%%%%%%
phi_on=1;  %%
%%%%%%%%%%%%%
if phi_on==1
%% Load data from chemical reaction mechanism file 'chem.inp'
[elements,spcs,reaction,atomic_mass,unbalanced_react_matrix,...
unbalanced_prod_matrix,rev,third_body,third_body_matrix,low,troe,sri,...
reaction_order_matrix,reverse_reaction_order_matrix,global_reactants,...
fuel_matrix,oxidizer_matrix,oxidizer_divide,fuel_divide]=rxndata(phi_on);

%% Calculate mole fraction from given Equivalence Ratio

%%%%%%%%%%%%%%%%%%%
phi=1.0;          %%%
%%%%%%%%%%%%%%%%%%%
X_i=equivalance_ratio(global_reactants,fuel_matrix,oxidizer_matrix,phi,oxidizer_divide,fuel_divide);
X_i=X_i';
else
%% Load reaction mechanism file  
[elements,spcs,reaction,atomic_mass,unbalanced_react_matrix,...
unbalanced_prod_matrix,rev,third_body,third_body_matrix,low,troe,sri,...
reaction_order_matrix,reverse_reaction_order_matrix]=rxndata(phi_on);

%% Mass fraction of the species
%According to species arrangement in chem.inp file
Y_i=[0.0283;0.22642;0;0.74528;0;0;0;0;0];
end
%% Load thermochemical data from 'thermo.dat' file
[a,temp,element_composition]=thermodat(spcs,elements);

%% Balancing Chemical Reaction and Obtaining Stoichiometric Coefficients for each reaction
[react_matrix,prod_matrix]= stoic_balance(unbalanced_react_matrix,unbalanced_prod_matrix,element_composition);
disp('Balancing Done')

%% Calculate Atomic Weight of Species
M_i=element_composition*atomic_mass;
[no_of_reactions,no_of_species]=size(react_matrix);
%% Managing Reaction Order for reactants and products

order_check=(reaction_order_matrix==0);
vr=zeros(size(react_matrix));
for i=1: no_of_reactions
    vr(i,:)=react_matrix(i,:);
    for j=1:no_of_species
    if order_check(i,j)==0
    vr(i,j)=reaction_order_matrix(i,j);
    end
    end
end
order_check=(reverse_reaction_order_matrix==0);
vp=zeros(size(prod_matrix));
for i=1: no_of_reactions
    vp(i,:)=prod_matrix(i,:);
    for j=1:no_of_species
    if order_check(i,j)==0
    vp(i,j)=reverse_reaction_order_matrix(i,j);
    end
    end
end
    
vr=vr';     %% matrix of forward reaction order
vp=vp';     %% matrix of reverse reaction order


%% Convert mass fraction to molar concentration or mole fraction to molar
%concentration as required

if phi_on==1
 X_iconc=moletoconc(P,X_i,R,Tt);
 Y_i=X_i.*M_i./sum(X_i.*M_i);
else
[X_iconc,X_i]=masstoconc(Tt,P,Y_i,M_i);
end

end




