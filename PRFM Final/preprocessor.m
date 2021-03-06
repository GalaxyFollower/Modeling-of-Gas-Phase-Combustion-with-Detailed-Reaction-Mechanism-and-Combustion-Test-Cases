function [Y_i,M_i,vr,vp,Arf_k,nf_k,Ef_k,a,temp,reaction,react_matrix,prod_matrix,rev,third_body,third_body_matrix,...
reaction_order_matrix,reverse_reaction_order_matrix,no_of_reactions,spcs] = preprocessor()

%% Comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to load reaction and species data                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            DEFINITIONS                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% i in parameters denote species and k in parameters denote reactions
% phi                   = equivalence ratio
% elements              = elements in chemical reaction mechanism file
% spcs                  = species in chemical reaction mechanism file
% reaction              = reactions in chemical reaction mechanism file
% atomic mass           = Atomic mass of the 'elements'
% third_body            = third body present or not (if 1 third body present
% third_body_matrix     = enhanced third-body efficiencies matrix for each reaction
% rev                   = reversible (if 1 reversible,if 0 irreversible
% unbalanced_react_matrix/unbalanced_prod_matrix = coefficients extracted
% from reaction mechanism file for reactants and products respectively
% global reactants      = stoichiometric coeffients for reactants of global reaction 
% no_of_reactions       = No of reactions
% no_of_species         = No of species in complete reaction
% reaction_order_matrix = Arbitary reaction order from chemical reaction
% mechanism file for forward reaction
% reverse_reaction_order_matrix= Arbitary reaction order from chemical
% reaction mechanism file for reverse reaction
% X                     = Species molar fraction
% X_i                   = ith specie molar fraction
% react_matrix          = stoichiometric coefficient reactant 
% vr                    = forward reaction order of reaction 
% vp                    = reverse reaction order of reaction
% prod_matrix           = stochiometric coefficient product
% Tt                    = Temperature at time t
% Arf_k                 = Frequency factor in Arrhenius equation of a reaction i
% Ef_k                  = Activation Energy
% nf_k                  = Temperature exponent
% R                     = Universal Gas constant
% k(f/b)                = Reaction rate constant( f= forward, b= backward)
% q_i                   = Rate of progress
% q_if                  = Rate of reaction progress in forward direction
% q_ir                  = Rate of reaction progess in reverse direction
% w_i                   = Species Production Rate( f= forward, b= backward)
% RoF(f/b)=             Rate of Formation of specise at each reaction step( f= forward, b= backward)
% 
% T_RoF   = Sum of Rate of Formation of specise all reaction step
% dX_i    = Change of species molar concentration
% dt                    = time step
% t                     = time period
% X_it                  = Species molar concenration over time period
% Hf_i                  = Heat of formation of specie i (kJ/mol)
% Y_i                   = Species mass fraction of species
% M_i                    = Species molecular weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% if INPUT is PHI phi_on=1 else 0 
phi_on=0;
if phi_on==1
%% Load reaction mechanism file
[elements,spcs,reaction,atomic_mass,unbalanced_react_matrix,...
unbalanced_prod_matrix,rev,third_body,third_body_matrix,...
reaction_order_matrix,reverse_reaction_order_matrix,global_reactants,...
fuel_matrix,oxidizer_matrix,oxidizer_divide]=rxndata;

%% Equivalence Ratio
phi=1.2;
X_i=equivalance_ratio(global_reactants,fuel_matrix,oxidizer_matrix,phi,oxidizer_divide);
X_i=X_i';
else
%% Load reaction mechanism file  
[elements,spcs,reaction,atomic_mass,unbalanced_react_matrix,...
unbalanced_prod_matrix,rev,third_body,third_body_matrix,...
reaction_order_matrix,reverse_reaction_order_matrix]=rxndata;
end

Arf_k=[reaction(:).Arf]';
nf_k=[reaction(:).nf]';
Ef_k=[reaction(:).Ef]'.*4.184;

%% Load thermodynamic file
[a,temp,element_composition]=thermodat(spcs,elements);

%% Balancing Chemical Reaction
[react_matrix,prod_matrix]= stoic_balance(unbalanced_react_matrix,unbalanced_prod_matrix,element_composition);
disp('Balancing Done')

%% Calculate Atomic Weight of Species
M_i=element_composition*atomic_mass/1000;
[no_of_reactions,no_of_species]=size(react_matrix);

%% Mass fraction of the species for phi<>1
%According to species arrangement in chem.inp file
%  Y_i=[2/62;32/62;0.0;28/62];  
% Y_i=[0.0283;0.22642;0;0.74528];
% Y_i=[0.029960;0;0.149396;0;0;0.820643];  % CH4 validation
 % X_i = [0;0;.21;0;0;0;0;0;.79]
% MWmix=sum(X_i.*M_i);
% Y_i=X_i.*M_i./MWmix

% Y_i = [0;0;2.3292e-001;0;0;0;0;0;7.6708e-001]; %lateral flue mass fraction
% Y_i = [0;0;2.3292e-001;0;0;0;0;0;7.6708e-001]
% Y_i=[0.0028;0;0.00034;0;0;0;0;0.2557;0.7381];

X_i = [0;0;0;0;0;0;0;0;1]
% X_i = [0.0803;0.0278;0.0072;0.0053;0.0242;6.1894e-6;1.2691e-6;0.2682;0.5870]; %dillution
% X_i =[0.125843184;0.026189357;0.095630741;0.004662982;0.002995971;7.80146e-05; 9.91e-06;0.1370427;0.607547138]; %validation
MWmix=sum(X_i.*M_i);
Y_i=X_i.*M_i./MWmix;

% FINAL OUTPUT MUST BE Y_i

%% Convert to mass fraction for phi=1
if phi_on==1
 MWmix=sum(X_i.*M_i);
 Y_i=X_i.*M_i./MWmix;
end

%% Managing Stoichiometric Matrix for reactants and products
order_check=reaction_order_matrix==0;
vr=zeros(size(react_matrix));
for i=1: no_of_reactions
    vr(i,:)=react_matrix(i,:);
    for j=1:no_of_species
    if order_check(i,j)==0
    vr(i,j)=reaction_order_matrix(i,j);
    end
    end
end
order_check=reverse_reaction_order_matrix==0;
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

end



