
clc; clear;
%% This is the main program.
%% Comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to compute & plot the change in mass concentration in reactions %
% with respect to time.                                                   %
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
% Y_i                   = Species mass fraction of species
% M_i                   = Species molecular weight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial Conditions
dt=1e-8;            % each incremental time
t=1e-5;             %total reaction time
Tt=1500;            %temperature K
P=2.95*101325;           %input pressure in pascal
R=8.3144598;        %J/mol K
tspan=[0:dt:t];

%%Default input is equivalence ratio.To change default input open
%%"input_and_preprocessor.m"

%% Preprocessor and Input
[reaction,species,X_i,X_iconc,a,temp,react_matrix,prod_matrix,...
rev,third_body,third_body_matrix,low,troe,sri,vr,vp,M_i,Y_i]...
=input_and_preprocessor(P,R,Tt);

%% Thermodynamic Data and Reaction Rate Plots
preplots=0;
if preplots==1
Rcgs=1.998154; %data in cgs units to compare with CHEMKIN
preprocessor_plots(reaction,species,a,temp,react_matrix,prod_matrix,rev,Rcgs);
end

%% Initializing values  
M_i=M_i/1000;
M_avg=sum(X_i.*M_i);
[no_of_reactions,no_of_species]=size(react_matrix);

%% For mass flow rate calcualation only
    %Initial enthalpy and specific heat
    [cp_i,H_i,S_i]=nasa(a,temp,Tt,R);
    cp_total=sum(cp_i.*X_i./M_avg);
    mixture_enthalpy_initial=sum(Y_i.*(H_i)./M_i);
    
%% Start of CHCPR Calculation
y0 = zeros(1,no_of_species+1);
y0(1,1:no_of_species)=X_iconc;           %Species Molar Concentration
y0(1,no_of_species+1) = Tt;              %Temperature                

%% Main Solver
options = odeset('RelTol',1e-8,'AbsTol',1e-20);
[t,y] = ode15s(@(t,y) odefunction(t,y,reaction,react_matrix,prod_matrix,...
rev,third_body,third_body_matrix,vr,vp,a,temp,R,P,M_i),tspan,y0,options);

%% Post Processing     
RE = cell(numel(t),5);          % number of postprocessing variables = 3
for ii = 1:numel(t)    
[~,RE(ii,:)] = odefunction(t(ii),y(ii,:)',reaction,react_matrix,...
prod_matrix,rev,third_body,third_body_matrix,vr,vp,a,temp,R,P,M_i);
end
%%  
X_it = molfrctn2(y(:,1:no_of_species));
W_i = cell2mat(RE(:,1));
Q_v = cell2mat(RE(:,4));
density=cell2mat(RE(:,5));
%%
Tt_final=y(numel(t),no_of_species+1); %final temperature
Xi_final=X_it(numel(t),1:no_of_species)'; %final mole fraction
%% Ignition and Reaction Time and Temperature
% 
% Tig=Tt+0.05*(y(numel(t),no_of_species+1)-Tt); % Ignition Temperature
% sprintf('Ignition Temperature= %f',Tig)
% 
% Trt=Tt+0.95*(y(numel(t),no_of_species+1)-Tt); % Reaction Temperature
% for i=1:numel(t)
%     if Tig<y(i,no_of_species+1)
%         tig=t(i);
%         sprintf('Ignition time= %f',tig)
%         break;
%     end
% end
% for i=1:numel(t)
%     if Trt<y(i,no_of_species+1)
%         trt=t(i);
%         break;
%     end
% end
% 
%  sprintf('Reaction time= %f',trt-tig)

%% Plots  
%% Molar Fraction
figure(50)
plot1=plot(t,X_it,'LineWidth',2);
for i=1:no_of_species
set(plot1(i),'DisplayName',species(i).species);
end
xlabel('Reaction Time(s)');
ylabel('Species Molar Fraction');
title('Finite Rate Reaction Progression');
legend1 = legend('show');
set(legend1,'Location','best');
%% Temperature
figure(51)
plot(t,y(:,no_of_species+1),'k','LineWidth',2)
xlabel('Reaction Time(s)')
ylabel('Temperature(K)')
title('Temperature during Reaction progression')
%% Volumetric Heat Release
figure(52)
plot(t,Q_v,'r','LineWidth',2)
xlabel('Reaction Time(s)')
ylabel('Heat (J/m^3-s)')
title('Finite Rate Volumetric Reaction Heat Release')
%% Production Rate
figure(53)
plot2=plot(t,W_i,'LineWidth',2);
for i=1:no_of_species
set(plot2(i),'DisplayName',species(i).species);
end
xlabel('Reaction Time(s)')
ylabel('Rate (mol/cm^3-s)')
title('Producation Rate')
legend1 = legend('show');
set(legend1,'Location','best');
%% Density
figure(94)
plot(t,density,'LineWidth',2)
xlabel('Reaction Time(s)')
ylabel('Density (kg/m3)')
title('Density change')

%% Mass Flow rate calculation


[cp_final,Hi_final,Si_final]=nasa(a,temp,Tt_final,R);
Yi_final=(Xi_final.*M_i)./sum(Xi_final.*M_i)
mixture_enthalpy_final=sum((Yi_final.*Hi_final./M_i));

Q=mixture_enthalpy_final-mixture_enthalpy_initial % heat release J/kg-s

m_dot=(272190)/Q % mass flow rate
v=7.3 % velocity

density_initial=density(1,1) % Initial density
area=m_dot/(density_initial*v) % Area
d=sqrt(area*4/pi)       %diameter
mf=m_dot*Y_i


%% For data Extraction only
T=table(t,y(:,no_of_species+1),Q_v,X_it);