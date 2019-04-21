function preprocessor_plots(reaction,spcs,a,temp,react_matrix,prod_matrix,rev,R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to compute & plot the reaction rates and thermodynamic data with 
% temperature.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                            DEFINITIONS                            %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Same as in main file
% T=temperature


%% For Plotting
[no_of_reactions,no_of_species]=size(react_matrix);

Arf_k=[reaction(:).Arf]';
nf_k=[reaction(:).nf]';
Ef_k=[reaction(:).Ef]';
k=1;

%% Compute forward rate, backward rate and equilibrium constant (Kc) from
% temperature range in thermodata file
for T=temp(1,1):100:temp(1,2)
    % Forward Reaction Rate
    kf=(Arf_k.*power(T,nf_k).*exp(-(Ef_k/(R*T))));
    
    % thermodynamic datas at temperature T
    [cp_i,H_i,S_i]=nasa(a,temp,T,R);
    
    enthalpy(:,k)=(H_i)/1000;
    entropy(:,k)=S_i;
    specific_heat(:,k)=cp_i/R;

    % Calculates gibb's free energy of reaction
    del_G=(prod_matrix-react_matrix)*(H_i-T*S_i);
    
    % Equilibrium Constant
    Kc=(101325/(8.3144598e6*T)).^(sum(prod_matrix-react_matrix,2)).*exp(-del_G/(R*T));
    
    % Backward Reaction Rate
    kb=rev'.*(kf./Kc);

    T1(1,k)=T;
    Te(1,k)=1/T;
    kff(:,k)=log10(kf);
    kbb(:,k)=log10(kb);
    kc(:,k)=Kc;
    k=k+1;
end


%% Plotting Thermochemical Data
for i=1:no_of_species
name=spcs(i).species;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'on');

% Create multiple lines using matrix input to plot
% plot1=plot(T1,enthalpy(i,:),T1,entropy(1,:),T1,specific_heat(1,:));
plot1 = plot(T1,enthalpy(i,:),T1,entropy(i,:),T1,specific_heat(i,:),'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','H (kcal/mol)');
set(plot1(2),'DisplayName','S (cal/mol)');
set(plot1(3),'DisplayName','cp (cal/mol-K)');

% Create xlabel
xlabel('T (K)');

% Create ylabel
ylabel('H,S,cp');

% Create title
title(name);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','best');

end

%% Plotting Reaction Rates
for i=1:no_of_reactions
 name=reaction(i).mechanism;
% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'on');
grid(axes1,'on');
% Create multiple lines using matrix input to plot
plot1 = plot(Te,kff(i,:),Te,kbb(i,:),'LineWidth',2);
set(plot1(1),'DisplayName','Forward');
set(plot1(2),'DisplayName','Reverse');

% Create xlabel
xlabel('1/T');

% Create ylabel
ylabel('log_1_0(k)');

% Create title
title(name);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','best');
end
end