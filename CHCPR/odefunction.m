function [dydt,RE] = odefunction(t,y,reaction,react_matrix,...
prod_matrix,rev,third_body,third_body_matrix,vr,vp,a,temp,R,P,M_i)
%% Start of Loop
[no_of_reactions,no_of_species]=size(prod_matrix);
y(1:no_of_species);
% X_iconc=y(1:no_of_species);
%% Specific heat, enthalpy and entropy at temperature T_t (J-mol-sec-K)
[cp_i,H_i,S_i]=nasa(a,temp,y(no_of_species+1),R);
%       h_mix=sum((Y_i.*(H_i-H_298)./M_i));
%       s_mix=sum(Y_i.*(S_i-S_298)./M_i);
   
Arf_k=[reaction(:).Arf]';
nf_k=[reaction(:).nf]';
Ef_k=[reaction(:).Ef]'.*4.1839999686;

% Arf_k_low=[reaction(:).low_Arf]';
% nf_k_low=[reaction(:).low_nf]';
% Ef_k_low=[reaction(:).low_Ef]'; 

    

%% Forward Reaction Rate Constant

% Default or high limit reaction rate coefficient
    kf=(Arf_k.*power(y(no_of_species+1),nf_k).*exp(-(Ef_k/(R*y(no_of_species+1)))));
% Low limit reaction rate coefficients    
%   kf_low=(Arf_k_low.*power(Tt,nf_k_low).*exp(-(Ef_k_low/(R*Tt))));


%% For Backward Reaction Rate Constant 
%     S_i=S_i-R*log(P/Patm)*ones(no_of_species,1);

%     S_i_mix=log(X_i);
%     cond=isinf(S_i_mix);
%     for i=1:no_of_species
%         if cond(i)==1
%             S_i_mix(i)=0;
%         end
%     end

    del_G=(prod_matrix-react_matrix)*(H_i-y(no_of_species+1)*S_i);
    
    % Equilibrium constant
    Kc=((101325)/(8.3144598e6*y(no_of_species+1))).^(sum(prod_matrix-react_matrix,2)).*exp(-del_G/(R*y(no_of_species+1)));
    kb=rev'.*(kf./Kc);
    
    kf=kf';
    kb=kb';
    
%%  Rate of Progress Variable(q_i)
    c_third=zeros(1,no_of_reactions);
    progress_i_react=zeros(size(react_matrix))';
    progress_i_prod=zeros(size(prod_matrix))';
    for i=1:no_of_reactions 
      progress_i_react(:,i)=power(y(1:no_of_species),vr(:,i));
      progress_i_prod(:,i)=power(y(1:no_of_species),vp(:,i));
      
      % Third Body Concentration Calculation
      if third_body(i)==1
         c_third(1,i)=third_body_matrix(i,:)*y(1:no_of_species);
      else
         c_third(1,i)=1;
      end
%       % When there is pressure dependance
%       if low(i)==1
%           Pr_i=c_third(1,i)*kf_low(i)/kf(i); % reduced pressure
%       if troe(i)==1
%           % For TROE Formuation 
%           F_i=troe_calculation(Tt,reaction(i).troe,Pr_i);
%        
%       elseif sri(i)==1
%           % For SRI Formulation
%           F_i=sri_calculation(Tt,reaction(i).sri,Pr_i);
%       else
%           % For Lindemann Formulation
%           F_i=1;
%       end
%       c_i=(Pr_i/(1+Pr_i))*F_i;
%       c_third(1,i)=c_i;
%       end
    end
    
    
    q_if=kf.*prod(progress_i_react).*c_third;
    q_ir=kb.*prod(progress_i_prod).*c_third;

    %% Production Rate for forward and reverse reaction
     w_if=(q_if*(prod_matrix-react_matrix))' ; % mol/cm^3-s Forward
     w_ir=(q_ir*(prod_matrix-react_matrix))'; % Reverse 
     c_i=(c_third.*kb.*prod(progress_i_prod)*react_matrix+c_third.*kf.*prod(progress_i_react)*prod_matrix);
     d_i=(c_third.*kf.*prod(progress_i_react)*react_matrix+c_third.*kb.*prod(progress_i_prod)*prod_matrix);
     
     %% Producation Rate (w_i)
     w_i=w_if-w_ir;
%       w_i=c_i-d_i;
   
     %% Volumetric Heat Rate   
      Q_vol =-sum(H_i.*w_i*10^6);
      
     %% Density  
       density=sum(y(1:no_of_species).*M_i*10^6);
     %% System of equations
      % Temperature
         dydt = zeros(no_of_species+1,1);
         delT = (-sum(H_i.*w_i)/sum(y(1:no_of_species).*cp_i));
%        delT_f = (-sum(H_i.*c_i')/sum(y(1:no_of_species).*cp_i));
%        delT_r = (-sum(H_i.*d_i')/sum(y(1:no_of_species).*cp_i));
      % Concentration  
         delX=w_i-y(1:no_of_species)*((sum(w_i)/sum(y(1:no_of_species)))+(delT/y(no_of_species+1)));
%        delX_f=c_i'-y(1:no_of_species)*((sumd(c_i)/sum(y(1:no_of_species)))+(delT_f/y(no_of_species+1)));
%        delX_r=d_i'-y(1:no_of_species)*((sum(d_i)/sum(y(1:no_of_species)))+(delT_r/y(no_of_species+1)));
%        dydt(1:no_of_species) =delX_f-delX_r;
%        dydt(no_of_species+1) =delT_f-delT_r;
         dydt(1:no_of_species)=delX;
         dydt(no_of_species+1)=delT;
      
     %% Post Processing
      RE = {w_i' kf' kb' Q_vol density};    
end
      