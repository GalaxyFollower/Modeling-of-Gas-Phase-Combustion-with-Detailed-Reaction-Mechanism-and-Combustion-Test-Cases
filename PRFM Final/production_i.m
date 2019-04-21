
function [w_i,cp_i,H_i,kf] = production_i(Tt,P,rho_i,M_i,vr,vp,Arf_k,nf_k,Ef_k,a,temp,react_matrix,prod_matrix,rev,third_body,third_body_matrix,~)


 R=8.310;                        %J/mol K

X_iconc = (rho_i./M_i)/1e6;

%% Specific heat, enthalpy and entropy at temperature T_t (J-mol-sec-K)
    [cp_i,H_i,S_i]=nasa(a,temp,Tt,R);
    

%% Forward Reaction Rate Constant
    kf=(Arf_k.*power(Tt,nf_k).*exp(-(Ef_k/(R*Tt))));
%   kf=exp(log(Arf)+nf*log(Tt)-(Ef/(R*Tt)));
    [no_of_reactions,no_of_species]=size(react_matrix);

%% For Backward Reaction Rate Constant 
%     S_i=S_i-R*log(P/Patm)*ones(no_of_species,1);

%     S_i_mix=log(X_i);
%     cond=isinf(S_i_mix);
%     for i=1:no_of_species
%         if cond(i)==1
%             S_i_mix(i)=0;
%         end
%     end

    del_G=(prod_matrix-react_matrix)*(H_i-Tt*S_i);
    
    % Equilibrium constant
    Kc=((101325)/(8.3144598e6*Tt)).^(sum(prod_matrix-react_matrix,2)).*exp(-del_G/(R*Tt));
    kb=rev'.*(kf./Kc);
    
    
%     Kc=(Patm/(P*R*Tt)).^sum(prod_matrix-react_matrix,2).*exp(-(Ef_k-del_H_i)/(R*Tt));
    
    kf=kf';
    kb=kb';
    
%%  Rate of Progress Variable(q_i)
 
    progress_i_react=zeros(size(react_matrix))';
    progress_i_prod=zeros(size(prod_matrix))';
    for i=1:no_of_reactions 
      progress_i_react(:,i)=power(X_iconc,vr(:,i));
      progress_i_prod(:,i)=power(X_iconc,vp(:,i));
      if third_body(i)==1
%           p=third_body_matrix(i,:)*X_iconc
         c_third(1,i)=third_body_matrix(i,:)*X_iconc;
     else
         c_third(1,i)=1;
     end
    end
    q_if=kf.*prod(progress_i_react).*c_third;
    q_ir=rev.*kb.*prod(progress_i_prod).*c_third;
%     q_i=kf.*prod(progress_i_react)-rev.*kb.*prod(progress_i_prod);
%% Production Rate for forward and reverse reaction
     w_if=(q_if*(prod_matrix-react_matrix))' ; % mol/cm^3-s Forward
     w_ir=(q_ir*(prod_matrix-react_matrix))'; % Reverse 
     
     %% Producation Rate (w_i)
     w_i=w_if-w_ir;
     w_i=real(w_i*1e6)
     
end
     