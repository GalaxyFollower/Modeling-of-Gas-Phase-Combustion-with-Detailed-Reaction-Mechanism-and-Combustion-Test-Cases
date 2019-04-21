%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATES density OF EACH SPECIES FOR GIVEN CONDITIONS.
%% Paramaters
% X_i=mole fraction
% M_i=Moleular weight in gram
% R=universal gas costant j/mol-K
% P= pressure (pa)

%% Calculation
function density=densityvolume(X_i,M_i,R,Tt,P,X_iconc)
density=P*sum(X_i.*M_i)/(R*Tt);
% de=sum(X_iconc.*M_i*10^6)
end