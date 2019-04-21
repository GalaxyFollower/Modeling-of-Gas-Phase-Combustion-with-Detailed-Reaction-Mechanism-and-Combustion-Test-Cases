function [stoic_react_matrix,stoic_prod_matrix]= stoic_balance(react_matrix,prod_matrix,no_of_elements)

%%
m=size(react_matrix);
el=size(no_of_elements);
stoic_react_matrix=zeros(m(1),m(2));
stoic_prod_matrix=zeros(m(1),m(2));

for i=1:m(1)   % no.of reactions
%% Create matrix consisting total no of elements in a reaction 'elemental_matrix'.
% %total no of species in reaction
c1_matrix=repmat(react_matrix(i,:),el(2),1);
c2_matrix=repmat(prod_matrix(i,:),el(2),1);
c1_react=c1_matrix.*no_of_elements';
c2_prod=c2_matrix.*no_of_elements';
compare=sum(c1_react,2)==sum(c2_prod,2);
if compare==ones(el(2),1)
    stoic_react_matrix(i,:)=react_matrix(i,:);
    stoic_prod_matrix(i,:)= prod_matrix(i,:);
else
% %Create Repeated copies of array ex. if c_matrix=[1 1 1 0] and e_l(2)=3
%then c_matrix=[1 1 1 0;1 1 1 0; 1 1 1 0]. Done for elementwise
%multiplication with 'no_of_elements' matrix without using loop.
comp1=(react_matrix(i,:)-prod_matrix(i,:))>0;
comp=sum(comp1,2);
elemental_matrix =c1_react-c2_prod;

% % Gives the total no of elements of species present in given reaction
% ex. for reaction H2+O2+N2=OH+N2 with  input species H2 O2 OH N2
%elemental matrix=[2 0 1 0;0 2 1 0;0 0 0 4]
% elemental_matrix=no_of_elements'.*c_matrix
%% Balancing the reaction

% %Gives matrix after Gauss Jordan Elimination
% for above reaction gauss_eli_matrix=[1 0 0.5 0; 0 1 0.5 0; 0 0 0 0]
gauss_eli_matrix=rref(elemental_matrix);

% %Gives the stoichiometric balance reaction matrix with +sign for products
% and negaative for reactants
before_stoic_matrix=null(gauss_eli_matrix,'r')*2;
mm=size (before_stoic_matrix);

for nn=1:mm(2)    
B=1;
ccd=before_stoic_matrix(:,nn)>0;
if sum(ccd)>=comp
 kk = max(abs(before_stoic_matrix(:,nn)));% start at the end
 while B~=0 && kk>=0
     tmp = mod(before_stoic_matrix(:,nn),kk);
     B = sum(tmp(:));
     kk = kk - 1;
 end
kk = kk+1;
stoic_matrix = before_stoic_matrix(:,nn)/kk;
    
stoic_react_matrix(i,:)=abs(react_matrix(i,:).*stoic_matrix');
stoic_prod_matrix(i,:)=abs(prod_matrix(i,:).*stoic_matrix');
end
end

end
end