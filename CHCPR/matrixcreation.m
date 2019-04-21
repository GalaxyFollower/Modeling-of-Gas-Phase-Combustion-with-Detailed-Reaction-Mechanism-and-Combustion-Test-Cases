function [unbalanced_react_or_prod_matrix,third_body]=matrixcreation(react_or_prod,spcs)


%%% THIS FUNCTION SEPARATES THE COEFFICIENTS OF REACTANTS OR PRODUCTS FROM
%%% REACTANT OR PRODUCT SPECIES.


%% EXTRACTS THE COEFFICIENTS WEATHER INTEGER OR DECIMAL
react_or_prod_coeffs=regexp(react_or_prod,'^(\d*\.\d*|\d*)', 'match', 'once');

%% EXTRACTS THE SPECIES FROM REACTION
react_or_prod_species_matrix=strtrim(regexp(react_or_prod,'[A-Z]\w+|[A-Z]','match'));


%% IF THE COEFFICIENT IS NOT GIVEN ASSUMES '1' ELSE THE GIVEN COEFFICENT IN USED.
n=size(spcs);
p=size(react_or_prod_species_matrix);

for ff=1:p(2)
mpp=react_or_prod_coeffs{1,ff};
if strcmp(mpp,'')==1
    react_or_prod_coeff_matrix(1,ff)=1;
else
    react_or_prod_coeff_matrix(1,ff)=str2double(react_or_prod_coeffs{1,ff});
end
end

%%
unbalanced_react_or_prod_matrix=zeros(1,n(2));

%% CREATES  A COEFFICIENT MATRIX FOR REACTANTS OR PRODUCT
%%ARRANGED ACCORDING TO SUPPLIED SPECIES
third_body=0;
    p=size(react_or_prod_species_matrix);
for j=1: n(2) % no.of species
    for k= 1 : p(2) % no of reactant species of current reaction
        if strcmp(spcs(j).species,react_or_prod_species_matrix{k}{1})==1
            unbalanced_react_or_prod_matrix(1,j)=unbalanced_react_or_prod_matrix(1,j)+react_or_prod_coeff_matrix(1,k);
        end
        if strcmp('M',react_or_prod_species_matrix{k}{1})==1
            third_body=1;
        end
    end

end
    
    