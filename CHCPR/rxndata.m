function [elm,spcs,reaction,atomic_weight,unbalanced_react_matrix,...
unbalanced_prod_matrix,rev,third_body,third_body_matrix,low,troe,sri,...
reaction_order_matrix,reverse_reaction_order_matrix,global_reactants,...
fuel_matrix,oxidizer_matrix,oxidizer_divide,fuel_divide]=rxndata(phi_on)

%% Extract the ELEMENTS, reacting SPECIES and REACTION MECHANISMS from  a
%%chemical reaction file. 'chem.inp'

% For extraction, the layout of the file should be as
%   ELEMENTS H O N END
%   SPECIES H2 O2 N2 OH END
%   REACTIONS 
%   H2+O2=OH+OH 1.73e+13 0 24070
%   END

%%   NOTE:THEY CAN VERTICALLY OR HORIZONALLY OR UNEVENLY ARRANGED BUT MUST
%   BE CONTAINED WITHIN RESPECTIVE HEADINGS with 'END' at end of each
%   headings and at least one space ' ' separating the entities.

%%Open 'chem.inp'file and extract all data into 'rdata'.
fid = fopen('chem.inp');


%% File Opening Check
if fid == -1
    disp('chem.inp File open not successful')
    error('chem.inp not found')
else
    disp('chem.inp File open successful')
    
%% Start scanning data
    rdata=textscan(fid,'%s','CommentStyle','%%');
    

%% For elements:
% PUTTING THE SPECIES FROM CELL 'rdata' TO STRUCTURE 'elm'
elm=struct('elements',zeros(120,1));
m=1;
k=1;
while strcmpi(rdata{1}{m},'END')~=1  
    if strcmpi(rdata{1}{m},'ELEMENTS')==1  
    m=m+1;
    end
    elm(k).elements=upper(strtrim(rdata{1}{m}));
    m=m+1;
    k=k+1;
end


%% FOR SPECIES:
%PUTTING THE SPECIES FROM CELL 'rdata' TO STRUCTURE 'spcs'
m=m+1;
k=1;
spcs=struct('species',zeros(500,1));
while strcmpi(rdata{1}{m},'END')~=1 
     if strcmpi(rdata{1}{m},'SPECIES')==1 
     m=m+1;
     end
     spcs(k).species=upper(rdata{1}{m});
     m=m+1;
k=k+1;
end 

%FOR REACTION:
%PUTTING THE REACTION MECHANISM, A FACTOR, N AND ENERGY FROM CELL 'rdata'
%TO STRUCTURE 'reaction' CONTAINING 'mechanism','Arf','nf','Ef'.
m=m+1;
k=1;
no_of_species=size(spcs);
reaction=struct('mechanism',zeros(1000,1),'Arf',zeros(1000,1),'nf',zeros(1000,1));
while strcmpi(rdata{1}{m},'END')~=1  
    if strcmpi(rdata{1}{m},'REACTIONS')==1
        m=m+1;
    end
    if strcmpi(rdata{1}{m},'DUP')==1 || strcmpi(rdata{1}{m},'DUPLICATE')==1
        m=m+1;
    end
    reaction(k).mechanism=upper(rdata{1}{m});
    reaction(k).Arf=str2double(rdata{1}{m+1});
    reaction(k).nf=str2double(rdata{1}{m+2});
    reaction(k).Ef=str2double(rdata{1}{m+3});
    

%% SPLIT REACTION MECHANISM INTO REACTANTS AND PRODUCTS FROM '='SIGN AND 
% STORE IN 'splitted_rxn'.Reactants are on (:,1) and in form 'H2+O2'
%and Products are on (:,2) and in form 'OH+OH' of cell 'spitted_rxn'.

% splitted_rxn{k} =strsplit(reaction(k).mechanism,{'=>','='}); 
splitted_rxn{k} =regexp(reaction(k).mechanism,'=','split');
irrev_match=regexp(reaction(k).mechanism,'=>', 'match');
if strcmp(irrev_match,'=>')==1
rev(k)=0;
else
rev(k)=1;
end


% THUS OBTAINED REACTANTS ARE FURTHER SPLITTED FROM '+' SIGN TO OBTAIN
% 'H2' AND 'O2' IN (:,1) AND (:,2) OF 'reactants'. 

% reactants{k}=strtrim(strsplit(splitted_rxn{k}{1},'+'));
reactants{k}=strtrim(regexp(splitted_rxn{k}{1},'+','split'));

%SIMILAR TO 'reactants'

% products{k}=strtrim(strsplit(splitted_rxn{k}{2},'+'));
products{k}=strtrim(regexp(splitted_rxn{k}{2},'+','split'));

%% For Matlab older than 2013 instead of above as strsplit is not defined.
% splitted_rxn{i} =regexp(reaction(k).mechanism,'=','split'); 
% products{i}=strtrim(regexp(splitted_rxn{k}{2},'+','split'));
% reactants{i}=strtrim(regexp(splitted_rxn{k}{1},'+','split'));
%% 


%% CREATES  A COEFFICIENT MATRIX FOR REACTANTS
%%ARRANGED ACCORDING TO SUPPLIED SPECIES
[unbalanced_react_matrix(k,:),third_body(k)]=matrixcreation(reactants{k},spcs);

%% CREATES  A COEFFICIENT MATRIX FOR PRODUCTS
%%ARRANGED ACCORDING TO SUPPLIED SPECIES    

unbalanced_prod_matrix(k,:)=matrixcreation(products{k},spcs);
m=m+4;

%%Arbitary reaction order "Forward"
if strcmpi(rdata{1}{m},'FODR')==1
    reaction_order_matrix(k,:)=zeros(1,no_of_species(2));
    m=m+1;
    reaction_order_input{k}=rdata{1}{m};
    reaction_order_split=strtrim(regexp(reaction_order_input{k},'/','split'));
    sff=size(reaction_order_split);
    for i=1:2:sff(2)
                for j=1:no_of_species(2)
                    if strcmpi(spcs(j).species,reaction_order_split(1,i))==1
                        reaction_order_matrix(k,j)=str2double(reaction_order_split(1,i+1));
                    end
                end
     end
m=m+1;
else
    reaction_order_matrix(k,:)=zeros(1,no_of_species(2));
end

%%Arbitary reaction order "Forward"
if strcmpi(rdata{1}{m},'RODR')==1
    reverse_reaction_order_matrix(k,:)=zeros(1,no_of_species(2));
    m=m+1;
    reverse_reaction_order_input{k}=rdata{1}{m};
    reverse_reaction_order_split=strtrim(regexp(reaction_order_input{k},'/','split'));
    sff=size(reaction_order_split);
    for i=1:2:sff(2)
                for j=1:no_of_species(2)
                    if strcmpi(spcs(j).species,reverse_reaction_order_split(1,i))==1
                        reverse_reaction_order_matrix(k,j)=str2double(reverse_reaction_order_split(1,i+1));
                    end
                end
     end
m=m+1;
else
    reverse_reaction_order_matrix(k,:)=zeros(1,no_of_species(2));
end

%% LOW Extraction
if strcmpi(rdata{1}{m},'LOW')==1
    low(k)=1;
    m=m+2;
    reaction(k).low_Arf= str2double(rdata{1}{m});
    reaction(k).low_nf=str2double(rdata{1}{m+1});
    reaction(k).low_Ef=str2double(rdata{1}{m+2});
    m=m+4;
else
    low(k)=0;
    reaction(k).low_Arf=0;
    reaction(k).low_nf=0;
    reaction(k).low_Ef=0;
end

%% Troe Extraction
if strcmpi(rdata{1}{m},'TROE')==1
    troe(k)=1;
    m=m+2;
    reaction(k).troe.alpha=str2double(rdata{1}{m});
    reaction(k).troe.T3=str2double(rdata{1}{m+1});
    reaction(k).troe.T1=str2double(rdata{1}{m+2});
    reaction(k).troe.T2=0;
    m=m+3;
    if strcmpi(rdata{1}{m+1},'/')==1
    reaction(k).troe.T2=str2double(rdata{1}{m});
    m=m+1;
    end
    m=m+1;
else
    troe(k)=0;
    reaction(k).troe=0;
end

%% SRI Extraction
    if strcmpi(rdata{1}{m},'SRI')==1
    sri(k)=1;
    m=m+2;
    reaction(k).sri.a=str2double(rdata{1}{m});
    reaction(k).sri.b=str2double(rdata{1}{m+1});
    reaction(k).sri.c=str2double(rdata{1}{m+2});
    reaction(k).sri.d=1;
    reaction(k).sri.e=0;
    m=m+3;
    if strcmpi(rdata{1}{m+2},'/')==1
    reaction(k).sri.d=str2double(rdata{1}{m});
    reaction(k).sri.e=str2double(rdata{1}{m+1});
    m=m+2;
    end
    m=m+1;
else
    sri(k)=0;
    reaction(k).sri=0;
    end

%% third_body_matrix=struct;

if third_body(k)==1
    third_body_on=1;
    third_body_matrix(k,:)=ones(1,no_of_species(2));
  while strcmpi(regexp(rdata{1}{m},'=','match','once'),'=')~=1 && strcmpi(rdata{1}{m},'END')~=1
    if strcmpi(regexp(rdata{1}{m},'/','match','once'),'/')==1
        third_body_combined{k}=rdata{1}{m};
        third_body_split=strtrim(regexp(third_body_combined{k},'/','split'));
        sf=size(third_body_split);
            for i=1:2:sf(2)
                for j=1:no_of_species(2)
                    if strcmpi(spcs(j).species,third_body_split(1,i))==1
                        third_body_matrix(k,j)=str2double(third_body_split(1,i+1));
                    end
                end
            end
    end
    m=m+1;
  end
else
    third_body_matrix(k,:)=zeros(1,no_of_species(2))  ;
end
    
k=k+1;
end


%% GLOBAL REACTION
size_rdata=size(rdata{1});
% global_reaction_on=0;
% for i=1:size_rdata(1) % no of data on chemical reaction mechanism file
%     if strcmpi(rdata{1}{i},'GLOBAL')==1
%         global_reaction_on=1;
%     end
% end

if phi_on==1
    count=0;
    for i=1:size_rdata(1)
    if strcmpi(rdata{1}{i},'GLOBAL')==1
        count=count+1;
    end
    end
    if count==0
        error('No global reaction defined');
    end
m=m+1;
while strcmpi(rdata{1}{m},'END')~=1  
    if strcmpi(rdata{1}{m},'GLOBAL')==1
        m=m+1;
    end
    globalreac=upper(rdata{1}{m});
    fuel=upper(rdata{1}{m+1});
    oxidizer=upper(rdata{1}{m+2});
    m=m+3;

end



%% Split Global Reaction

global_splitted_rxn =regexp(globalreac,'=','split'); 
global_reactants=strtrim(regexp(global_splitted_rxn{1,1},'+','split'));
global_reactants=matrixcreation(global_reactants,spcs);

global_fuel=strtrim(regexp(fuel,'/','split'));
global_oxidizer=strtrim(regexp(oxidizer,'/','split'));

sf=size(global_fuel);
fuel_matrix=zeros(1,no_of_species(2));
% fuel_matrix
for i=1:2:sf(2)
    count=0;
    for j=1:no_of_species(2)
        if strcmpi(spcs(j).species,global_fuel(1,i))==1
            fuel_matrix(1,j)=str2double(global_fuel(1,i+1));
            count=count+1;
        end
    end
%     if count==0
%         sprintf('Species - %s  ',spcs(i).species)
%         error('Species not found. Check for above species chemical mechanism file');
%     end
end

so=size(global_oxidizer);
oxidizer_matrix=zeros(1,no_of_species(2));
for i=1:2:so(2)
    count=0;
    for j=1:no_of_species(2)
        if strcmpi(spcs(j).species,global_oxidizer(1,i))==1
            oxidizer_matrix(1,j)=str2double(global_oxidizer(1,i+1));
            count=count+1;
        end
    end
%         if count==0
%         sprintf('Species - %s  ',spcs(i).species)
%         error('Species not found. Check for above species chemical mechanism file');
%         end
end
fuel_divide=[global_fuel(1,1) global_fuel(1,2)];
oxidizer_divide=[global_oxidizer(1,1) global_oxidizer(1,2)];
end


%% EXAMPLE IF SUPPLIED SPECIES ARE H H2 O2 OH 
%if H2+O2=OH+OH AND H+H=H2 ARE THE REACTION
%REACT MATRIX IS  [0 1 1 0;2 0 0 0]
%PROD_MATRIX IS [0 0 0 2; 0 1 0 0

%% File Close Check
closeresult = fclose(fid);
if closeresult == 0
disp('chem.inp File close successful')
else
disp('chem.inp File close not successful')
end

%% Extract the ELEMENTS and respective ATOMIC MASS from file
fid = fopen('mass_table.txt');
mdata=textscan(fid,'%s %s');
no_of_elements=size(elm);
p=size(mdata{1});
atomic_weight=zeros(no_of_elements(2),1);
for i=1:no_of_elements(2)  %no. of elements in input file  
    count=0;
    for j=1:p(1) %total no. of elements on mass_table file
        temp_data=strtrim(mdata{1,1}{j,1});
        if strcmpi(temp_data,elm(i).elements)==1 %%compare elements
          atomic_weight(i,1)=str2double(mdata{1,2}{j,1});
          count=count+1;
        end
    end
    if count==0
        sprintf('Element - %s  ',elm(i).elements)
        error('Element not found. Check for above element in atomic_mass file and chemical mechanism file');
    end

end  
clear ans count fid i j k m n p splitted_rxn rdata temp_data mdata oxidizer fuel


end



