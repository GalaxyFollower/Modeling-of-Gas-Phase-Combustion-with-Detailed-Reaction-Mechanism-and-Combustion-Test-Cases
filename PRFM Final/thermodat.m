function [coeff,ttemp,no_of_elements] = thermodat(spcs,elements)

%%opening file 'thermo.dat'%%
fid = fopen('thermo.dat', 'r' );


%% File Opening Check
if fid == -1
disp('thermo.dat File open not successful')
else
disp('thermo.dat File open successful')


%% Extract the thermodynamic curve fitting coefficients, temperature range
%for different elements from thermodynamic data file

% For extraction, the layout of the file should be as
%   LINE 1:
%   name =upto 16 Alphanumeric and character 1 to 16
%   low temperature = character 46-55
%   low temperature = character 56-65
%   low temperature = character 46-75
%   LINE 2:
%   first 5 coefficient for range [mid,high] = each of length 15 character 
%i.e from 1 to 75
%   LINE 3:
%   last 2 coefficient for range [mid,high] and first 3 coefficient for
%   range [low,mid] = length 15 character each from 1 to 75
%   LINE 4:
%   last 2 coefficient for range [low,mid]= length 15 character each from 1
%   to 60
% 'END' at the end of thermo.dat file from character 1:3

%% defining structure 'ttdata' with elements 'name' having size [1000,1]%%
ttdata=struct('name',zeros(1000,1));

format short e;
i=1;

%% Reading first 4 line of 'thermo.dat'%%
tdata1=fgetl(fid);
tdata2=fgetl(fid);
tdata3=fgetl(fid);
tdata4=fgetl(fid);

%% assigning the values from thermo.dat to structure 'ttdata' until program
%finds 'END'. 

while strcmp(tdata1(1:3),'END')~=1
    
ttdata(i).name=strtok(strtrim(tdata1(1:19)));      % for element name
ttdata(i).elems=tdata1(25:44);
ttdata(i).phase=tdata1(45);
ttdata(i).temp(1)=str2double(tdata1(46:55));   % low temperature limit
ttdata(i).temp(2)=str2double(tdata1(56:65));  % high temperature limit
ttdata(i).temp(3)=str2double(tdata1(66:75));  %mid temperature  %mid temperature


%% assigining first 5 coefficients for range(mid,high)
n=1;
for k=1:5
        ttdata(i).a(k)=str2double(tdata2(n:n+14));
        n=n+15;
end


%% assigining last 2 coefficients for range(mid,high) and first 3
%coefficients for temperature range (low,mid)
n=1;
for k=6:10
    ttdata(i).a(k)=str2double(tdata3(n:n+14));
    n=n+15;  
end


%% assigining last 4 coefficients for range (low,mid)n=1;
n=1;
for k=11:14
    ttdata(i).a(k)=str2double(tdata4(n:n+14));
    n=n+15;
end  

%% reading next 4 lines
tdata1=fgetl(fid);
tdata2=fgetl(fid);
tdata3=fgetl(fid);
tdata4=fgetl(fid);
i=i+1;

end

%% File Close Check
closeresult = fclose(fid);
if closeresult == 0
disp('thermo.dat File close successful')
else
disp('thermo.dat File close not successful')
end
end


% spcs={'H2','O2','N2','OH'};
%% Compare the species with the datas in 'ttdata.name' and extract the 
%fitting coefficients,temperature range,phase,elemental data.
n=size(spcs);
m=size(ttdata);
p=size(elements);
coeff=zeros(n(2),14);
phase=char(n(2),1);
ttemp=zeros(n(2),3);

for i=1:n(2)
    for j=1:m(2)
        if strcmp(spcs(i).species,ttdata(j).name)==1
                ttemp(i,:)=ttdata(j).temp;
                coeff(i,:)=ttdata(j).a;
                elemns(i,:)=ttdata(j).elems;
                phase(i,:)=ttdata(j).phase;
        end
    end
end

%% Arrange the elemental data collected above according to input elements
% also find the no.of each elements on the input species.

no_of_elements=zeros(n(2),p(2));
    for i=1:n(2)
        for j=1:p(2)
            s=1;
            te=0;
            for k=1:4 
                elements_th=strtrim(elemns(i,s:s+1));
                if strcmp(elements(j).elements,elements_th)==1   
                     te= str2double(elemns(i,s+2:s+4));
                end
                no_of_elements(i,j)=te;
                s=s+5;
            end
        end
    end         
    
    
    
clear ttdata tdata1 tdata2 tdata3 tdata4 i j k l m s fid elems te p n 
clear elements_th



%% Example : Species H H2 O2 N2 
% Coefficient matrix is [5*15]matrix
%H: [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a13,a14,a15]
%H2: same as above
%02:
%N2:
%OH:
end
