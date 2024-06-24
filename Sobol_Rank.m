% This script will Rank the top 5 Sobol sensativity indicies (Sii or STii) generated from the
% LHS simulations. 

% Input: Sobol sensativity (Sii or STii) cell array from main script

% Output: Sensativity_Rank = rank all parameters for each substance during
% simulation
% Top_Rank = rank top 5 parameters for each substance during
% simulation

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

tic

% Identify negative sensativities and assume sesnsativity value is insignificant or equal to 0 if 95% CI interval of sensativity values includes 0.
% Sensativity matricies should be checked for negative values prior to ranking. Refer to Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., Gatelli, D., ... & Tarantola, S. (2008).
% Global sensitivity analysis: The primer, to troubleshoot the origin of large negative indicie values. 

for idx =2:numel(Sii)
   Sii{idx}(Sii{idx}<0)=0; %uncomment if ranking total sensativity values from model 
%  STii{idx}(STii{idx}<0)=0; %uncomment if ranking total sensativity values from model  
  
end

Sensativity_Rank = cell(1,10);
Top_Rank = cell(1,10);
Top_Sensativity = cell(1,10);

% Rank Input Sensativities for each substances (highest rank = 49 , lowest
% rank =1)
for i=2:numel(Sii)
    for k=1:100
    Sensativity_Rank{i-1}(k,:)=tiedrank(Sii{i}(k+1,:));  
    %Find Top 5 Ranks from each time increment
    [~,I]=maxk(Sensativity_Rank{i-1}(k,:),5);
    [B,~]=maxk(Sii{i}(k+1,:),5);
    
    %Ranks stored in cell array for each molecule (column: Rank x row:
    %Scaled Time)
    Top_Rank{i-1}(k,:)= I;
    Top_Sensativity{i-1}(k,:)= B;
    end  
end

toc
save('Ranked_Sensativity_Analysis')