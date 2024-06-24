% For each simulation day (0-30), conduct Rank Based Covariate (RBC) analysis on model results from  "Female
% Peak E," "Female Low E," and "Male" differ.  

% To run this script, you must import data files all_sols_t, all_sols_female_day15,
% and all_sols_day0 after running from main_script.m for each respected
% simulation.

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2024
% mail to: conner.hutcherson@utsouthwestern.edu

post_hoc{1} = ones(101,4); % IL-1
post_hoc{2} = ones(101,4); % TNF-a
post_hoc{3} = ones(101,4); % IL-10
post_hoc{4} = ones(101,4); % TGFb
post_hoc{5} = ones(101,4); % MMP-9
post_hoc{6} = ones(101,4); % MMP-1
post_hoc{7} = ones(101,4);% TIMP-1
post_hoc{8} = ones(101,4);% IL-6
post_hoc{9} = ones(101,4);% MMP-13
post_hoc{10} = ones(101,4);% MMP-3

alpha = 0.05;

inds = [4:13];
count=0;

for i=1:length(tspan) % Loop through days (each column of all_sols corresponds to a time)
    count = count+1;
    % Extract data from solution cells (all_sols_...) and covariates (hormones) 
    
estrogenMale = e2_vary_male(:,i); %  Estrogen concentrations for the Male model
progesteroneMale = pr_vary_male(:,i); % Progesterone concentrations for the Male model
testosteroneMale = test_vary_male(:,i); % Testosterone for the Male model 

estrogenFemaleMenses = e2_vary_female_menses(:,i); %  Estrogen concentrations for the Female model at cycle day 1
progesteroneFemaleMenses = pr_vary_female_menses(:,i); % Progesterone concentrations for the Female model at cycle day 1
testosteroneFemaleMenses = test_vary_female_menses(:,i);% Testosterone for for the Female model at cycle day 1

estrogenFemalePreOV = e2_vary_female_preov(:,i);%  Estrogen concentrations for the Female model at cycle day 15
progesteronefemalePreOV = pr_vary_female_preov(:,i);% Progesterone concentrations for the Female model at cycle day 15
testosteroneFemalePreOV = test_vary_female_preov(:,i);% Testosterone for for the Female model at cycle day 15

% Combine covariate values from all groups into a single array for analysis

Estrogen = [estrogenFemalePreOV; estrogenFemaleMenses; estrogenMale]; % Combine values from all groups into a single array
Progesterone = [progesteronefemalePreOV; progesteroneFemaleMenses; progesteroneMale];
Testosterone = [testosteroneFemalePreOV; testosteroneFemaleMenses; testosteroneMale];

    for j=1:10 % Loop through substances for each group
        female_preov_output = all_sols_female_preov{inds(j)}(:,i);
        female_menses_output = all_sols_female_menses{inds(j)}(:,i);
        male_output = all_sols_male{inds(j)}(:,i);
        
        % Combine model output values for each substance from all groups into a single array for analysis
        Substance_Output = [female_preov_output; female_menses_output; male_output];
        
        % Assign numerical value for Groups
        Female_Preov_Group = ones(1000,1);% Pre-Ov Group 1
        Female_Menses_Group = ones(1000,1)*2; % Menses Group 2
        Male_Group = ones(1000,1)*3; %Male Group 3
        
        % Combine group values into a single array for analysis
        Group = [Female_Preov_Group; Female_Menses_Group; Male_Group];
        
        % Function to run ranked based covariate adjustment (RBC) analysis
        [p, F, df1, df2, residuals] = RBCtest(Substance_Output, Estrogen, Progesterone, Testosterone, Group);
 
        if p<0.05

%           % Perform Dunn procedure for multiple non parametric comparisons.
            % Q-values for the post-hoc tests can be found in the array
            % "post_hoc" 
            
            % Female Pre OV (Peak E) vs Menses (Low E) (column 2)
            % Female Pre OV (Peak E) vs Male (High T) (column 3)
            % Female Menses (Low E) vs Male (High T) (column 4)
% 
            [Q] = dunn2(residuals',Group');
%             [post_hoc{j}{count}] = dunn2(residuals',Group');
% 
           for k = 1:3
              
               if strcmp(Q{k,1}, '1-2') || strcmp(Q{k,1}, '2-1')
                  post_hoc{j}(count,2) = 2*(1-normcdf(Q{k,2})); 

               elseif strcmp(Q{k,1}, '1-3') || strcmp(Q{k,1}, '3-1')
                      post_hoc{j}(count,3) = 2*(1-normcdf(Q{k,2}));   

               else  
                       post_hoc{j}(count,4) = 2*(1-normcdf(Q{k,2})); 
               end
              
           end

        end

    end 

 end  
        








