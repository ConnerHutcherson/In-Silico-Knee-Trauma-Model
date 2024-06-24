% For each simulation day (0-30), do Rank Based Covariate (RBC) analysis on
% resuts from female cycle simulations. Day 0, Day 7, Day 15, Day 18, and
% Day 26
% 
% To run this script, you must import data files  after running from main_script.m for each respected
% simulation.

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2024
% mail to: conner.hutcherson@utsouthwestern.edu

post_hoc{1} = ones(10,4); % IL-1
post_hoc{2} = ones(10,4); % TNF-a
post_hoc{3} = ones(10,4); % IL-10
post_hoc{4} = ones(10,4); % TGFb
post_hoc{5} = ones(10,4); % MMP-9
post_hoc{6} = ones(10,4); % MMP-1
post_hoc{7} = ones(10,4);% TIMP-1
post_hoc{8} = ones(10,4);% IL-6
post_hoc{9} = ones(10,4);% MMP-13
post_hoc{10} = ones(10,4);% MMP-3

alpha = 0.05;
Model_Iterations = 1000;
inds = [4:13];
count=0;

female_Day0_output = zeros(Model_Iterations,1);
female_Day7_output = zeros(Model_Iterations,1);
female_Day15_output = zeros(Model_Iterations,1);
female_Day18_output = zeros(Model_Iterations,1);
female_Day26_output = zeros(Model_Iterations,1);

estrogenDay0 = zeros(Model_Iterations,1);
estrogenDay7 = zeros(Model_Iterations,1);
estrogenDay15 = zeros(Model_Iterations,1);
estrogenDay18 = zeros(Model_Iterations,1);
estrogenDay26 = zeros(Model_Iterations,1);

progesteroneDay0 = zeros(Model_Iterations,1); % Progesterone concentrations for each female group
progesteroneDay7 = zeros(Model_Iterations,1);
progesteroneDay15 = zeros(Model_Iterations,1);
progesteroneDay18 = zeros(Model_Iterations,1);
progesteroneDay26 = zeros(Model_Iterations,1);

testosteroneDay0 = zeros(Model_Iterations,1); % Testosterone for for each female group
testosteroneDay7 = zeros(Model_Iterations,1);
testosteroneDay15 = zeros(Model_Iterations,1);
testosteroneDay18 = zeros(Model_Iterations,1);
testosteroneDay26 = zeros(Model_Iterations,1);


for j=1:10
    for i = 1: Model_Iterations % Loop through substances for each group
        [female_Day0_output(i,:), time1] = max(all_sols_female_day0{inds(j)}(i,:));
        [female_Day7_output(i,:), time2] = max(all_sols_female_day7{inds(j)}(i,:));
        [female_Day15_output(i,:),time3] = max(all_sols_female_day15{inds(j)}(i,:));
        [female_Day18_output(i,:),time4] = max(all_sols_female_day18{inds(j)}(i,:));
        [female_Day26_output(i,:),time5] = max(all_sols_female_day26{inds(j)}(i,:));
        
        estrogenDay0(i,:) = e2_vary_day0(i,time1); %  Estrogen concentrations for each female group
        estrogenDay7(i,:) = e2_vary_day7(i,time2);
        estrogenDay15(i,:) = e2_vary_day15(i,time3);
        estrogenDay18(i,:) = e2_vary_day18(i,time4);
        estrogenDay26(i,:) = e2_vary_day26(i,time5);

        progesteroneDay0(i,:) = pr_vary_day0(i,time1); % Progesterone concentrations for each female group
        progesteroneDay7(i,:) = pr_vary_day7(i,time2);
        progesteroneDay15(i,:) = pr_vary_day15(i,time3);
        progesteroneDay18(i,:) = pr_vary_day18(i,time4);
        progesteroneDay26(i,:) = pr_vary_day26(i,time5);

        testosteroneDay0(i,:) = test_vary_day0(i,time1); % Testosterone for for each female group
        testosteroneDay7(i,:) = test_vary_day7(i,time2);
        testosteroneDay15(i,:) = test_vary_day15(i,time3);
        testosteroneDay18(i,:) = test_vary_day18(i,time4);
        testosteroneDay26(i,:) = test_vary_day26(i,time5);
        
    end 
 Substance_Output = [female_Day0_output; female_Day7_output; female_Day15_output; female_Day18_output; female_Day26_output];
 Estrogen = [estrogenDay0; estrogenDay7; estrogenDay15; estrogenDay18; estrogenDay26]; % Combine values from all groups into a single array
 Progesterone = [progesteroneDay0; progesteroneDay7; progesteroneDay15; progesteroneDay18; progesteroneDay26];
 Testosterone = [testosteroneDay0; testosteroneDay7; testosteroneDay15; testosteroneDay18; testosteroneDay26];
        
 % Assign numerical value for Groups
        Female_Day0_Group = ones(Model_Iterations,1);% Day 0
        Female_Day7_Group = ones(Model_Iterations,1)*2; % Day 7
        Female_Day15_Group = ones(Model_Iterations,1)*3; %
        Female_Day18_Group = ones(Model_Iterations,1)*4; %
        Female_Day26_Group = ones(Model_Iterations,1)*5; %
        
        % Combine group values into a single array for analysis
        Group = [Female_Day0_Group; Female_Day7_Group; Female_Day15_Group; Female_Day18_Group;Female_Day26_Group];
 
 % Function to run ranked based covariate adjustment (RBC) analysis
 [p, F, df1, df2, residuals] = RBCtest(Substance_Output, Estrogen, Progesterone, Testosterone, Group);
 
   if p<0.05
 % P-values for the post-hoc tests can be found in the array
            % "post_hoc" 
            
            % Female Pre OV (Peak E) and Menses (Low E) (columns 1 and 2)
            % Perform Dunn procedure for multiple non parametric comparisons.
            [post_hoc{j}] = dunn2(residuals',Group');
%             post_hoc{j}= 2*(1-normcdf(Q{:,2}));  % convert q value to p value      
   end
                
end
        
        
       
    






