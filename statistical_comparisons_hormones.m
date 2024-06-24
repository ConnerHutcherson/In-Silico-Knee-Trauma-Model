% For each simulation day (0-30), do Kruskal-Wallis test to see if "Female
% Peak E," "Female Low E," and "Male" differ.  

% To run this script, you must import data files all_sols_t, all_sols_female_day15,
% and all_sols_day0 after running from main_script.m for each respected
% simulation.

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

post_hoc{1} = ones(101,4); % IL-1
post_hoc{2} = ones(101,4); % TNF-a
post_hoc{3} = ones(101,4); % IL-10
post_hoc{4} = ones(101,4); % MMP-9
post_hoc{5} = ones(101,4); % MMP-1
post_hoc{6} = ones(101,4); % MMP-3


for h=1:4
    post_hoc{h}(:,1)=1:101;
end

alpha = 0.0167;

inds = [4,5,6,8,9,13];
count=0;

for i=1:101 % Loop through days (each column of all_sols corresponds to a time)
    count = count+1;
    % Extract data from solution cells (all_sols_...)

    for j=1:6 % Loop through substances (IL-1, TNF-a, IL-10, MMP-9, MMP-1, MMP-13, MMP-3) 
        e2p = all_sols_female_day15{inds(j)}(:,i);
        e2l = all_sols_female_day0{inds(j)}(:,i);
        te = all_sols_t{inds(j)}(:,i);
        
       
        all_sols_sex_test = [e2p, e2l, te];

%         Kruskal-Wallis
%   
        [p,tbl,stats] = kruskalwallis(all_sols_sex_test,[],'off');
        
        % If significant differences, do pairwise Mann-Whitney U tests to see which
        % differences are significant (Bonferroni correction: alpha =
        % 0.0167) One-sided tests
       if p<0.05
%             
            % P-values for the post-hoc tests can be found in the array
            % "post_hoc" 
            
            % Female Pre OV (Peak E) and Menses (Low E) (columns 1 and 2)
            % Perform Mann-Whitney U Test
            post_hoc{j}(count,2)=ranksum(e2p,e2l,'alpha',alpha);

            % Female Menses (Low E) and Male (columns 2 and 3)
            % Perform Mann-Whitney U Test
            post_hoc{j}(count,3)=ranksum(e2l,te,'alpha',alpha);%,'tail','right');


            % Female Pre OV and Male (columns 1 and 3)
            % Perform Mann-Whitney U Test
            post_hoc{j}(count,4)=ranksum(e2p,te,'alpha',alpha);%,'tail','right');
%         

        end        
       
    end
    

    
end





