function [prmod] = pr_mod(sampling_matrix, i, P4_min_conc,P4_max_conc,P4_median_conc, t_hormone)

% Function to perturb the input pr-values according to the lhs matrix

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: Sampled progesterone parameter values (prmod)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


for j=1:length(t_hormone)
    prmod_below(j) = -sampling_matrix(i)*(P4_max_conc(j)-P4_min_conc(j))/2+P4_median_conc(j);
    prmod_above(j) = sampling_matrix(i)*(P4_max_conc(j)-P4_min_conc(j))/2+P4_median_conc(j);
    pr_combined = cat(1,prmod_above, prmod_below);
        

end

prmod = pr_combined(randi(2),:);
   
end


