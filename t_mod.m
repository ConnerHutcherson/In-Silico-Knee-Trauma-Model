function [tmod] = t_mod(sampling_matrix, i ,T_min_conc,T_max_conc,T_median_conc, t_hormone)

% Function to perturb the input testosterone values according to the lhs matrix

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: Sampled testosterone  values (tmod)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

for j=1:length(t_hormone)    
    tmod_below(j) = -sampling_matrix(i)*(T_max_conc(j)-T_min_conc(j))/2+T_median_conc(j);
    tmod_above(j) = sampling_matrix(i)*(T_max_conc(j)-T_min_conc(j))/2+T_median_conc(j);
    t_combined = cat(1,tmod_above, tmod_below);
        

end

tmod = t_combined(randi(2),:);
    
   
end

