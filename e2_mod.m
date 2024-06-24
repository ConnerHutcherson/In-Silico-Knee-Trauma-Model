function [e2mod] = e2_mod(sampling_matrix, i, E2_min_conc, E2_max_conc, E2_median_conc, t_hormone)

% Function to perturb the input estrogen values according to the lhs matrix

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: Sampled estrogeb values (e2mod)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

for j=1:length(t_hormone)
    e2mod_below(j) = -sampling_matrix(i)*(E2_max_conc(j)-E2_min_conc(j))/2 + E2_median_conc(j);
    e2mod_above(j) = sampling_matrix(i)*(E2_max_conc(j)-E2_min_conc(j))/2 + E2_median_conc(j);
    e2_combined = cat(1,e2mod_above, e2mod_below);
        

end

e2mod = e2_combined(randi(2),:);