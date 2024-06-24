function [abcs] = abcs_mod(sampling_matrix, i,range_factor)

% Function to perturb the input abc feedback-values according to the lhs matrix

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: Sampled feedback parameter values (abcs)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


abcs = all_parameter_fits;

abcs_min = abcs-abcs.*range_factor;
abcs_max = abcs+abcs.*range_factor;

for j=1:length(abcs)
    abcs(j) = sampling_matrix(i,j)*(abcs_max(j)-abcs_min(j))+abcs_min(j);
end
end