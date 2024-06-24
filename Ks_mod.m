function [Ks P] = Ks_mod(sampling_matrix, i,range_factor)

% Function to perturb the input K-values according to the lhs matrix

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: Sampled input parameter values (Ks and P)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


Ks = rate_coeffs;

Ks_min = Ks-Ks.*range_factor;
Ks_max = Ks+Ks.*range_factor;

for j=1:length(Ks)
    Ks(j) = sampling_matrix(i,j)*(Ks_max(j)-Ks_min(j))+Ks_min(j);
    P(:,j) = sampling_matrix(i,j)*(Ks_max(j)-Ks_min(j))+Ks_min(j);
end