% This script is used to generate plots for LHS model outputs that incorporate sex hormones. 

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

subs = zeros(1,num_substances);

n = [4:13];

subs(n) = 1;

plot_concentration_bands(subs,median_sol_female,high_sol_female,low_sol_female,tspan,'Fem')

plot_concentration_bands(subs,median_sol_male,high_sol_male,low_sol_male,tspan,'Male')

