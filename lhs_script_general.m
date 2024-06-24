function [sol_nominal,high_sol,low_sol,median_sol] = lhs_script_general(num_substances,tspan,num_coeffs)

% This script performs latin hypercube sampling for the rate coefficients. 
% It calls "Ks_mod" to recalculate parameter sets with the sampling matrix
% It stores n_sets of solutions, calculates the median and 1st and 3rd
% quartiles

% Use with "plot_confidence_bands" 

% Input: sampling matrix for model inputs (sampling_matrix)
% Range of variation of inputs (range_factor)

% Output: model outputs for lower sample range (low_sol, 25% quartile), higher range
% (high_sol, 75% quartile) and median (median_sol)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

global Ks e2 estr0 pr test

% Don't vary hormones for K perturbations
% if isequal(type,'Ks')
    e2 = estr0;
    pr = 0;
    test = 0;
% end

sol_nominal = dde23(@ode_ftn, 14, @cytokine_hist_zero, tspan);

n_sets = 100;

sampling_matrix = lhsdesign(n_sets,num_coeffs);

sol_cell=cell(1,num_substances);
for i=1:num_substances
    sol_cell{i}=zeros(n_sets, length(tspan));
end

for i=1:n_sets
    % Vary Ks and number of SFs 
%     if isequal(type,'Ks')
        [Ks,P] = Ks_mod(sampling_matrix,i,0.6);
        Param(i,:)=P;
      
%     end
    
    
    sol = dde23(@ode_ftn, 14, @cytokine_hist_zero, tspan); % T = time, C = cytokine data
    sol1 = deval(sol,tspan);
    if min(sol1) < 0
        disp(min(sol1))
    end
    for j=1:num_substances
        sol_cell{j}(i,:)=sol1(j,:);
    end
  
%     [Gsen{i} O_r{i} O_v{i} O_b{i}]= Param_var_local_CWH(tspan,sol_nominal); % logarithmic sensativity analysis
    
end

% % -----------Plotting the critical trigger of chronic inflammation ---------*** CHANGE P_param TO CHECK FOR OTHER OUTPUTS ***
%   rnk = 3; % User defined number of most sensitizing mechanisms to plot; Currently the program plots the top three mechanisms as Rank 1, Rank 2 and Rank3 at every simulation time point
%   P_param = [8 9 12 13]; % The output identifying numbers for MMP-1, MMP-13 and MMP-3; For a complete list of all molecular mediator output identifying numbers check "cytokine_hist_CWH.m"
%   analysis_CWH(O_b,P_param,tspan,n_sets,rnk)
%   [h]=Sensativity_Analysis_Heatmap(Gsen,P_param);
% 
% median_sol_cell = cell(1,num_substances);
% low_sol_cell = cell(1,num_substances);
% high_sol_cell = cell(1,num_substances);

for i=1:num_substances
    median_sol_cell{i} = median(sol_cell{i},1);
    low_sol_cell{i} = prctile(sol_cell{i},25,1);
    high_sol_cell{i} = prctile(sol_cell{i},75,1);
end


median_sol = zeros(length(tspan),num_substances);
high_sol = zeros(length(tspan),num_substances);
low_sol = zeros(length(tspan),num_substances);

for i=1:length(median_sol_cell)
    median_sol(:,i)=median_sol_cell{i}';
    high_sol(:,i)=high_sol_cell{i}';
    low_sol(:,i)=low_sol_cell{i}';
end