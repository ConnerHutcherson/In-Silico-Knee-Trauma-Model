function [high_sol,low_sol,median_sol,sol_cell,Param_vary,e2_vary_corrected,pr_vary_corrected,test_vary_corrected,abcs_Param] = lhs_script_ks_and_et_prototype(num_substances,tspan,num_coeffs,num_fcoefs,E2_min_conc,E2_max_conc,E2_median_conc,P4_min_conc,P4_max_conc,P4_median_conc,T_min_conc,T_max_conc,T_median_conc,pct_vary,t_injury,t_hormone)
% This script performs latin hypercube sampling for the rate coefficients and hormone concentrations. 
% It calls "Ks_mod" to recalculate parameter sets with the sampling matrix
% It stores n_sets of solutions, calculates the median and 1st and 3rd
% quartiles

% Use with "plot_confidence_bands" 

% Output: model outputs for lower sample range (low_sol, 25% quartile), higher range
% (high_sol, 75% quartile) and median (median_sol)

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


global Ks e2 pr test initial_conditions Ksf

initial_conditions = zeros(15,1);
% 
input_e2 = (E2_max_conc-E2_min_conc)/2+E2_min_conc;
input_pr = (P4_max_conc-P4_min_conc)/2+P4_min_conc;
input_test = (T_max_conc-T_min_conc)/2+T_min_conc;

%Estrogen Input
e2=zeros(length(t_hormone),2);
e2(:,1) = t_hormone;
e2(:,2) = input_e2;
    
%Progesterone Input
pr=zeros(length(t_hormone),2);
pr(:,1) = t_hormone;
pr(:,2) = input_pr;

%Testosterone Input
test=zeros(length(t_hormone),2);
test(:,1) = t_hormone;
test(:,2) = input_test;

%Solve initial condition for inputted start point
sol_female_cycle_intial = dde23(@ode_ftn, 12, @cytokine_hist_zero, tspan);
vec_initial = size(sol_female_cycle_intial.x);
initial_conditions = sol_female_cycle_intial.y(:,(vec_initial(:,2)));
 
%% LHS Modeling
n_sets = 10;

% Nominal model run
% sol_nominal = dde23(@ode_ftn, 13, @cytokine_hist, tspan);

% Generate LHS sampling matrix
sampling_matrix = lhsdesign(n_sets,num_coeffs+num_fcoefs+3);

% sampling_matrix = lhsdesign(n_sets,num_coeffs+61);


% Generate Cell to store model results
sol_cell=cell(1,num_substances);

for i=1:num_substances
    sol_cell{i}=zeros(n_sets, length(tspan));
end
 
 e2_vary = zeros(n_sets,101);
 e2_vary_corrected = zeros(n_sets,101);
 pr_vary = zeros(n_sets,101);
 pr_vary_corrected = zeros(n_sets,101);
 test_vary = zeros(n_sets,101);
 test_vary_corrected = zeros(n_sets,101);
 Param_vary = zeros(n_sets,num_coeffs);
 abcs_Param = zeros(n_sets,num_fcoefs);

 
for i=1:n_sets
    % Vary Ks
    [Ks,P] = Ks_mod(sampling_matrix(:,1:num_coeffs),i,pct_vary);
    Param_vary(i,:)=P;

    % Feedback Parameter Modification
    [abcs] = abcs_mod(sampling_matrix(:,num_coeffs:num_coeffs+num_fcoefs),i,pct_vary);
    abcs_Param(i,:) = abcs;
    Ksf = abcs;
%     
% Ksf = sampling_matrix(i,num_coeffs+1:num_coeffs+58);
    
         
    % Vary hormone profiles
    e2=zeros(length(t_hormone),2);
    e2(:,1) = t_hormone;
    e2mod = e2_mod(sampling_matrix(:,num_coeffs+num_fcoefs+1),i,E2_min_conc, E2_max_conc, E2_median_conc, t_hormone);
%     e2mod = e2_mod(sampling_matrix(:,num_coeffs+59),i,E2_min_conc, E2_max_conc, E2_median_conc, t_hormone);
    e2mod_corrected = max(e2mod,0); % Makes negative  values equal 0
    e2(:,2) = e2mod_corrected';
    [e2,~] = change_injury_timing(e2,t_injury);

    e2_vary(i,:) = e2mod;
    e2_vary_corrected (i,:) = e2mod_corrected;
    
    pr=zeros(length(t_hormone),2);
    pr(:,1) = t_hormone;
    prmod = pr_mod(sampling_matrix(:,num_coeffs+num_fcoefs+2),i,P4_min_conc,P4_max_conc,P4_median_conc, t_hormone);
%     prmod = pr_mod(sampling_matrix(:,num_coeffs+60),i,P4_min_conc,P4_max_conc,P4_median_conc, t_hormone);
    prmod_corrected = max(prmod,0); % Makes negative  values equal 0
    pr(:,2)= prmod_corrected';
    [pr,~] = change_injury_timing(pr,t_injury);
    
    pr_vary(i,:) = prmod;
    pr_vary_corrected (i,:) = prmod_corrected;
    
    
    test=zeros(length(t_hormone),2);
    test(:,1) = t_hormone;
    tmod = t_mod(sampling_matrix(:,num_coeffs+num_fcoefs+3), i ,T_min_conc,T_max_conc,T_median_conc, t_hormone);
%   tmod = t_mod(sampling_matrix(:,num_coeffs+61), i ,T_min_conc,T_max_conc,T_median_conc, t_hormone);
    tmod_corrected = max(tmod,0); % Makes negative  values equal 0
    test(:,2) = tmod_corrected';
    [test,~] = change_injury_timing(test,t_injury);
    
    test_vary(i,:) = tmod;
    test_vary_corrected (i,:) = tmod_corrected;
 
    % Run Trauma Model
    sol = dde23(@ode_ftn_Sobol, 13, @cytokine_hist, tspan); % T = time, C = cytokine data
    x = deval(sol,tspan);
    
    if min(x) < 0
        disp(min(x))
    end
    
    % Store Model Results in Cell
    for j=1:num_substances
        sol_cell{j}(i,:)=x(j,:);
    end
end

% Summarize model results
median_sol_cell = cell(1,num_substances);
low_sol_cell = cell(1,num_substances);
high_sol_cell = cell(1,num_substances);

% Store summarized model results
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