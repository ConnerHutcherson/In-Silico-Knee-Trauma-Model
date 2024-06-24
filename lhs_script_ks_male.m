function [high_sol_male,low_sol_male,median_sol_male,sol_cell_male, Param_vary, e2_vary_male, pr_vary_male, test_vary_male,abcs_Param] = lhs_script_ks_male(num_substances,tspan,num_coeffs,num_fcoefs,e_min,e_max,p_min,p_max,t_min,t_max,pct_vary,~,t_hormone)
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

%% Establish Initiation Conditions

h(1)=2.0e6; %Resident M0 density (cells/ml)
h(2)=0; % M1 initial value
h(3)=0; % M2 initiatl value
h(4)=initial_conditions(4,:); %IL-1
h(5)=initial_conditions(5,:);%TNF
h(6)=initial_conditions(6,:); %IL-10
h(7)=initial_conditions(7,:);%TGF
h(8)=initial_conditions(8,:); %MMP-9
h(9)=initial_conditions(9,:);%MMP-1
h(10)=initial_conditions(10,:);%TIMP
h(11)=initial_conditions(11,:); %IL-6
h(12)= initial_conditions(12,:); % MMP-13
h(13)=initial_conditions(13,:); %MMP-3
h(14)= 5e5; %SF initial value (cells/ml)
h(15)= 5e5; %CH initial value (cells/ml)
h(16)= 0;

% 
input_e2 = (e_max-e_min)/2+e_min;
input_pr = (p_max-p_min)/2+p_min;
input_test = (t_max-t_min)/2+t_min;

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

 
%% LHS Modeling
n_sets = 1000;

% Nominal model run
% sol_nominal = dde23(@ode_ftn, 13, @cytokine_hist, tspan);

% Generate LHS sampling matrix
sampling_matrix = lhsdesign(n_sets,num_coeffs+num_fcoefs+3);

% Generate Cell to store model results
sol_cell_male=cell(1,num_substances);

for i=1:num_substances
    sol_cell_male{i}=zeros(n_sets, length(tspan));
end

 e2_vary_male = zeros(n_sets,101);
 pr_vary_male = zeros(n_sets,101);
 test_vary_male = zeros(n_sets,101);
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

    % Vary hormone profiles
    e2=zeros(length(t_hormone),2);
    e2(:,1) = t_hormone;
    e2(:,2) = sampling_matrix(i,num_coeffs+num_fcoefs+1)*(e_max-e_min)+e_min; %sampling_matrix(i,1)
    
    e2_vary_male(i,:) = e2(:,2)';
    
    pr=zeros(length(t_hormone),2);
    pr(:,1) = t_hormone;
    pr(:,2) = sampling_matrix(i,num_coeffs+num_fcoefs+2)*(p_max-p_min)+p_min;
    
    pr_vary_male(i,:) = pr(:,2)';
    
    test=zeros(length(t_hormone),2);
    test(:,1) = t_hormone;
    test(:,2) = sampling_matrix(i,num_coeffs+num_fcoefs+3)*(t_max-t_min)+t_min; % sampling_matrix(i,2)
    
    test_vary_male(i,:) = test(:,2)';
    
    % Run Trauma Model
    sol = dde23(@ode_ftn_Sobol, 13, @cytokine_hist, tspan); % T = time, C = cytokine data
    x = deval(sol,tspan);
    
    if min(x) < 0
        disp(min(x))
    end
    
    % Store Model Results in Cell
    for j=1:num_substances
        sol_cell_male{j}(i,:)=x(j,:);
    end
end

% Summarize model results
median_sol_cell = cell(1,num_substances);
low_sol_cell = cell(1,num_substances);
high_sol_cell = cell(1,num_substances);

% Store summarized model results
for i=1:num_substances
    median_sol_cell{i} = median(sol_cell_male{i},1);
    low_sol_cell{i} = prctile(sol_cell_male{i},25,1);
    high_sol_cell{i} = prctile(sol_cell_male{i},75,1);
end

median_sol_male = zeros(length(tspan),num_substances);
high_sol_male = zeros(length(tspan),num_substances);
low_sol_male = zeros(length(tspan),num_substances);

for i=1:length(median_sol_cell)
    median_sol_male(:,i)=median_sol_cell{i}';
    high_sol_male(:,i)=high_sol_cell{i}';
    low_sol_male(:,i)=low_sol_cell{i}';
end