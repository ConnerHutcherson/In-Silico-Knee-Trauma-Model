function [Y,Sii,STii,median_sol,low_sol,high_sol] = lhs_script_sobol(num_substances,tspan,num_coeffs,num_fcoefs,e_min,e_max,p_min,p_max,t_min,t_max,range_factor,t_injury,t_hormone)
    
% This script will guide through the application of variance based
% sensativity analysis for the time-varying outputs of the Hutcherson et al 2023 trauma model. 
% The time-varying model output is the model prediction
% (estimated molecule concentration) at each time step of the simulation.

% REFERENCES

% Reusser, D. E., Zehe, E., 2011. Inferring model structural deficits 
% by analyzing temporal dynamics of model performance and parameter 
% sensitivity. Water Resources Research 47 (7).

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

tic

global Ks e2 estr0 test pr initial_conditions


%% Step 1: Set current directory to 'my_dir' and add path to sub-folders:    

my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2 Setup model and define input variability space

initial_conditions = zeros(15,1);


input_e2 = (e_max-e_min)/2+e_min;
input_pr = (p_max-p_min)/2+p_min;
input_test = (t_max-t_min)/2+t_min;

% Estrogen Input
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

% Function to rearrange hormone profiles and time vector   
[e2,~] = change_injury_timing(e2,t_injury);
[pr,~] = change_injury_timing(pr,t_injury);
[test,~] = change_injury_timing(test,t_injury);
% 
% %Solve initial condition for inputted start point
sol_female_cycle_intial = dde23(@ode_ftn_CWH, 12, @cytokine_hist_zero_CWH, tspan);
vec_initial = size(sol_female_cycle_intial.x);
initial_conditions = sol_female_cycle_intial.y(:,(vec_initial(:,2)));

%% Step 3 Define choices for sampling
n_sets = 4000;
warmup = 2;% lenght of model warmup period (days)
% (sensitivity indices will not be computed for the warmup period)

%% Step 4 Sample inputs space

%LHS Sampling Matrix
sampling_matrix = lhsdesign(n_sets,num_coeffs+3);


%Parameter Modification
for i=1:n_sets
    % Vary Ks and number of SFs 
    
    [Ks,P] = Ks_mod_CWH(sampling_matrix(:,1:num_coeffs),i,range_factor);    
    Param(i,:)= P;
    
end

%Hormone Modification
for k=1:n_sets
    e2=zeros(length(t_hormone),2);
    estrogen(:,1) = t_hormone;
    estrogen(:,2) = sampling_matrix(k,num_coeffs+1)*(e_max-e_min)+e_min; %sampling_matrix(i,1)
    [estrogen,~] = change_injury_timing(estrogen,t_injury);
    
    progesterone=zeros(length(t_hormone),2);
    progesterone(:,1) = t_hormone;
    progesterone(:,2) = sampling_matrix(k,num_coeffs+2)*(p_max-p_min)+p_min;
    [progesterone,~] = change_injury_timing(progesterone,t_injury);
    
    testosterone=zeros(length(t_hormone),2);
    testosterone(:,1) = t_hormone;
    testosterone(:,2) = sampling_matrix(k,num_coeffs+3)*(t_max-t_min)+t_min; % sampling_matrix(i,2)
    [testosterone,~] = change_injury_timing(testosterone,t_injury);

    E2(:,k)= estrogen(:,2);
    P4(:,k)= progesterone(:,2);
    T(:,k)= testosterone(:,2);
end

hormone_param = [Param,e2(2,:)',pr(2,:)',test(2,:)']; %combine parameter and hormone inputs

%Resample Parameter Space for Sobol Analysis
 [ XA, XB, XC ] = vbsa_resampling(Param) ;
%  [ E2A, E2B, E2C ] = vbsa_resampling(E2');
%  [ P4A, P4B, P4C ] = vbsa_resampling(P4');
%  [ TA, TB, TC ] = vbsa_resampling(T');
%  

%% Step 5 evaluate the model 

% Set up output matricies
 YA=cell(1,num_substances);
 YB=cell(1,num_substances);
 YC=cell(1,num_substances);
 
 for i=1:num_substances
    YA{i}=zeros(size(XA,1), length(tspan));
    YB{i}=zeros(size(XB,1), length(tspan));
    YC{i}=zeros(size(XC,1), length(tspan));
end

Ks = []; %clear global paramater array
e2 = [];
pr = [];
test = [];

% Evaluate model using input matricies XA and XB

for i=1:size(XA,1)
    
    Ks = XA(i,1:49);
        
 % Vary hormone profiles
    e2=zeros(length(t_hormone),2);
    e2(:,1) = t_hormone;
    e2(:,2) = E2A(i,:)';
    
    pr=zeros(length(t_hormone),2);
    pr(:,1) = t_hormone;
    pr(:,2) = P4A(i,:)';

    
    test=zeros(length(t_hormone),2);
    test(:,1) = t_hormone;
    test(:,2) = TA(i,:)';

    Ya = dde23(@ode_ftn_CWH, 13, @cytokine_hist_CWH, tspan) ; % size (N,T)
    Yaa = deval(Ya,tspan);
        for j=1:num_substances
        YA{j}(i,:)=Yaa(j,:);
    
        end

    Ks = XB(i,:);
    e2(:,2) = E2B(i,:)';
    pr(:,2) = P4B(i,:)';
    test(:,2) = TB(i,:)';

%    

    Yb = dde23(@ode_ftn_CWH, 13, @cytokine_hist_CWH, tspan) ; % size (N,T)
    Ybb = deval(Yb,tspan);
         for j=1:num_substances
         YB{j}(i,:)=Ybb(j,:);
    
         end

end

% Generate Median and IQR of nominal model outputs

median_sol_cell = cell(1,num_substances);
low_sol_cell = cell(1,num_substances);
high_sol_cell = cell(1,num_substances);

for i=1:num_substances
    median_sol_cell{i} = median(YA{i},1);
    low_sol_cell{i} = prctile(YA{i},25,1);
    high_sol_cell{i} = prctile(YA{i},75,1);
end


median_sol = zeros(length(tspan),num_substances);
high_sol = zeros(length(tspan),num_substances);
low_sol = zeros(length(tspan),num_substances);

for i=1:length(median_sol_cell)
    median_sol(:,i)=median_sol_cell{i}';
    high_sol(:,i)=high_sol_cell{i}';
    low_sol(:,i)=low_sol_cell{i}';
end

%% Evaluate model using recomination matrix XC 

Ks=[]; %clear global paramater array
e2 = [];
pr = [];
test = [];

 for k=1:size(XC,1)
    
    Ks = XC(k,:);
      
 % Vary hormone profiles
    e2=zeros(length(t_hormone),2);
    e2(:,1) = t_hormone;
    e2(:,2) = E2C(k,:)';
 
    pr=zeros(length(t_hormone),2);
    pr(:,1) = t_hormone;
    pr(:,2) = P4C(k,:)';
     
    test=zeros(length(t_hormone),2);
    test(:,1) = t_hormone;
    test(:,2) = TC(k,:)';     
    
    Yc = dde23(@ode_ftn_CWH, 12, @cytokine_hist_CWH, tspan); % size (N*M,T)
    Ycc = deval(Yc,tspan);
   
         for j=1:num_substances
         YC{j}(k,:)=Ycc(j,:);
    
         end        
end

%% Store all model outputs in a single cell

Y = [ YA; YB; YC ] ; % (N*(2+M),T) ... only needed for the next plot!


%% Step 6 Compute time-varying Sensitivity Indices (Sobol)

[Sii,STii] = Sobol_Analysis(YA,YB,YC,num_coeffs,num_substances,tspan,warmup);

%Select sensitivity index to be plotted in the figure:
Si_plot = 1;
STi_plot = 0;

%plot sensativity indicies
Sobol_Plots(Sii,STii,Si_plot,STi_plot,num_coeffs,num_substances)
end