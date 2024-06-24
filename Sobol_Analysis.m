function [Sii,STii] = Sobol_Analysis(YA,YB,YC,num_coeffs,num_fcoefs,num_substances,tspan,warmup)

% Description: Function that facilitates the calculation of Sobol
% sensitivities for the model inputs by the SAFE VBSA toolbox. (Pianosi F
% et al 2015)

% Input:  Model output values for parameters effects and combined parameter
% effects (YA, YB, YC)

% number of model coefficients (num_coeffs)
% number of model substances (num_substances)
% model time vector (tspan)
% delay time for sensativity analysis (warmup) 

% Output: Sensativity indicy output for first and total order terms  (Sii,
% STii) 

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

T = size(YA{1},2) ;

 for i=2:num_substances
     
    Sii{i}=zeros(num_coeffs+num_fcoefs, length(tspan));
    STii{i}=zeros(num_coeffs+num_fcoefs, length(tspan));

%
 end

    for k=2:num_substances
   Si  = nan(T,num_coeffs+num_fcoefs) ;
   STi = nan(T,num_coeffs+num_fcoefs) ;


    for t=warmup:T
        [ Si(t,:), STi(t,:)] = vbsa_indices(YA{:,k}(:,t),YB{:,k}(:,t),YC{:,k}(:,t) );
          
          Sii{k} =Si;
          STii{k} = STi;
         

    end  
    end 
end

