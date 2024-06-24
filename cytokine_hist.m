function h = cytokine_hist(t)

% Function that sets the initial conditions for the differential equations and specifies the trigger of the inflammatory process

% Input: time (t)

% Output: initiatl conditions for modeled substances (h)

% This script prepared by Bethany Luke and Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

% Specify values of the functions before and at t=0

global initial_conditions
  
 if t < 0
    h = zeros(15,1);%
elseif t==0
    h = zeros(15,1);%
    h(1) = 2e8;
end

     
% %ng/ml concentrations (normal run)
% h(4)=0.0400/1000;%IL-1
% h(5)=0.1796/1000;%TNF
% h(6)=0.439/1000;%IL-10
% h(7)=5.830/1000;%TGF
% h(8)= 0; %MMP-9
% h(9)=0.4216;%MMP-1
% h(10)=254.8;%TIMP
% h(11)=19.91/1000; %IL-6
% h(12)= 3.071; % MMP-13
% h(13) = 5.251; %MMP-3
% h(14)=5e5; %SF initial value cells/ml
% h(15)=5e5; %CH initial value cells/ml


% Uncomment to run for Initial Condition Values, %Comment actual values
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
h(14)=initial_conditions(14,:); %SF initial value cells/ml
h(15)=initial_conditions(15,:); %CH initial value cells/ml



