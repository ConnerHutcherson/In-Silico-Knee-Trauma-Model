function h = cytokine_hist_zero(t)

% Function called prior to cytokine_hist that sets the initial conditions at 0 for the calculation of the initial conditions for the specified injury time (during the menstrual cycle).

% Input: time (t)

% Output: initiatl conditions for modeled substances (h)

% This script prepared by Bethany Luke and Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

% Specify values of the functions before and at t=0

 if t < 0
    h = zeros(15,1);%
elseif t==0
    h = zeros(15,1);%
end

     
% Uncomment to run for Initial Condition Values, %Comment actual values
h(1)=0; %Platelets
h(2)=0; %M1 Macrophage/ml
h(3)=0; %M2 Macrophage/ml
h(4)=0; %IL-1
h(5)=0; %TNF
h(6)=0; %IL-10
h(7)=0; %TGF
h(8)= 0; %MMP-9
h(9)=0; %MMP-1
h(10)=0;%TIMP
h(11)=0; %IL-6
h(12)= 0; % MMP-13
h(13)= 0; %MMP-3
h(14)=5e5; %SF initial value cells/ml
h(15)=5e5; %CH initial value cells/ml


