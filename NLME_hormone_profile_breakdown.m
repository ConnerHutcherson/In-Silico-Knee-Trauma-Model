function [P4_min_conc,P4_max_conc,P4_mean_conc,P4_median_conc,P4_bottom10,P4_top10,E2_min_conc,E2_max_conc,E2_mean_conc,E2_median_conc,E2_bottom10,E2_top10,T_min_conc,T_max_conc,T_mean_conc,T_median_conc] = NLME_hormone_profile_breakdown(t_hormone,E2_profiles,P4_profiles,T_profiles)
% The purpose of this function is to summarize the female sex hormone profiles generated
% by the NLME model

% Inputs:
% t_hormone: hormone profile time values (hr)
% E2_profiles: Corresponding Estrogen concentration profiles fitted in the NLME Model (pg/ml)
% P4_profiles: Corresponding Progesterone concentration profiles fitted in the NLME Model (ng/ml)
% T_profiles: Corresponding Testosterone concentration profiles fitted in the NLME Model (ng/dl)


% Outputs:
% P4_min_conc: minimum menstrual Progesterone concentration values determined across female human population (uM/L)
% P4_max_conc: maximum menstrual Progesterone concentration values determined across female human population (uM/L)
% P4_mean_conc: average menstrual Progesterone concentration values determined across female human population (uM/L)
% P4_bottom10: lower bound on the top 10% of Progesterone concentration values determined across female human population (uM/L)
% P4_top10: upper bound on the top 10% of Progesterone concentration values determined across female human population (uM/L)
% E2_min_conc: minimum menstrual Estradiol concentration values determined across female human population (pM/L)
% E2_max_conc: maximum menstrual Estradiol concentration values determined across female human population (pM/L)
% E2_mean_conc: average menstrual Estradiol concentration values determined across female human population (pM/L)
% E2_bottom10: lower bound on the top 10% of Estradiol concentration values determined across female human population (pM/L)
% E2_top10: upper bound on the top 10% of Estradiol concentration values determined across female human population (pM/L)
% T_min_conc: minimum menstrual Testosterone concentration values determined across female human population (nM/L)
% T_max_conc: maximum menstrual Testosterone concentration values determined across female human population (nM/L)
% T_mean_conc: average menstrual Testosterone concentration values determined across female human population (nM/L)

% This script prepared by Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

% Molecular weights for estrogen, progesterone and testosterone needed to
% convert to input values to nM/L units
e2_wt=272.38;% g/mol
p4_wt=314.46; % g/mol
t_wt=288.4; %g/mol

for q = 1:size(t_hormone,2)
    
% Male/post-menopausal estrogen serum concentration
P4_min_conc(q) = (min(P4_profiles(:,q))/p4_wt); % uM from nOC Females
P4_min_conc = max(P4_min_conc,0); % Makes negative P4 concentration values equal 0
P4_max_conc(q) = (max(P4_profiles(:,q))/p4_wt); % uM from nOC Females
P4_mean_conc(q) = (mean(P4_profiles(:,q))/p4_wt); %uM from nOC Females
P4_median_conc(q) = (median(P4_profiles(:,q))/p4_wt); %uM from nOC Females
P4_bottom10(q) = P4_min_conc(:,q)+(P4_max_conc(:,q)-P4_min_conc(:,q))*0.20; % Upper bound on the lowest 10% of progesterone concentrations
P4_top10(q) = P4_max_conc(:,q)-(P4_max_conc(:,q)-P4_min_conc(:,q))*0.20; % Lower bound on the top 10% of progesterone concentrations
E2_min_conc(q) = (min((E2_profiles(:,q))/1000)/e2_wt)*1e3; % nM from nOC Females
E2_min_conc = max(E2_min_conc,0); % Makes negative E2 concentration values equal 0
E2_max_conc(q) = (max((E2_profiles(:,q))/1000)/e2_wt)*1e3; % nM from nOC Females
E2_mean_conc(q) = (mean((E2_profiles(:,q))/1000)/e2_wt)*1e3; %nM from nOC Females
E2_median_conc(q) = (median((E2_profiles(:,q))/1000)/e2_wt)*1e3; %nM from nOC Females
E2_bottom10(q) = E2_min_conc(:,q)+(E2_max_conc(:,q)-E2_min_conc(:,q))*0.20; % Upper bound on the lowest 10% of estrogen concentrations
E2_top10(q) = E2_max_conc(:,q)-(E2_max_conc(:,q)-E2_min_conc(:,q))*0.20; % Lower bound on the top 10% of estrogen concentrations
T_min_conc(q) = (min((T_profiles(:,q))/(t_wt/10))); % nM from nOC Females
T_max_conc(q) = (max((T_profiles(:,q))/(t_wt/10))); % nM from nOC Females
T_mean_conc(q) = (mean((T_profiles(:,q))/(t_wt/10))); %nM from nOC Females
T_median_conc(q) = (median((T_profiles(:,q))/(t_wt/10))); %nM from nOC Females

end
end

