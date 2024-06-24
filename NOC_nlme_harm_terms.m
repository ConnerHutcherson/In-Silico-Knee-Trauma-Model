function [xrange,Po_E2,Po_P4,Po_T] = NOC_nlme_harm_terms()

% ID               - subject number
% Subject_Code     - subject code
% age_group        - group based on age 20-25,25-30,30-35, and 35-40
% cyc_length_group - group based on the length of menstrual cycle
%                    see Mumford paper
%                      1 = cycle length <26 days, short
%                      2 = cycle length 26 - 35 days, normal
%                      3 = cycle length > 35 days, long
% Norm_cyc_day     - normalized and centered cycle day.
%                     0 = first day of cycle
%                     1 = last day of cycle
%                   0.5 = day when E2 reach its peak.
% E2               - estradiol contrations in pg/mL
% P4               - progesterone contrations in ng/mL
% Normalized_E2    - estradiol contrations normalized to is max (%)
% Normalized_P4    - progesterone contrations normalized to is max (%)
% HRT              - half relaxation time in ms
% normalized_HRT   - half relaxation time normalized to HRT at peak E2

% This script prepared by Subaryani Soedirdjo
% UT Southwestern Medical Center, 2023
% mail to: subaryani.soedirdjo@utsouthwestern.edu


% Hormonal Data
file_path = 'C:\Users\chutc\OneDrive\Documents\MATLAB\Knee_Trauma_Model_Prototype\'; % change to current file path

[~, ~,raw] = xlsread([file_path 'Hormone_profile_NOC_v2'],'Sheet1');


data = cell2mat(raw(2:end,[1,7,8,9,10]));


ID = data(:,1);
Norm_cyc_day = data(:,2);
E2 = data(:,3);
P4 = data(:,4);
T = data(:,5);


[~,~,stats_E2,~,~,xrange,Po_E2] = mixed_effect_hormones_10harms(ID(~isnan(E2)),Norm_cyc_day(~isnan(E2)),E2(~isnan(E2)));
[~,~,stats_P4,~,~,~,Po_P4] = mixed_effect_hormones_5harms(ID(~isnan(P4)),Norm_cyc_day(~isnan(P4)),P4(~isnan(P4)));
[~,~,stats_T,~,~,~,Po_T] = mixed_effect_hormones_4harms(ID(~isnan(T)),Norm_cyc_day(~isnan(T)),T(~isnan(T)));


end


function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_3harms(group,time,response)

options   = statset( 'MaxIter',10000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 3 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, .01*ones(1,7),...
                                'FEParamsSelect', [1 3 4 6 7 8 9],...
                                'REParamsSelect', [1 2 5], 'Options', options);                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1)];
    
    
xrange = 0:.01:1;
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end

end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_2harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 3 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, .01*ones(1,5),...
                                'FEParamsSelect', [1 3 4 6 7],...
                                'REParamsSelect', [1 2 5], 'Options', options);                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)];
    
    
xrange = 0:.01:1;
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end

end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_4harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, .01*ones(1,9),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1)                            ];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_5harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,11),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_6harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,13),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_7harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 16)  * cos(14 * pi* (t + PHI(:, 17) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,15),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15 16 17],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1),...
        beta(14)*ones(kk,1)        , beta(15)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_8harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 16)  * cos(14 * pi* (t + PHI(:, 17) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 18)  * cos(16 * pi* (t + PHI(:, 19) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,17),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1),...
        beta(14)*ones(kk,1)        , beta(15)*ones(kk,1), beta(16)*ones(kk,1),...
        beta(17)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_9harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 16)  * cos(14 * pi* (t + PHI(:, 17) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 18)  * cos(16 * pi* (t + PHI(:, 19) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 20)  * cos(18 * pi* (t + PHI(:, 21) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,19),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1),...
        beta(14)*ones(kk,1)        , beta(15)*ones(kk,1), beta(16)*ones(kk,1),...
        beta(17)*ones(kk,1)        , beta(18)*ones(kk,1), beta(19)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_10harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 16)  * cos(14 * pi* (t + PHI(:, 17) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 18)  * cos(16 * pi* (t + PHI(:, 19) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 20)  * cos(18 * pi* (t + PHI(:, 21) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 22)  * cos(20 * pi* (t + PHI(:, 23) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,21),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1),...
        beta(14)*ones(kk,1)        , beta(15)*ones(kk,1), beta(16)*ones(kk,1),...
        beta(17)*ones(kk,1)        , beta(18)*ones(kk,1), beta(19)*ones(kk,1),...
        beta(20)*ones(kk,1)        , beta(21)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_1harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 3 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, 0.1.*ones(1,3),...
                                'FEParamsSelect', [1 3 4],...
                                'REParamsSelect', [1 2 5], 'Options', options);                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'];
    
    
xrange = 0:.01:1;
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end

end

function [beta,PSI,stats,B,PHI,xrange,Po] = mixed_effect_hormones_11harms(group,time,response)

options   = statset( 'MaxIter',100000,'RobustWgtFun','logistic','TolFun',1e-10,'TolX',1e-10);
options.Display = 'iter';

% 4 harmonics mixed effect  model
nlme_mixedmodel = @(PHI, t)...
                  (PHI(:, 1) + exp(PHI(:, 2)) .*...
                   (PHI(:, 3)   * cos(2 * pi* (t + PHI(:, 4)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 6)   * cos(4 * pi* (t + PHI(:, 7)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 8)   * cos(6 * pi* (t + PHI(:, 9)   + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 10)  * cos(8 * pi* (t + PHI(:, 11)  + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 12)  * cos(10 * pi* (t + PHI(:, 13) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 14)  * cos(12 * pi* (t + PHI(:, 15) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 16)  * cos(14 * pi* (t + PHI(:, 17) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 18)  * cos(16 * pi* (t + PHI(:, 19) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 20)  * cos(18 * pi* (t + PHI(:, 21) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 22)  * cos(20 * pi* (t + PHI(:, 23) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))  +  ...
                    PHI(:, 24)  * cos(22 * pi* (t + PHI(:, 25) + exp(PHI(:,5))./(1 + exp(PHI(:,5)))))));
               
[beta, PSI, stats, B] = nlmefit(time, response, group, [], ...
                                nlme_mixedmodel, ones(1,23),...
                                'FEParamsSelect', [1 3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25],...
                                'REParamsSelect', [1 2 5], 'Options', options);
                            

% harmonics fluctutation is fixed effect. a nominal mean across subjects is expected - fixed
% mean, amplitude and phase-shift are random effect - subject dependent
% % fixed effects ( 5 subjects ) 
% beta
% % random effects (3 parameters x 11 subjects) 
% B
% % plotting the fixed effect model : 

kk = max(size(B(1,:)'));
PHI = [ beta(1)*ones(kk,1)+B(1,:)' , B(2,:)'            , beta(2)*ones(kk,1), ...
        beta(3)*ones(kk,1)         , B(3,:)'            , beta(4)*ones(kk,1), ...
        beta(5)*ones(kk,1)         , beta(6)*ones(kk,1) , beta(7)*ones(kk,1), ...
        beta(8)*ones(kk,1)         , beta(9)*ones(kk,1) , beta(10)*ones(kk,1),...
        beta(11)*ones(kk,1)        , beta(12)*ones(kk,1), beta(13)*ones(kk,1),...
        beta(14)*ones(kk,1)        , beta(15)*ones(kk,1), beta(16)*ones(kk,1),...
        beta(17)*ones(kk,1)        , beta(18)*ones(kk,1), beta(19)*ones(kk,1),...
        beta(20)*ones(kk,1)        , beta(21)*ones(kk,1), beta(22)*ones(kk,1),...
        beta(23)*ones(kk,1)];

xrange = min(time):.01:max(time);
Po = [];
for i = 1 : size(PHI,1)
   Po =  [Po; nlme_mixedmodel(PHI(i,:),xrange)];
end


end