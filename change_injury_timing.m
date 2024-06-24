function [rearranged_h,cycle_t] = change_injury_timing(h,t_injury)
% The purpose of this code is to restructure the estrogen array so that it
% starts at different times in the menstrual cycle

% Inputs: 
% h: hormone array with a column that has ascending time values in the first
% column and corresponding hormone (e2 or p) values in the second column
% (current format in the model code)
% t_injury: normalized time in the menstrual cycle where injury should
% occur (takes on a value between 0 and 1)
% cycle_length: the number of days in a simulated cycle (e.g., 29)


% Outputs: 
% rearranged_h: array with time values for the ode solver (time still 
% starts at 0, for the injury simulation, even when a person is on cycle
% day n) and the rearranged hormone values (e2 or p). In other words the
% first column of the h array corresponds to time after injury.
% cycle_t: time array that correpsonds to the cycle days (starts on cycle
% day n). Might not be necessary

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

% Convert normalized t_injury to the cycle day at injury
t_injury_abs = t_injury*24;


% Use array logic to identify elements of the e2 array where time is
% greater than or equal to t_injury_abs
inds = h(:,1)>=t_injury_abs;


% Use the indices to pull the values of the hormone and time from times at and  after
% t_injury_abs
h_after = h(inds,2);
t_after = h(inds,1);


% Use the indices to pull the values of hormone and time before t_injury_abs
h_before = h(~inds,2);
t_before = h(~inds,1);

% Combine the arrays
rearranged_h = zeros(size(h));
rearranged_h(:,1) = h(:,1); % The time column remains the same as the input (time since injury, which the ode solver needs)
rearranged_h(:,2) = [h_after; h_before]; % The hormone column gets rearranged so that the hormone value at the time of the injury is the value it would take on the specified day of the cycle

% Might not be necessary - I don't think this time array will work in the
% ode solver. Gives the day of the menstrual cycle that corresponds to the
% current hormone value. 
cycle_t = [t_after; t_before];

