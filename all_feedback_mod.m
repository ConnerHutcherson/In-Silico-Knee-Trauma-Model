function [f] = all_feedback_mod(f_nom)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Ksf

f_min = f_nom-f_nom.*0.2;
f_max = f_nom+f_nom.*0.2;

for j=1:length(f_nom)
    f(j) = Ksf(:,j)*(f_max(j)-f_min(j))+f_min(j);
end
end