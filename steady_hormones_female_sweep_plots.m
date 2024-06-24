subs = zeros(1,num_substances);
n = [2:13];
subs(n) = 1;

for k =1:cycle_length
    
    plot_concentrations_sweep_CWH(rows,cols,subs,sol_female_sweep_cell{k}',tspan);
%     plot_concentrations_sweep_CWH(rows,cols,subs,median_sol_female_cell{k}',tspan);
    hold on
    
end

