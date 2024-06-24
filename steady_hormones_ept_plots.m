subs = zeros(1,num_substances);
n = [2:13];
subs(n) = 1;
plot_concentrations_CWH(rows,cols,subs,sol_female_menses,sol_female_menses_t,'-r')
