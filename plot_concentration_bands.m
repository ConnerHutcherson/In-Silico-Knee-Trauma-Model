function plot_concentration_bands(subs,median_sol,high_sol,low_sol,tspan,type)
% Plot confidence bands for concentration profiles
% "hold on" is used after every plot so that multiple lines can be plotted
% together. 
% To plot multiple lines, add a new figure, call this function for the
% first scenario, then (without adding another figure), call the function
% again as many times as needed
% The subplot indices (rows, cols, and subs) are defined in "main_script"

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu
% 

C = median_sol;
T = tspan';
xmax=20;%max(T)/24;

if isequal(type,'Ks')
    line_color='k';
elseif isequal(type,'Fem')
    line_color='r-';
elseif isequal(type,'FemPeak')
    line_color='g:';
elseif isequal(type,'Male')
    line_color='b--';
else
    line_color='k';
end

% if length(get(gca,'children'))>=4
%     line_color = 'm--';
% end
rows = 3;
cols = 4;
% M1 subplot
if subs(2) > 0
   sp1=figure(1);
%     subplot(rows,cols,1); %gca; %
    %axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,2)./1000,high_sol(:,2)./1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,2)./1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('M1 cell concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (*10^3 cells/mL)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 400])
    
    axis square
    hold on
end
 
% M2 subplot
if subs(3) > 0
     sp2=figure(2);
%     subplot(rows,cols,2);
%      axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,3)./1000,high_sol(:,3)./1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,3)./1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('M2 cell concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (*10^3 cells/mL)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%    axis([0 xmax 0 15])
    axis square
    hold on
end

% TNF-a subplot
if subs(5) > 0
     sp4=figure(4);
%     subplot(rows,cols,3);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,5).*1000,high_sol(:,5).*1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,5).*1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    title('TNFa','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('pg/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 12])
    axis square
    hold on
end

% IL-10 subplot
if subs(6) > 0
    sp5=figure(5);
     sp5=figure(5);
%     subplot(rows,cols,4);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,6).*1000,high_sol(:,6).*1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,6).*1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    title('IL-10','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('pg/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 40])
    axis square
    hold on
end

% IL-1b subplot
if subs(4) > 0
    sp3 = figure(3);
%     subplot(rows,cols,5);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,4).*1000,high_sol(:,4).*1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,4).*1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    title('IL-1b','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('pg/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 30])
    axis square
    hold on
end

% TGF-b subplot
if subs(7) > 0
    sp6=figure(6);
%     subplot(rows,cols,6);
%      axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,7),high_sol(:,7),tspan./24,line_color);
    hold on
    plot(T./24,C(:,7),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('TGFb','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 0.6])
    axis square
    hold on
end

% MMP-9 subplot
if subs(8) > 0
    sp7=figure(7);
%     subplot(rows,cols,7);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,8),high_sol(:,8),tspan./24,line_color);
    hold on
    plot(T./24,C(:,8),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('MMP-9 ','FontSize',16,'FontWeight','Bold','Fontname','Arial')
   xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 25])
    axis square
    hold on
end

% MMP-1 subplot
if subs(9) > 0
    sp8=figure(8);
%     subplot(rows,cols,8);
%      axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,9),high_sol(:,9),tspan./24,line_color);
    hold on
    plot(T./24,C(:,9),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('MMP-1 ','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 15])
    axis square
    hold on
end

% TIMP-1 subplot
if subs(10) > 0
    sp9=figure(9);
%     subplot(rows,cols,9);
%      axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,10),high_sol(:,10),tspan./24,line_color);
    hold on
    plot(T./24,C(:,10),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('TIMP-1','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 400])
    axis square
    hold on
end

% IL-6 subplot
if subs(11) > 0
    sp10=figure(10);
%     subplot(rows,cols,10);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,11).*1000,high_sol(:,11).*1000,tspan./24,line_color);
    hold on
    plot(T./24,C(:,11).*1000,line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    title('IL-6 ','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('pg/ml','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 800])
    axis square
    hold on
end

% MMP-13 subplot
if subs(12) > 0
    sp11=figure(11);
%     subplot(rows,cols,11);
%      axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,12),high_sol(:,12),tspan./24,line_color);
    hold on
    plot(T./24,C(:,12),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('MMP-13 ','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 3])
    axis square
    hold on
end

%MMP-3 subplot
if subs(13) > 0;
sp12=figure(12);
%     subplot(rows,cols,12);
%     axes('Position',[0.2 0.2 0.7 0.7])
    ciplot_mod(low_sol(:,13),high_sol(:,13),tspan./24,line_color);
    hold on
    plot(T./24,C(:,13),line_color,'LineWidth',2);%,'Color',[0.5,0.5,0.5])% 
    axis tight
    title('MMP-3','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Post-Injury(days)','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    ylabel('ng/mL','FontSize',10,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
%     axis([0 xmax 0 60])
    axis square
    hold on
end