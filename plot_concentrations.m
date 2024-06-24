function plot_concentrations(rows,cols,subs,sol,t_new,linetype)
% Plot line graphs for concentration profiles
% "hold on" is used after every plot so that multiple lines can be plotted
% together. 
% To plot multiple lines, add a new figure, call this function for the
% first scenario, then (without adding another figure), call the function
% again as many times as needed
% The subplot indices (rows, cols, and subs) are defined in "main_script"

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


T = t_new';
C = sol;
xmax=20;%max(T)/24;


% Platelet subplot
if subs(1) > 0
%     sp5=figure; 
% subplot(rows,cols,subs(1));
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,1).*1000,linetype,'LineWidth',4) 
    title('Platelet concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    axis square
    hold on
end

% M1 subplot
if subs(2) > 0
%     sp1=figure(1);
subplot(rows,cols,1);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,2),linetype,'LineWidth',4)
    title('M1 cell concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (cells/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    axis square
    hold on
end

% M2 subplot
if subs(3) > 0 
%     sp2=figure(2);
subplot(rows,cols,2);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,3),linetype,'LineWidth',4)
    title('M2 cell concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (cells/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    axis square
    hold on
end



% TNF-a subplot
if subs(5) > 0
%     sp4=figure(4);
subplot(rows,cols,3);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,5).*1000,linetype,'LineWidth',4)
    title('TNF-\alpha concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
     xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 12])
    axis square
    hold on
end

% IL-10 subplot
if subs(6) > 0
%     sp5=figure(5);
subplot(rows,cols,4);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,6).*1000,linetype,'LineWidth',4)
    title('IL-10 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 30])
    axis square
    hold on
end

% TGF-b subplot
if subs(7) > 0
%     sp6=figure(6);
subplot(rows,cols,5);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,7).*1000,linetype,'LineWidth',4)
    title('TGF-\beta concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    axis square
    hold on
end

% MMP-9 subplot
if subs(8) > 0
%     sp7=figure(7);
   subplot(rows,cols,6);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,8),linetype,'LineWidth',4)
    title('MMP-9 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (ng/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    axis square
    hold on
end

% MMP-1 subplot
if subs(9) > 0
%     sp8=figure(8);
    subplot(rows,cols,7);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,9),linetype,'LineWidth',4)
    title('MMP-1 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (ng/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
   %axis([0 xmax 0 3000])
    axis square
    hold on
end

% TIMP-1 subplot
if subs(10) > 0
%     sp9=figure(9);
    subplot(rows,cols,8);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,10),linetype,'LineWidth',4)
    title('TIMP-1 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (ng/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 4.5e5])
    axis square
    hold on
end

% IL-6 subplot
if subs(11) > 0
%     sp10=figure(10);
    subplot(rows,cols,9);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,11).*1000,linetype,'LineWidth',4)
    title('IL-6 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 25])
    axis square
    hold on
end

% IL-1b subplot
if subs(4) > 0
%     sp3 = figure(3);
    subplot(rows,cols,10);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,4).*1000,linetype,'LineWidth',4)
    title('IL-1\beta concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (pg/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 10])
    axis square
    hold on
end

% MMP-13 subplot
if subs(12) > 0
%     sp11=figure(11);
    subplot(rows,cols,11);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,12),linetype,'LineWidth',4)
    title('MMP-13 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (ng/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 3000])
    axis square
    hold on
    
%MMP-3 subplot
    if subs(13)>0
%        sp12=figure(12);
subplot(rows,cols,12);
%     axes('Position',[0.2 0.2 0.7 0.7])
    plot(T./24,C(:,13),linetype,'LineWidth',4)
    title('MMP-3 concentration','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    xlabel('Time (days)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    ylabel('Concentration (ng/mL)','FontSize',16,'FontWeight','Bold','Fontname','Arial')
    set(gca,'FontSize',16,'FontWeight','Bold','Fontname','Arial');
    %axis([0 xmax 0 3000])
    axis square
    hold on 
    end
    
end
