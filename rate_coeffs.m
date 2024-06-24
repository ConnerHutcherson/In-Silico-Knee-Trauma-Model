function Ks = rate_coeffs()

% Rate coefficients for the system of differential equations
% These are called in "ode_ftn" 

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu


% Platelets
Ks(1) = 0.69; % hr^-1 Platelet degradation, Nagaraja SI ref. 13

% Macrophages
Ks(2) = 400; % mL^-1 M1 migration (assumed by Nagaraja et al. 2014 see 2-2-17 notes for comments)
Ks(3) = 8.3e-3; % hr^-1 Macrophage lymphatic removal (M1 & M2 Assumption), Nagaraja SI ref. 3
Ks(4) = 1/12;  % hr^-1 Macrophage phenotype transition parameter (M1 to M2) (ask Bethany why 1/12 rather 0.1 in Nagaraja SI ref. 3)

% IL-1b
Ks(5) = 1.23e-6; % ng/(cell*hr) M1 production of IL-1b, Nagaraja SI ref. 7
Ks(6) = 2.45e-7; % ng/(cell*hr) M2 production of IL-1b, Nagaraja SI ref. 7
Ks(7) = 1.028e-9; % ng/(cell*hr) SF production of IL-1b, see 2-24-17 notes
Ks(8)=  5.33e-11; % ng/(cell*hr) Chondrocyte production of IL-1b, Kaneva et al. 2012 
Ks(9) = 0.693/(3/60); %  hr^-1 IL-1b degradation, Nagaraja SI ref. 8 


% TNF-a
Ks(10) = 3.46e-7; % ng/(cell*hr) M1 production of TNF-a, Nagaraja SI ref. 7    
Ks(11) = 4.29e-8; % ng/(cell*hr) M2 production of TNF-a, Nagaraja SI ref. 7
Ks(12) = 2.583e-9; % ng/(cell*hr) SF production of TNF-a Huang et al. 2011 see 2-24-2017 notes
Ks(13)= 4.045e-10; % ng/(cell*hr) Chondrocyte production of TNF-a,(Parker et al. Arthritis Research & Therapy 2013)
Ks(14) = 0.693/(5/60);% hr^-1 TNF-a degradation, Kaneda et al. 2004 reports 5 minute half life of native TNFa 

% IL-10
Ks(15) = 7.60e-8; % ng/(cell*hr) M1 production of IL-10, Nagaraja SI ref. 7
Ks(16) = 1.55e-7; % ng/(cell*hr) M2 production of IL-10, Nagaraja SI ref. 7
Ks(17) = 3.0556e-10; % ng/(cell*hr) SF production of IL-10 Huang et al. 2011 see 2-24-17 notes
Ks(18)= 1.52e-9; % ng/(cell*hr) Chondrocyte production of IL-10,R.D. Müller et al. 2008 
Ks(19) = 0.693/(20/60); % hr^-1 IL-10 degradation, Nagaraja SI ref. 1

% TGF-b
Ks(20) = 7.447e-10; % ng/(cell*hr) Chondrocyte production of TGF-b,. Homandberg et al. 1997
Ks(21) = 1.88e-6; % ng/(cell*hr) M1 production of TGF-b, Nagaraja SI ref. 4
Ks(22) = 1.6e-8; % ng/(cell*hr) M2 production of TGF-b, Nagaraja SI ref. 4
Ks(23) = 3.15e-8; % ng/(cell*hr) SF production of TGF-b, Li et al. 2011 see 2-24-2017 notes
Ks(24) = 0.693/(15/60); % hr^-1 TGF-b degradation, Tarant 2010: half-life of TGF is generally < 1 hr in blood
Ks(46) = 1.25e-8; % ng/(cell*hr) Platelet production of TGF-b, Nagaraja SI ref. 14, 15

% MMP-9
Ks(25)=6.7708e-6; % ng/(cell*hr) M1 production of MMP-9, From Jager et al. 2016
Ks(26)=1.4323e-5; % ng/(cell*hr) M2 production of MMP-9, From Jager et al. 2016
% Ks(27)= 6.36667e-07; % ng/(cell*hr) Chondrocyte production of MMP-9, Meszaros et al. 2015 
t_half_MMP9 = 3.26; % hr, MMP-9 clearance from joint (assumed) T Mwang et al. 2018
Ks(28) = 0.693/t_half_MMP9; % hr^-1

% MMP-1
Ks(29) = 2e-8;% ng/(cell*hr) M1 production of MMP-1 Serra et al. 2010 (Nagaraja 2017)
Ks(30) = 7.24e-8; % ng/(cell*hr) SF production of MMP-1 Cha et al. 2003
Ks(31)= 4.6333e-06; % ng/(cell*hr) Chondrocyte production of MMP-1, (RANKIN et al. 2020)
t_half_MMP1 = 3.26 ; % hr MMP-1 clearance from joint (assumed) T Mwang et al. 2018 
Ks(32) = 0.693/t_half_MMP1;

% TIMP-1
Ks(33) = 3.89e-7; % ng/(cell*hr) M1 production of TIMP-1, Russell et al. 2002 fig 3b, no stim.
Ks(34) = 2.594e-7; % ng/(cell*hr) M2 production of TIMP-1, Russell et al. 2002 fig 3b, max stim. 
Ks(35) = 3.05e-4; % ng/(cell*hr) SF production of TIMP-1, Asano et al. 2006
Ks(36)= 1.813e-5; % ng/(cell*hr) Chondrocyte production of TIMP-1, Fearon et al. 2006 
t_half_TIMP1 = 1.1; % hr TIMP-1 half life Batra et al. 2012
Ks(37) = 0.693/t_half_TIMP1; % hr^-1 TIMP-1 degradation


% IL-6
Ks(38) = 1.18e-6; % 5e-5; % ng/(cell*hr) M1 production of IL-6 Nagaraja et al. SI ref 10
Ks(39) = Ks(38)*0.1; % ng/(cell*hr) M2 production of IL-6 Nagaraja et al. assumed
Ks(40)= 5.08e-10; % ng/(cell*hr) Chondrocyte production of IL-6,Kaneva et al. 2012 
t_half_IL6 = 1 ; % hr IL-6 half life, May LT, et al 1995
Ks(41) =  0.693/t_half_IL6; % hr^-1 IL-6 degradation Nagaraja et al. SI ref 8
Ks(42) = 65/(5e3/0.1*48)/1000; % ng/(cell*hr) SF production of IL-6 Inuoe et al. 2001

% MMP-13
Ks(43)= 1.1521e-7; %ng/(cell*hr) Chondrocyte production of MMP-13 Hashizume et al. 2010)
Ks(44)= 1.16667e-06; %ng/(cell*hr) SF production of MMP-13, Asano et al. 2006
t_half_MMP13 = 3.26; %hr MMP-13 clearance from joint (assumed) T Mwang et al. 2018
Ks(45)= 0.693/ t_half_MMP13; % hr^1 MMP-13, 

% MMP-3
Ks(47)= 3.5000e-05; %ng/(cell*hr) Chondrocyte production of MMP-3, (Wu et al. 2014)
Ks(48)=4.89583e-07; %ng/(cell*hr) SF production of MMP-3,Fuchs, S. et al. 2004
t_half_MMP13 = 3.26 ; %hr MMP-3 clearance from joint (assumed) T Mwang et al. 2018
Ks(49)= 0.693/ t_half_MMP13; % hr^1 MMP-3, 


