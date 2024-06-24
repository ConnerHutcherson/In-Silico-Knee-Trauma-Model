
function a = all_parameter_fits()

% Values for all feedback coefficients. Each coefficient was determined
% with the experimental data directly above it and the curve fitting tool.
% These coeffiencts are called in "all_feedback_functions.m"


%% Synovial fibroblast feedback parameter fits
% TNF-alpha regulation of MMP-1 (from Asano et al. 2006 table 2)
TNF_fit = [0,20]; % ng/mL
MMP1_TNF = [25.7,93.8]; % ng/mL
MMP1_TNF = (MMP1_TNF-MMP1_TNF(1))/MMP1_TNF(1); % Normalize
a(1) = 2.782; 

% TNF-alpha regulation of TIMP-1 (from Asano et al. 2006 table 3)
TIMP_TNF = [146.8,228.7];
TIMP_TNF = (TIMP_TNF-TIMP_TNF(1))/TIMP_TNF(1);
a(2) = 0.5858; 

% TNF-A regulation of MMP-13 (from Asano et al. 2006 table 1 & 2)
MMP13_TNF = [0.56,1.6]; % ng/mL
MMP13_TNF = (MMP13_TNF-MMP13_TNF(1))/MMP13_TNF(1); 
a(3)= 1.950;

% TNF-alpha regulation of IL1b (Ganesan et al. 2012 table 3)
TNF_fit = [0, 1.6]; % ng/mL
IL1_TNF = [62.4,180.6]; % pg/mL
IL1_TNF = (IL1_TNF-IL1_TNF(1))/IL1_TNF(1); % Normalize
a(4) = 3.078; 

% TNF-alpha regulation of TGF-beta (Li et al. 2011 table 1)
TNF_fit = [0, 20]; % ng/mL
TGF_TNF = [409.7, 1305.4]; % pg/mL
TGF_TNF = (TGF_TNF-TGF_TNF(1))/TGF_TNF(1);
a(5) = 2.296; 

% IL-1b regulation of TIMP1 (Huang et al. 2010)
IL1_fit = [0, 1]; % ng/mL
TIMP1_IL1 = [46.5, 52]; % ng/mL
TIMP1_IL1 = (TIMP1_IL1-TIMP1_IL1(1))/TIMP1_IL1(1);
a(6) = 0.2366;

% IL-1b regulation of IL-1b (Huang et al. 2010)
IL1_IL1 = [3.7,22.4];
IL1_IL1 = (IL1_IL1-IL1_IL1(1))/IL1_IL1(1);
a(7) = 10.11; 

% IL-1b regulation of TNF-a (Huang et al. 2010)
TNF_IL1 = [9.3, 22.8];
TNF_IL1 = (TNF_IL1-TNF_IL1(1))/TNF_IL1(1);
a(8) = 2.903; 


% IL-1b regulation of MMP-1 (Yorifuji et al. 2016 fig 5b (Formerly, Moon et
% al. 2010 fig. 6a))
IL1_fit = [0,1]; % ng/mL
MMP1_IL1 = [4.167,101.6];
MMP1_IL1 = (MMP1_IL1-MMP1_IL1(1))/MMP1_IL1(1);
a(9) = 46.76; 

% IL-1b regulation of MMP-3 (Fuchs. S, et al. 2004)
IL1_fit =[0,10]; % ng/ml
MMP3_IL1= [12.78, 1233.23]; % ng/ml
MMP3_IL1 = (MMP3_IL1-MMP3_IL1(1))/MMP3_IL1(1);
a(65) = 105;

% IL-6 regulation of MMP-3 (Fuchs. S, et al. 2004)
IL6_fit =[0,10]; % ng/ml
MMP3_IL6= [12.78, 479.23] ;% ng/ml
MMP3_IL6 = (MMP3_IL6-MMP3_IL6(1))/MMP3_IL6(1);
a(66) = 40.15;

% TNF-a regulation of MMP-3 (Fuchs. S, et al. 2004)
TNF_fit =[0,10] ;% ng/ml
MMP3_TNF= [12.78, 562.3] ;% ng/ml
MMP3_TNF = (MMP3_TNF-MMP3_TNF(1))/MMP3_TNF(1);
a(67) = 47.3;

%% Hormone feedback parameter fits
% Kou et al. 2015, E2 effect on macrophage polarization in TMJ
E2 = [84.08,341.38]./1000; % Nanomolar (male/post-menopausal range, pre-menopausal range) 
% Calippe et al. 2010, figure 2b
IL1_e2 = [2.850336, 4.4998865];
IL1_e2 = (IL1_e2-IL1_e2(1))/IL1_e2(1);
a(10) =  2.616; % a_e2il1


% D'Agostino et al. 1999, figure 1
E2_free = [0, 100]./1000; % Nanomolar (Free concentration)
E2_serum = E2_free./0.02; 
IL10_e2 = [696.1123, 347.64676];
IL10_e2 = IL10_e2/IL10_e2(1);
a(11) = -0.1389; % a_e2il10


% D'Agostino et al. 1999, T effect on IL-10, figure 1
T_free = [0,0.1]; % Nanomolar (Free concentration)
T_serum = T_free./0.02;
IL10_T = [693.88324, 995.73865];
IL10_T = (IL10_T-IL10_T(1))/IL10_T(1);
a(12) = 0.522; % a_til10 

% D'Agostino et al. 1999, T effect on TNF, figure 2
TNF_t = [15.001847,10.349695];
TNF_t = TNF_t/TNF_t(1);
a(13) = -0.07424; % a_ttnf

% Lei et al. 2014, P effect on TNF, figure 2
P_free = [0,100]./1000; % MICROmolar (Free concentration)
P_serum = P_free./0.02;
TNF_p = [99.88827, 78.65922];
TNF_p = TNF_p/TNF_p(1);
a(14) = -0.04779; % a_ptnf


%% Macrophage feedback parameter fits
%IL-10 effect on TNF-a x
a(15) = 0.4666;
a(16) = -1.528;
a(17) = 0.5332;

% IL10 effect on IL1b, Nagaraja f3, SI ref. 22
a(18) = 0.6334;
a(19) = -1.794;
a(20) = 0.3667;

% TGFb effect on TNFa, Nagaraja f4, SI ref. 23
a(21) = 0.6211;
a(22) = -0.8305;
a(23) = 0.4466;

% TGFb effect on IL1b, Nagaraja f5, SI ref. 23
a(24) = 0.69;
a(25) = -20.37;
a(26) = 0.31;

% TGFb effect on IL10, Nagaraja f10, SI ref. 23
a(27) = 274.5;

% TGFb effect on monocyte chemotaxis (from Nagaraja 2014)
a(28) = -240.29;
a(29) = 298.93;
a(30) = -0.5926;
a(31) = 60.593;

% TNFa effect on monocyte chemotaxis (from Nagaraja 2014)
a(32) = -0.3164;
a(33) = 10.708;

% From Saren et al. 1996, figure 2
% IL-1 up-regulation of MMP-9
IL_fit = [0,40];
MMP9_IL1 = [1.218,3.852];
MMP9_IL1 = (MMP9_IL1-MMP9_IL1(1))/MMP9_IL1(1); % Normalize the concentration data
a(34) = 4.325;

% TNF upregulation of MMP9
TNF_fit = [0,10,20];
MMP9_TNF = [1.218, 2.747,3.153];
MMP9_TNF = (MMP9_TNF-MMP9_TNF(1))/MMP9_TNF(1);
a(35) = 1.531;

% TNF-a regulation of TIMP1 production by macrophages (Saren et al. 1996)
TNF_fit = [0, 2, 20, 40]; % ng/mL
TIMP1_TNF = [0.181, 0.071, 0.071, 0.071];
TIMP1_TNF = TIMP1_TNF/TIMP1_TNF(1);

a(36)=0.6077;
a(37)=-3.683;
a(38)=0.3923;

% IL-1b regulation of TIMP1 production by macrophages (Saren et al. 1996)
IL1_fit = [0, 10, 20];
TIMP1_IL1 = [0.181,0.0813, 0.0716];
TIMP1_IL1 = TIMP1_IL1/TIMP1_IL1(1);

a(39)=0.6102;
a(40)=-0.233;
a(41)=0.3898;

% IL10 regulation of TIMP1 production by macrophages (Jovanovic et al.
% 2000)
IL10_fit = [0, 10]; % ng/mL
TIMP1_IL10_M1 = [194, 94.8];
TIMP1_IL10_M1 = (TIMP1_IL10_M1)/TIMP1_IL10_M1(1);

a(42) = 1;
a(43) = -0.07161;
a(44) = 0;

% IL-10 effect on macrophage transformation 

IL10=[0,0.01,0.1,1,10,100]; % Data digitized from Kuwata et al. 2003 figure 6b
TNF_IL10 = [2725.3408, 2184.0051, 1499.6345, 781.18054, 485.02225, 386.37488];

% The feedback function increases the magnitude of the number subtracted
% from dM1/dt equation -> The function should increase as M1 concentration
% decreases and as M2 concentration increases
M2_IL10 = TNF_IL10.^-1;
M2_IL10_n = (M2_IL10-M2_IL10(1))/M2_IL10(1);

a(45) = 0.3213;

% IL-6 effect on TNF-a production by macrophages (Nagaraja si ref 24)
a(46) = 4.488; 
a(47) = 0.1541;

% IL-6 effect on IL-1b production by macrophages (Nagaraja si ref 25)
a(48) = 4.459;
a(49) = 0.1571; 


% IL-10 effect on IL-6 production by macrophages (Nagaraja si ref 22)
a(50) = 0.3298;
a(51) = -1.189;
a(52) = 0.6695;


% TGF-b effect on IL-6 production (Nagaraja si ref 26)
a(53) =  0.9821;

% E2 effect on IL-6 Liu et al. 2014
E2_fit = [0, 1]./0.02; % Nanomolar, converted to "total" serum concentration from "free"
IL6_E2 = [10.9407, 2.0869]; % ng/mL
IL6_E2 = (IL6_E2)/IL6_E2(1);

a(54) = -0.03314;

% Pr effect on IL-6 Sun et al. 2012 (figure 1b at 24 hours)
P_fit = [0, 1e-7*1e6]./0.02; % MICROmolar, converted to "total" from "free" concentration
IL6_Pr = [0.4077, 0.2372]; % ng/mL
IL6_Pr = IL6_Pr/IL6_Pr(1);

a(55) = -0.1083;


% IL-1b effect on SF IL-6 Inuoe et al. 2001
IL1_fit = [0, 1]; % ng/mL
IL6_IL1 = [0.065, 10.12]; % ng/mL
IL6_IL1 = (IL6_IL1-IL6_IL1(1))/IL6_IL1(1);

a(56) = 1.987;

% TNF-a effect on SF IL-6 Mrosewski et al. 2014
TNF_fit = [0,10]; % ng/mL
IL6_TNF = [8,259]; % Relative gene expression
IL6_TNF = (IL6_TNF-IL6_TNF(1))/IL6_TNF(1);

a(57)=1.066;

% IL-6 effect on M1 IL-10 Kothari et al. 
IL6_fit = [0,1, 10, 25];%, 50, 100];
IL6_IL10 = [9.54, 15.54, 182.19, 345.3];%, 635.1, 641.1];
IL6_IL10 = (IL6_IL10-IL6_IL10(1))/IL6_IL10(1);

a(58) = 0.1424;

% IL-10 effect on M1 MMP-9 Kothari et al. 
IL10_fit = [0, 100];
IL10_MMP9 = [1.25, 0.512];
IL10_MMP9 = IL10_MMP9/IL10_MMP9(1);

a(59) = -0.008926;

% IL-6 effect on SF TIMP-1
IL6_fit = [0, 50]; % ng/mL
IL6_TIMP1 = [126.55, 173.33]; % ng/mL
IL6_TIMP1 = (IL6_TIMP1-IL6_TIMP1(1))/IL6_TIMP1(1);

a(60) = 0.2753;

%% Chondrocyte Feedback Parameter Fits


% %IL-1B effect on CH MMP-1 (Fu et al.2016)
IL1_fit= [0,10]; % ng/ml
IL1_MMP1ch= [96.15, 839.74 ]; %pg/ml
IL1_MMP1ch= (IL1_MMP1ch-IL1_MMP1ch(1))/IL1_MMP1ch(1);

a(61)= 0.6122;

% IL-1b up-regulation of CH MMP-13 (Dunn et al. 2014)
MMP13_IL1= [327.5, 1092]; % pg/ml
MMP13_IL1 = (MMP13_IL1-MMP13_IL1(1))/MMP13_IL1(1);
a(62) = 2.568 ;

%TNF-A effect on CH IL-6 (Kaneva et al. 2012)
TNF_fit= [0,20,40,60];% ng/ml
TNF_IL6= [30.48, 52.38, 62.86, 154.29]; % pg/ml
TNF_IL6= (TNF_IL6-TNF_IL6(1))/TNF_IL6(1);

a(63)= 0.0267 ;

% IL-1B effect on CH IL-6 (Kloesch et al. 2012)
IL1_fit= [0, 5];% ng/ml
IL1_IL6= [0.77, 1.93]; % 
IL1_IL6= (IL1_IL6-IL1_IL6(1))/IL1_IL6(1);

a(64)= 1.808;

% IL-1b up-regulation of CH MMP-3 (Fu et al. 2016)
IL1_fit =[0,10] ;% ng/ml
MMP3_IL1= [1, 4.36]; % pg/ml
MMP3_IL1 = (MMP3_IL1-MMP3_IL1(1))/MMP3_IL1(1);
a(68) = 3.696 ; 

% TGF-B up-regulation of CH TIMP (Siliacci et al. 1998)
TGF_fit= [0, 10];% ng/ml
TGF_TIMP= [0.28, 1.35]; %ug/ml
TGF_TIMP= (TGF_TIMP-TGF_TIMP(1))/TGF_TIMP(1);

a(69)= 4.204;

% IL_1 Suppression of TIMP-1 (Siliacci et al. 1998)
IL1_fit = [0,1];% ng/ml
IL1_TIMP = [0.28,0.2]; % ug/ml
IL1_TIMP = IL1_TIMP/IL1_TIMP(1);

a(70)= -0.337;

% IL-6 Upregulates MMP-13 (Hashizume et al. 2010)
IL6_fit = [0 , 100]; % ng/ml
IL6_MMP13 = [3.824 , 43.235]; % ng/ml
IL6_MMP13 = (IL6_MMP13-IL6_MMP13(1))/IL6_MMP13(1);

a(71)= 10.41;


% TNFa Upregulates IL-1 (Kaneva et al .2012)
TNF_fit= [0,20,40,60];% ng/ml
TNF_IL1= [3.2, 4.0, 2.68, 14.5]; % pg/ml
TNF_IL1= (TNF_IL1-TNF_IL1(1))/TNF_IL1(1);

a(72)= 0.019;

% Il-6 Upregulates TGFb (Villiger et al. 1993)
IL6_fit = [0 , 10]; % ng/ml
IL6_TGF = [1.2 , 6.0]; % ng/ml
IL6_TGF = (IL6_TGF-IL6_TGF(1))/IL6_TGF(1);

a(73)= 4.4;

% IL-1 Upregulates TGFb (Villiger et al. 1993)
IL1_fit = [0 , 10]; % ng/ml
IL1_TGF = [1.2 , 1.70]; % ng/ml
IL1_TGF = (IL1_TGF-IL1_TGF(1))/IL1_TGF(1);

a(74)= 0.4583 ;

% Il-6 Upregulates MMP3 (Hashizume et al. 2010)
IL6_fit = [0 , 100]; % ng/ml
IL6_MMP3 = [6.47 , 47.06]; % ng/ml
IL6_MMP3 = (IL6_MMP3-IL6_MMP3(1))/IL6_MMP3(1);

a(75)= 6.336 ;

% Il-6 Upregulates MMP1 (Hashizume et al. 2010)
IL6_MMP1 = [2.09 , 17.54]; % ng/ml
IL6_MMP1 = (IL6_MMP1-IL6_MMP1(1))/IL6_MMP1(1);

a(76)= 7.466 ;

