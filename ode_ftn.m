function [X] = ode_ftn(t,x,Z,Ks)
% This function solves the system of differential equations. 
% It receives input from either "main_script" and/or "lhs_script_general." 
% Initial conditions are specified in "cytokine_hist."
% It calls "all_feedback_functions" to calculate feedback regulation
% Some variables are passed globally because they need to be manipulated
% outside this function

% This script prepared by Bethany Luke & Conner Hutcherson
% UT Southwestern Medical Center, 2023
% mail to: conner.hutcherson@utsouthwestern.edu

global Ks e2 estr0 test pr 

% Assign meaningful names to the model variables

P = x(1);
M1 = x(2);
M2 = x(3);
IL1b = x(4);
TNFa = x(5);
IL10 = x(6);
TGFb = x(7);
MMP9 = x(8);
MMP1 = x(9);
TIMP1 = x(10);
IL6 = x(11);
MMP13 = x(12);
MMP3 = x(13);
SFs = x(14);
CHs = x(15);
TGFb_del = Z(7);

%% Calculate hormone concentrations

% Current estrogen concentration
if length(e2) > 1
estr = ppval(spline(e2(:,1),e2(:,2)),t); % estr0; %
% subplot(2,1,1);
% plot(t,estr,'ok'), hold on;

    if estr < estr0
        estr=estr0;
    end
else
    estr = e2;
end

% Current progesterone concentration
if length(pr) > 1
    prog = ppval(spline(pr(:,1),pr(:,2)),t); % estr0; 
%     subplot(2,1,2);
%     plot(t,prog,'ok'), hold on;
else
    prog = 0;
end
% Current testosterone concentration
if length(test) > 1
    testosterone = ppval(spline(test(:,1),test(:,2)),t); % estr0; 
%     subplot(2,1,2);
%     plot(t,test,'ok'), hold on;
else
    testosterone = 0;
end

%% Calculate the feedback functions
f = all_feedback_functions(IL1b,TNFa,IL10,TGFb,TGFb_del,estr,prog,testosterone,IL6);

% SF feedback functions
f_tnfmmp1 = f(1);
f_tnftimp1s = f(2);
f_tnfmmp13 = f(3);
f_tnfil1 = f(4);
f_tnftgf = f(5);
f_il1timp1s = f(6);
f_il1il1 = f(7);
f_il1tnf = f(8);
f_il1mmp1 = f(9);
f_il1il6 = f(34); 
f_tnfil6 = f(35);
f_il1mmp3s = f(56);
f_il6mmp3s = f(57);
f_tnfmmp3s = f(58);

% Hormone feedback functions
f_e2il1 = f(10);
f_e2il10 = f(11);
f_til10 = f(12);
f_ttnf = f(13);
f_ptnf = f(14);

% Macrophage feedback functions
f_il10tnf = f(15);
f_il10il1 = f(16);
f_tgftnf = f(17);
f_tgfil1 = f(18);
f_tgfil10 = f(19);
f_tgfm1 = f(20);
f_tnfm1 = f(21);
f_il1mmp9 = f(22);
f_tnfmmp9 = f(23);
f_tnftimp1 = f(24);
f_il1timp1 = f(25);
f_il10timp1 = f(26);
f_m1m2 = f(27);
f_il6tnf = f(28);
f_il6il1 = f(29);
f_il10il6 = f(30);
f_tgfil6 = f(31);
f_e2il6 = f(32);
f_pril6 = f(33);
f_il6il10 = f(36);
f_il10mmp9 = f(37); 
f_il6timp1 = f(38);

% % Chondrocyte feedback functions
f_il1mmp1ch = f(43);
f_il1mmp13ch = f(44);
f_il1mmp3ch = f(45);
f_tnfil6ch = f(46);
f_il1il6ch = f(47);
f_tgftimpch = f(48);
f_il1timpch = f(49);
f_il6mmp13ch = f(50);
f_tnfil1ch = f(51);
f_il6tgfch = f(52);
f_il1tgfch = f(53);
f_il6mmp3ch = f(54);
f_il6mmp1ch = f(55);

%% Assign meaningful names to the rate coefficients

kd_P = Ks(1);
kM_in = Ks(2);
kd_M = Ks(3);
k_M1M2 = Ks(4);
kIL1b_M1 = Ks(5); 
kIL1b_M2 = Ks(6);
kIL1b_SF = Ks(7);
kIL1b_CH = Ks(8);
kIL1b_d = Ks(9);
kTNFa_M1 = Ks(10); 
kTNFa_M2 = Ks(11); 
kTNFa_SF = Ks(12);
kTNFa_CH = Ks(13);
kTNFa_d = Ks(14);
kIL10_M1 = Ks(15); 
kIL10_M2 = Ks(16); 
kIL10_SF = Ks(17);
kIL10_CH = Ks(18);
kIL10_d = Ks(19);
kTGFb_CH= Ks(20);
kTGFb_M1 = Ks(21); 
kTGFb_M2 = Ks(22); 
kTGFb_SF = Ks(23);
kTGFb_d = Ks(24); 
kMMP9_M1 = Ks(25);
kMMP9_M2 = Ks(26);
kMMP9_d = Ks(28);
kMMP1_M1 = Ks(29);
kMMP1_SF = Ks(30); 
kMMP1_CH = Ks(31);
kMMP1_d = Ks(32);
kTIMP1_M1 = Ks(33); 
kTIMP1_M2 = Ks(34);  
kTIMP1_SF = Ks(35); 
kTIMP1_CH = Ks(36);
kTIMP1_d = Ks(37);
kIL6_M1 = Ks(38);
kIL6_M2 = Ks(39);
kIL6_SF = Ks(40);
kIL6_d = Ks(41);
kIL6_CH = Ks(42);
kMMP13_CH = Ks(43);
kMMP13_SF = Ks(44);
kMMP13_d = Ks(45);
kTGF_P = Ks(46);
kMMP3_CH = Ks(47);
kMMP3_SF = Ks(48);
kMMP3_d = Ks(49);

%% Solve model equations

% Platelet equation
dP_dt = -kd_P*P;

% Don't let SFs produce inflammatory factors until macrophages start to
% invade
M_switch=P>10e-12;
% M_switch=0;

% M1 equation
dM1_dt = kM_in*M_switch*(f_tgfm1+f_tnfm1) - k_M1M2*f_m1m2*M1 - kd_M*M1;

% M2 equation
dM2_dt = k_M1M2*f_m1m2*M1-kd_M*M2;

% IL-1b equation
dIL1b_dt = kIL1b_M1*(1+f_il6il1)*f_il10il1*f_tgfil1*(1+f_e2il1)*M1 + kIL1b_M2*M2 + kIL1b_SF*(1+f_tnfil1)*(1+f_il1il1)*SFs + kIL1b_CH*(1+f_tnfil1ch)*CHs - kIL1b_d*IL1b; 

% TNFa equation
dTNFa_dt = kTNFa_M1*(1+f_il6tnf)*f_il10tnf*f_tgftnf*f_ttnf*f_ptnf*M1 + kTNFa_M2*M2 + kTNFa_SF*(1+f_il1tnf)*SFs + kTNFa_CH*CHs  - kTNFa_d*TNFa; 

% IL-10 equation
dIL10_dt = kIL10_M1*(1+f_tgfil10)*f_e2il10*(1+f_til10)*(1+f_il6il10)*M1 + kIL10_M2*M2 + kIL10_SF*SFs + kIL10_CH*CHs - kIL10_d*IL10;%

% TGFb equation
dTGFb_dt = kTGF_P*P+kTGFb_M1*M1 + kTGFb_M2*M2 + kTGFb_SF*(1+f_tnftgf)*SFs + kTGFb_CH*(1+f_il6tgfch)*(1+f_il1tgfch)*CHs - kTGFb_d*TGFb; 

% MMP9 equation
dMMP9_dt = kMMP9_M1*(1+f_il1mmp9)*f_il10mmp9*(1+f_tnfmmp9)*M1 + kMMP9_M2*M2 - kMMP9_d*MMP9;

% MMP-1 equation
dMMP1_dt = kMMP1_M1*M1 + kMMP1_SF*(1+f_tnfmmp1)*(1+f_il1mmp1)*SFs + kMMP1_CH*(1+f_il1mmp1ch)*(1+f_il6mmp1ch)*CHs - kMMP1_d*MMP1;

% TIMP-1 equation
dTIMP1_dt = kTIMP1_M1*f_tnftimp1*f_il1timp1*f_il10timp1*M1 + kTIMP1_M2*M2 + kTIMP1_SF*(1+f_tnftimp1s)*(1+f_il1timp1s)*(1+f_il6timp1)*SFs + kTIMP1_CH*(1+f_tgftimpch)*f_il1timpch*CHs - kTIMP1_d*TIMP1;% 

% IL-6 equation
dIL6_dt = kIL6_M1*f_il10il6*(1+f_tgfil6)*f_e2il6*f_pril6*M1+kIL6_M2*M2+kIL6_SF*(1+f_il1il6)*(1+f_tnfil6)*SFs + kIL6_CH*(1+f_il1il6ch)*(1+f_tnfil6ch)*CHs - kIL6_d*IL6; 

% MMP-13 equation
dMMP13_dt = kMMP13_SF*(1+f_tnfmmp13)*SFs + kMMP13_CH*(1+f_il1mmp13ch)*(1+f_il6mmp13ch)*CHs - kMMP13_d*MMP13;

% MMP-3 equation
dMMP3_dt = kMMP3_SF*(1+f_il1mmp3s)*(1+f_il6mmp3s)*(1+f_tnfmmp3s)*SFs + kMMP3_CH*(1+f_il1mmp3ch)*(1+f_il6mmp3ch)*CHs - kMMP3_d*MMP3;

% Synovial fibroblast equation
dSFs_dt = 0;

% Chondrocyte equation
dCHs_dt = 0;


X = zeros(15,1);

X(1) = dP_dt;
X(2) = dM1_dt;
X(3) = dM2_dt;
X(4) = dIL1b_dt;
X(5) = dTNFa_dt;
X(6) = dIL10_dt;
X(7) = dTGFb_dt;
X(8) = dMMP9_dt;
X(9) = dMMP1_dt;
X(10) = dTIMP1_dt;
X(11) = dIL6_dt;
X(12) = dMMP13_dt;
X(13) = dMMP3_dt;
X(14) = dSFs_dt;
X(15) = dCHs_dt;






