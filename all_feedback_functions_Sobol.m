function f = all_feedback_functions_Sobol(IL1b,TNFa,IL10,TGFb,TGFb_del,estr,pr,test,IL6)
% Calculation of feedback functions. This function calls the necessary
% coefficients from "all_parameter_fits" and is called by "ode_ftn." These
% get recalculated for every integration step. 

global Ksf

coeffs = Ksf;

% Assign meaningful names to the coefficients

% SF coefficients
a_tnfmmp1=coeffs(1);
a_tnftimp1s=coeffs(2);
a_tnfmmp13=coeffs(3);
a_tnfil1=coeffs(4);
a_tnftgf=coeffs(5);
a_il1timp1s=coeffs(6);
a_il1il1=coeffs(7);
a_il1tnf=coeffs(8);
a_il1mmp1sf=coeffs(9);
a_il1mmp3sf=coeffs(65);
a_il6mmp3sf=coeffs(66);
a_tnfmmp3sf=coeffs(67);

% Hormone coefficients
a_e2il1 = coeffs(10);
a_e2il10 = coeffs(11);
a_til10 = coeffs(12);
a_ttnf = coeffs(13);
a_ptnf = coeffs(14);

% Macrophage coefficients
a_il10tnf = coeffs(15);
b_il10tnf = coeffs(16);
c_il10tnf = coeffs(17);
a_il10il1 = coeffs(18);
b_il10il1 = coeffs(19);
c_il10il1 = coeffs(20); 
a_tgftnf = coeffs(21);
b_tgftnf = coeffs(22);
c_tgftnf = coeffs(23);
a_tgfil1 = coeffs(24);
b_tgfil1 = coeffs(25);
c_tgfil1 = coeffs(26);
a_tgfil10 = coeffs(27);
r1q = coeffs(28);
r2q = coeffs(29);
r1l = coeffs(30);
r2l = coeffs(31);
r1 = coeffs(32);
r2 = coeffs(33);
a_il1mmp9 = coeffs(34);
a_tnfmmp9 = coeffs(35);
a_tnftimp1=coeffs(36);
b_tnftimp1=coeffs(37);
c_tnftimp1=coeffs(38);
a_il1timp1=coeffs(39);
b_il1timp1=coeffs(40);
c_il1timp1=coeffs(41);
a_il10timp1=coeffs(42);
b_il10timp1=coeffs(43);
c_il10timp1=coeffs(44);
a_m1m2 = coeffs(45);
a_il6tnf = coeffs(46);
b_il6tnf = coeffs(47);
a_il6il1 = coeffs(48);
b_il6il1 = coeffs(49);
a_il10il6 = coeffs(50);
b_il10il6 = coeffs(51);
c_il10il6 = coeffs(52);
a_tgfil6 = coeffs(53);
a_e2il6 = coeffs(54);
a_pril6 = coeffs(55);
a_il1il6 = coeffs(56);
a_tnfil6 = coeffs(57); 
a_il6il10 = coeffs(58);
a_il10mmp9 = coeffs(59);
a_il6timp1 = coeffs(60);

%Chondrocyte Coeefecients 
a_il1mmp1ch = coeffs(61);
a_il1mmp13ch = coeffs(62);
a_tnfil6ch = coeffs(63);
a_il1il6ch = coeffs(64);
a_il1mmp3ch = coeffs(68);
a_tgftimpch = coeffs(69);
a_il1timpch = coeffs(70);
a_il6mmp13ch = coeffs(71);
a_tnfil1ch = coeffs(72);
a_il6tgfch = coeffs(73);
a_il1tgfch = coeffs(74);
a_il6mmp3ch = coeffs(75);
a_il6mmp1ch = coeffs(76);

%% SF functions

% TNF effect on MMP-1
f(1) = a_tnfmmp1*TNFa/(1+TNFa);

% TNF effect on TIMP-1
f(2) = a_tnftimp1s*TNFa/(1+TNFa);

%TNF effect on MMP-13
f(3) = a_tnfmmp13*TNFa/(1+TNFa);

% TNF effect on IL-1
f(4) = a_tnfil1*TNFa/(1+TNFa);

% TNF effect on TGF-b
f(5) = a_tnftgf*TNFa/(1+TNFa);

% IL-1 effect on TIMP-1
f(6) = a_il1timp1s*IL1b/(1+IL1b);

% IL-1 effect on IL-1
f(7) = a_il1il1*IL1b/(1+IL1b);

% IL-1 effect on TNF-a
f(8) = a_il1tnf*IL1b/(1+IL1b);

% IL-1 effect on MMP-1
f(9) = a_il1mmp1sf*IL1b/(1+IL1b);

% IL1 effect on IL-6 in SF
f(34) = a_il1il6*IL1b/(1+IL1b);

% TNFa effect on IL-6 in SF
f(35) = a_tnfil6*TNFa/(1+TNFa);

% IL1 effect on MMP-3 in SF
f(56) = a_il1mmp3sf*IL1b/(1+IL1b);

% IL6 effect on MMP-3 in SF
f(57)= a_il6mmp3sf*IL6/(1+IL6);

% TNF effect on MMP-3 in SF
f(58) = a_tnfmmp3sf*TNFa/(1+TNFa);

%% Hormone functions

% Estrogen effect on M1 IL-1 production

f(10) = a_e2il1*estr/(1+estr); % 

% Estrogen effect on M1 IL-10 production
f(11) = exp(a_e2il10*estr);

% Testosterone effect on M1 IL-10 production
f(12) = a_til10*test/(1+test);
 
% Testosterone effect on M1 TNFa production
f(13) = exp(a_ttnf*test);

% Progesterone effect on M1 TNFa production
f(14) =exp(a_ptnf*pr);


%% Macrophage functions

% IL-10 effect on TNFa, Nagaraja f1, SI ref. 22
f(15) =  a_il10tnf*exp(b_il10tnf*IL10)+c_il10tnf; 

% IL10 effect on IL1b, Nagaraja f3, SI ref. 22
f(16) = a_il10il1*exp(b_il10il1*IL10)+c_il10il1;  

% TGFb effect on TNFa, Nagaraja f4, SI ref. 23
f(17) = a_tgftnf*exp(b_tgftnf*TGFb)+c_tgftnf; 

% TGFb effect on IL1b, Nagaraja f5, SI ref. 23
f(18) =a_tgfil1*exp(b_tgfil1*TGFb)+c_tgfil1; 

% TGFb effect on IL10, Nagaraja f10, SI ref. 27
f(19) = a_tgfil10*TGFb/(1+a_tgfil10*TGFb); 


% TGFb effect on monocyte chemotaxis
TGFb_del=TGFb_del*1e3; % Convert to pg/mL (because that's what Nagaraja did in their code)
if TGFb_del <= 1
    f(20) = r1q*TGFb_del^2+r2q*TGFb_del;
elseif TGFb_del > 1 && TGFb_del <= 10
    f(20) = r1l*TGFb_del+r2l;
else
    f(20) = 0;
end

% TNFa effect on monocyte chemotaxis
f(21) = r1*TNFa^2+r2*TNFa;  

% IL1b effect on MMP-9 (from Jager et al. 2016)
f(22) = a_il1mmp9*IL1b/(1+IL1b);

% TNFa effect on MMP-9
f(23) = a_tnfmmp9*TNFa/(1+TNFa);

% TNFa effect on TIMP-1 production by M1
f(24) = a_tnftimp1*exp(b_tnftimp1*IL1b)+c_tnftimp1; 

% IL-1 effect on TIMP-1 production by M1
f(25) = a_il1timp1*exp(b_il1timp1*IL1b)+c_il1timp1; 

%IL-10 effect on TIMP-1 production by M1
f(26) = a_il10timp1*exp(b_il10timp1*IL10)+c_il10timp1;

% IL-10 effect on macrophage transformation
f(27) = a_m1m2*IL10/(1+a_m1m2*IL10);

% IL-6 effect on TNF-a production by M1
f(28) = a_il6tnf*TNFa/(a_il6tnf+TNFa^b_il6tnf);

% IL-6 effect on IL-1b production by M1
f(29) = a_il6il1*IL1b/(a_il6il1+IL1b^b_il6il1);

% IL-10 effect on IL-6 production by M1
f(30) = a_il10il6*exp(b_il10il6*IL10)+c_il10il6;

% TGFb effect on IL-6 production by M1
f(31) = a_tgfil6*TGFb/(1+TGFb);

% E2 effect on IL-6 production by M1
f(32) = exp(a_e2il6*estr);

% Pr effect on IL-6 production by M1
f(33) = exp(a_pril6*pr);

% IL-6 effect on M1 IL-10 *
f(36) = a_il6il10*IL6/(1+a_il6il10*IL6);

% IL-10 effect on M1 MMP-9
f(37) = exp(a_il10mmp9*IL10);

% IL-6 effect on SF TIMP-1
f(38) = a_il6timp1*IL6/(1+IL6);

%% Chondrocyte functions

% IL-1B effect on CH MMP-1
f(43)= a_il1mmp1ch*IL1b/(1+IL1b);
 
%IL-1b effect on CH MMP-13
f(44)= a_il1mmp13ch*IL1b/(1+IL1b);

%IL-1b effect on CH MMP-3
f(45)= a_il1mmp3ch*IL1b/(1+IL1b);
 
%TNFa effect on CH IL-6
f(46)= a_tnfil6ch*TNFa/(1+TNFa);
 
% IL-1B effect on CH IL-6
f(47)= a_il1il6ch*IL1b/(1+IL1b);
 
%TGF-B effect on CH TIMP *
f(48)= a_tgftimpch*TGFb/(1+a_tgftimpch*TGFb);
 
%IL-1b effect on CH TIMP-1
f(49)= exp(a_il1timpch*IL1b);
 
%IL-6 effect on CH MMP-13
f(50) = a_il6mmp13ch*IL6/(1+IL6);
 
%TNFa effect on CH IL-1
f(51)= a_tnfil1ch*TNFa/(1+TNFa);
 
%IL-6 effect on CH TGFb
f(52) = a_il6tgfch*IL6/(1+IL6);
 
% IL-1B effect on CH TGFb
f(53)= a_il1tgfch*IL1b/(1+IL1b);
 
%Il6 effect on CH MMP3
f(54) = a_il6mmp3ch*IL6/(1+IL6);

%Il6 effect on CH MMP1
f(55) = a_il6mmp1ch*IL6/(1+IL6);

