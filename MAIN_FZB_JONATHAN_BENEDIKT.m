%Flugzeugauslegung Software 0.0.19
%              17.06.2018
%Jonathan Nicolai & Benedikt Kniebel
%       Copyright 0000-2099
%--------------------------------------%
%               MAIN
%           TLARs EINGABE
%--------------------------------------%
%--------------------------------------%
%--------------------------------------%
%!! NUR FÜR FZ MIT WENIGER ALS 240 PAX !!
%       Beispieldaten B737-500 

PAX=140;                %[Passagiere]
Reichweite=5200;        %[Reichweite in km]
M_cargo=8000;           %[Fracht in kg]
v_cruise=250;           %[Cruisegeschwindigkeit in m/s]
h_cruise=10000;         %[Cruisehoehe in m]
n_Aisle=1;              %[Anzahl Gaenge]
n_Sa=6;                 %[Anzahl Sitze pro Reihe (2x3)]
p_a=0.8;                %[Sitzabstand in m]
FK=2;                   %[FlugzeugKlasse 1:ATR 2:A319 3:A321]
c_D_cruise0=0.05;       %[Cruise-Widerstandsbeiwert des Flugzeugs] 

%Triebwerk
AK=1;                   %[Antiebskonzept 1:Turbofan, 2:Dispursal]
SFC_konv=(1.4665*10^(-5));%[specific fuel cosumption]

%FLUEGEL
lambda_FL_InitialSizing=9.16; %[Fluegelstreckung angenommen]
b=28.88;                %[Spannweite]
SF=1.0;                 %[Sizing Factor SF*MTOW_Vergleich]
SF_MLW=0.79;            %[Sizing Landing Factor: MLW= SF_MLW*MTOW_Vergleich] 
c_lProfil=0.5;          %[Design-Auftriebsbeiwert des Profils] 
c_L_cruise=0.6;         %[Cruise-Auftriebsbeiwert des Profils] 
c_l_maxto=2.2;          %[Maximaler Auftriebskoeffizienz TO]
c_l_maxland=2.7;        %[Maximaler Auftriebskoeffizienz Landing]
k_r_HTP=0.55;           %[Abstand NP-Fluegel Rumpf in %]
BFL=2200;               %[Startstrecke in m]
BPR=5;                  %[Bypass ratio]

t_lam_wing_HTP_VTP=0.2; %[Laminaranteil]
t_lam_fslg_nacelle=0;   %[Laminaranteil]

%------------------------------------------------
%x_lmag_wing=13;         %[??? Schwerpunkt]
%l_gamma_wing=10;        %[???]
%x_mac=1;                %[Abstand vorne <-> mittlere Fluegeltiefe]

t_c=0.2;                %[Dickenverhältnis]
x_t_wing=1;             %[relative Position der max. Dickenrücklage]
t_c_root_wing=0.2;      %[Dickenverhältnis Wurzel]
t_c_tip_wing=0.2;       %[Dickenverhältnis Spitze]
phi_m_wing=4;           %[Pfeilungswinkel an Position der max. Dickenrücklage]

%NACA0012 
x_t_HTP=0.4;            %[relative Position HTP der max. Dickenrücklage]
phi_m_HTP=36.55;        %[Pfeilungswinkel an Position HTP der max. Dickenrücklage]
t_c_HTP=0.2;            %[Dickenverhältnis HTP MAC]
t_c_root_HTP=0.2;       %[Dickenverhältnis Wurzel]
t_c_tip_HTP=0.2;        %[Dickenverhältnis Spitze]
lambda_zu_HTP=5;        %[Zuspitzung HTP]

%NACA0012 
x_t_VTP=0.4;            %[relative Position VTP der max. Dickenrücklage]
phi_m_VTP=36.55;        %[Pfeilungswinkel an Position VTP der max. Dickenrücklage]
t_c_VTP=0.2;            %[Dickenverhältnis VTP MAC]
t_c_root_VTP=0.2;       %[Dickenverhältnis Wurzel]
t_c_tip_VTP=0.2;        %[Dickenverhältnis Spitze]
lambda_zu_VTP=5;        %[Zuspitzung VTP]
%--------------------------------------------------------------------------

%Start/Landung
s_land=1800;            %[Landestrecke in m]
a_g=0.45;               %[Verzoegerung waehrend des Rollens in m/s2]
delta_y2=0;             %[Differenz zwischen dem gewünschten Steigwinkel im 2. Steigsegment und dem kleinsten zulässigen Winkel nach FAR 25]

%Sonstiges:
p_0=101325;       %[Druck am Boden in Pascal]
T_0=288.15;             %[Temperatur am Boden in Kelvin]
rho_0=1.225;            %[Dichte bei Hoehe 0 in kg/m3]
R=287.05;               %[J/(Kg*K)]
k=6.34*10^(-6);         %[Rauhigkeitskennwert]
Re_cruise=4500000;      %[Flugreynoldszahl]
e=0.85;                 %[Oswaldfaktor]
K_v=1.1;                %[Konstante velocity: M_crit3D = K_v_FA * v_c_FA]

%--------------------------------------%
%--------------------------------------%

       %Zu berechnende Variablen

%--------------------------------------%
%--------------------------------------%

%Massen
syms OEW;               %[Operational Empty Weight]
syms MTOW_vor;          %[Maximales Abfluggewicht vorlaeufige in kg]
syms MLW;               %[Maximales Landegewicht in kg]

syms m_wing;            %[Fluegelgewicht in kg]
syms m_fslg;            %[Rumpfgweicht in kg]
syms m_HTP;             %[Masse HTP in kg]
syms m_VTP;             %[Masse VTP in kg]
syms m_NG;              %[Masse Nose Gear in kg]
syms m_MG;              %[Masse Main Gear in kg]
syms m_prop;            %[Antrieb in kg]
syms m_sys;             %[Systemgewicht in kg]
syms m_furnish;         %[Einrichtungsgewicht in kg]
syms m_optItems;        %[Optinale Items Gewicht in kg]

syms m_fuel;            %[Fuelgewicht in kg]

%Schwerpunkte
syms x_cg;              %[Geamtschwerpunkt]
syms x_cg_wing;         %[Schwerpunkt Wing]
syms x_cg_HTP;          %[Schwerpunkt HTP]

%Cabin
syms d_cabin_innen;     %[Kabineninnendurchmesser]
syms d_cabin_aussen;    %[Kabinenaussendurchmesser]
syms l_fslg;            %[Kabinenlaenge]
syms l_front;           %[Kabinenlaenge front]
syms l_tail;            %[Kabinenlaenge tail]

%Fluegel
syms l_root;            %Fluegelwurzel
syms l_tip;             %Fluegelspitze
syms S_ref;             %[Referenzfluegelflaeche]
syms W_S;               %[Flächenbelastung]
syms phi_25;            %[Flügelpfeilung in Grad]
syms l_gamma;           %[Referenzfluegeltiefe]
syms y_mac;             %[Spannweitige Lage der Bezugs?ügeltiefe]
syms lambda_Fl;         %[Fluegelstreckung]
syms lambda_zu_FL;      %[Fluegelzuspitzung]
syms l_gamma;           %[Bezugsfluegeltiefe]
syms T_W_design;        %[]
syms S_HTP;             %[Hoehenleitwerkflaeche]
syms S_VTP;             %[HSeitenleitwerkflaeche]
syms v_HTP;             %[Volumenkoef?zient]
syms v_VTP;             %[Volumenkoef?zient]
syms l_gamma_HTP;       %[Hoehenleitwerk Referenzfluegeltiefe]
syms l_gamma_VTP;       %[Seitenleitwerk Referenzfluegeltiefe]

%Triebwerke
syms d_scaled;          %[Triebwerksdurchmesser]
syms T_W_to;            %[Thrust/Weight to]
syms T_W;               %[Thrust/Weight Design]
syms l_scaled;          %[Triebwekslaenge]
syms w_scaled;          %[Triebweksgewicht]
syms T_to_req_TA;       %[Benoenigter Triebweksschub to]
syms T_cruise_req_TA;   %[Benoenigter Triebweksschub cruise]
syms Faktor_over_req;   %[Faktor Schub Referenz/Berechnet]
syms laengenfaktor_TW;  %[Faktor Durchmesser Referenz/Berechnet]
syms duruchmesserfaktor_TW; %[Faktor Laenge Referenz/Berechnet]

%Hoehe/Druck ect.
syms rho_cruise;        %[Dichte Reiseflughoehe]
syms p_cruise;          %[Druck Reiseflughoehe]
syms Temp_cruise;       %[Temperatur Reiseflughoehe]

%Sonstiges
syms n_crew;            %[Anzahl Crewmitglieder]
syms n_toi;             %[Anzahl Toiletten]
syms Ma_cruise;         %[Berechnete Flugmachzahl]

% S_exposed_wing;       %[Fluegelflaeche ohne Rumpf]
% S_exposed_HTP;        %[HTP-Flaeche ohne Rumpf]
% S_exposed_VTP;        %[VTP-Flaeche ohne Rumpf]
% S_HTP_ref;               
% S_VTP_ref;                


%--------------------------------------%
%--------------------------------------%

%----------Programmablauf:-------------%

%--------------------------------------%
%--------------------------------------%

clc %Konsole aufraeumen

%--------Berechnung Konstanten--------

%Berechnung von Druck, Temperatur ect.
[p_cruise, Temp_cruise, rho_cruise]=Vorauslegung.Tempauslegung(p_0, T_0, rho_0, h_cruise);

%Berechnung von Flugmachzahl
[Ma_cruise]=Vorauslegung.Machzahl(v_cruise, R, Temp_cruise);

%--------Grundberechnung einmalig--------

%Rumpf
[d_cabin_innen,l_fslg, l_front, l_tail, n_toi, n_crew, d_cabin_aussen, V_Cabin, V_Cargo]=Vorauslegung.Rumpfauslegung(PAX, n_Aisle, n_Sa, p_a);     

%Gewichtsabschätzung
[MTOW_vor, MLW_vor]=Vorauslegung.VORGewichtsabschaetzung(FK,SF, SF_MLW);  

%Startstrecke
[W_S, T_W]=Vorauslegung.Startstrecke(BFL, rho_0, delta_y2, c_l_maxto, BPR, a_g, c_l_maxland, s_land, MTOW_vor, MLW_vor, e, lambda_FL_InitialSizing, v_cruise, rho_cruise, Ma_cruise);

%--------Schleifendurchlauf:--------
for index=1:100
MTOW_Ber = MTOW_vor(index);
MLW_Ber = MLW_vor(index);

%Fluegel
[S_ref, phi_25, l_gamma, y_mac, lambda_Fl,lambda_zu_FL, l_root, l_tip]=Vorauslegung.Fluegelauslegung(MTOW_Ber, W_S, b, K_v, Ma_cruise, t_c, c_lProfil);

%Leitwerk
[v_HTP,  v_VTP, S_VTP, S_HTP, lambda_Hlw, lambda_Slw, lambda_zu_Hlw, lambda_zu_Slw, l_gamma_HTP, l_gamma_VTP]=Vorauslegung.Leitwerksauslegung(k_r_HTP, l_fslg,lambda_Fl, lambda_zu_FL, b, l_gamma, S_ref);

%Triebwerk
[T_overall_req, T_Antrieb1 , T_Antrieb2, d_scaled, l_scaled, laengenfaktor_TW, duruchmesserfaktor_TW, Faktor_over_req, M_TW]=Vorauslegung.Triebwerksauslegung(AK, T_W, rho_cruise, v_cruise, lambda_Fl, e, p_cruise, p_0, T_0, Temp_cruise, c_L_cruise, MTOW_Ber, c_D_cruise0, S_ref);
 
%Fahrwerk
%Erst einmal nicht enthalten
 
%Aerodynamik
[L_D_max, C_L_opt, C_D_0]=Vorauslegung.Aerodynamik(lambda_Fl, e, Re_cruise, l_fslg, l_gamma, l_gamma_HTP, l_gamma_VTP, l_scaled, d_cabin_aussen, t_c, t_lam_wing_HTP_VTP, t_lam_fslg_nacelle, k, d_scaled, l_front, l_tail, Ma_cruise, phi_m_wing, x_t_wing, S_ref, x_t_HTP, x_t_VTP, phi_m_HTP, phi_m_VTP, lambda_zu_FL, t_c_root_wing, t_c_tip_wing, t_c_root_HTP, t_c_tip_HTP, S_HTP, S_VTP, lambda_zu_HTP, lambda_zu_VTP, AK,laengenfaktor_TW, duruchmesserfaktor_TW, t_c_root_VTP, t_c_tip_VTP, t_c_VTP, t_c_HTP);
 
%Massenberechnug
[OEW, ZFW, d_fslg_aussen, m_NG, m_MG, m_wing, m_fslg, m_HTP, m_VTP, m_prop, m_sys, m_optItems, m_furnish]=Vorauslegung.Gewichtsabschaetzung(d_cabin_innen, l_fslg, l_front, l_tail, lambda_Fl, MTOW_Ber, S_ref, lambda_zu_FL, phi_25, AK, S_HTP, S_VTP, MLW_Ber, n_crew, n_toi, PAX, V_Cargo, V_Cabin, h_cruise, M_TW, M_cargo);

%Fuelplanning
[m_fuel]=Vorauslegung.Fuelplanning(AK, T_Antrieb1, T_Antrieb2, Reichweite, SFC_konv, ZFW, v_cruise, L_D_max);

%[Schwerpunktsberechnung]
% [x_cg, x_cg_wing, x_cg_HTP]=Vorauslegung.Schwerpunktsberechnung(AK, S, b, x_lmag_wing, l_gamma_wing, l_fslg, x_mac, phi_25,m_wing, m_fslg, m_HTP, m_VTP, m_NG, m_MG, m_prop, m_sys, m_furnish, m_optItems, m_fuel, M_cargo);

%Vektor mit berechneten Differenzwerten in den Zeilen

Diff(index)=abs(MTOW_Ber-(OEW + m_fuel + M_cargo + PAX * 102));
disp(MTOW_Ber);
disp('Wert wird berechnet:');
disp(index);
end

%-------------------------------------------
 
%Minimierungsfunktion

[i, index_ausgewaehlt] = (min(Diff));
MTOW_ENDGUELTIG = MTOW_vor(index_ausgewaehlt);

%Konsolenausgabe
disp('Das entgueltige MTOW lautet: ');
disp(MTOW_ENDGUELTIG);
disp('Dies ist Wert:');
disp(index_ausgewaehlt);

%Differenz cd0= und c_DO Aero