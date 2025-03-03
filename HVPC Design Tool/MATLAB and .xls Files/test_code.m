clc; clear variables; close all;

%INPUTS
Po = 600; % [W]
fs = 500000; % [Hz]
Vsec = 7.5e3; % [V]
Vpri = 500; % [V]
kf = 4.44; % wavefor coeff, 4.44 for sine wave
Bm = 0.057; % Assumed
Jp = 500; % Current density, assumed 500 [A/cm^2]
kp = 0.1; % Packing factor, assumed

%Core Data
Ap = 19.07e-8; %m^4  Area product
Ac = 3.677e-4; %m^2  
Le = 13.9e-2; %m
Wa = 5.186e-4; %m^2
H = 4.5e-2; %m
W = 2.305e-2; %m
MLT = 12.9e-2; %m
Pcore = 3.9; %W
Wcore = 0.13; %kg

Np = round(Vpri/(kf*Bm*Ac*fs));

Ns = round(Np*(Vsec/Vpri));

% Wire

   % Calculate Copper Loss (neglegable compared to core loss)

    % PriKlayer = sqrt(pi.*Pri_Nstrands).*ds./2./(Pri_WireSize);
    % Pri_xp = ds./2./skindepth.*sqrt(pi.*PriKlayer);
    % SecKlayer = sqrt(pi.*Sec_Nstrands).*ds./2./(Sec_WireSize);
    % Sec_xp = ds./2./skindepth.*sqrt(pi.*SecKlayer);
    % Pri_Rdc = rou.*TLp./(pi.*Pri_WireSize.^2./4);
    % Sec_Rdc = rou.*TLs./(pi.*Sec_WireSize.^2./4);
    % Pri_Fr = Pri_xp.*((sinh(2.* Pri_xp) + sin(2.*Pri_xp))./(cosh(2.*Pri_xp) - cos(2.*Pri_xp)) ...
    %     + 2.*(Mlp.^2.*Pri_Nstrands - 1)./3.*(sinh(Pri_xp) - sin(Pri_xp))./(cosh(Pri_xp) + cos(Pri_xp)));
    % Pri_Rac = Pri_Rdc.*Pri_Fr;
    % Sec_Fr = Sec_xp.*((sinh(2.*Sec_xp) + sin(2.*Sec_xp))./(cosh(2.*Sec_xp) - cos(2.*Sec_xp)) ...
    %     + 2.*(Mls.^2.*Sec_Nstrands - 1)./3.*(sinh(Sec_xp) - sin(Sec_xp))./(cosh(Sec_xp) + cos(Sec_xp)));
    % Sec_Rac = Sec_Rdc.*Sec_Fr;
    % Pcopper = (Iprms.^2.*Pri_Rac + Isrms.^2.*Sec_Rac);


% Calculate temperature rise
   
Rth = 16.31.*1e-3.*(Ac.*Wa).^(-0.405);
Tafterloss = Rth.*(Pcore) + 25;  %assume Pcopper << Pcore

% Copper Weight





