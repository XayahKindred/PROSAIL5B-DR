%*************************************************************************
%                             main_PROSAIL
% This program simulates vegetation canopy reflectance based on the PROSAIL model:
% - Leaf optical properties: PROSPECT-5 (Feret et al., 2008);
% - Leaf inclination distribution:
%   · Ellipsoidal distribution (Campbell, 1990);
%   · Two-parameter LIDF distribution (Verhoef, 1998);
% - Canopy reflectance: 4SAIL model (Verhoef et al., 2007).
%
% Our contribution:
% - Introduced a new parameter, dust retention (Cdust, g/m²), into PROSAIL 
%   for simulating the spectral impact of leaf surface dust deposition;
% - The extended PROSAIL-DR model accurately captures the spectral features 
%   of urban vegetation, particularly suited for studying relationships between 
%   dust deposition and canopy reflectance.
%*************************************************************************

clc
TypeLidf=2;
% if 2-parameters LIDF: TypeLidf=1
if (TypeLidf==1)
    % LIDFa LIDF parameter a, which controls the average leaf slope
    % LIDFb LIDF parameter b, which controls the distribution's bimodality
    %	LIDF type 		a 		 b
    %	Planophile 		1		 0
    %	Erectophile    -1	 	 0
    %	Plagiophile 	0		-1
    %	Extremophile 	0		 1
    %	Spherical 	   -0.35 	-0.15
    %	Uniform 0 0
    % 	requirement: |LIDFa| + |LIDFb| < 1	
    LIDFa	=	-0.35;
    LIDFb	=	-0.15;
% if ellipsoidal LIDF: TypeLidf=2
elseif (TypeLidf==2)
    % 	LIDFa	= average leaf angle (degrees) 0 = planophile	/	90 = erectophile
    % 	LIDFb = 0
    LIDFa	=	30;
    LIDFb	=	0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEAF CHEM & STR PROPERTIES	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cab    = 31;    % Chlorophyll content (µg·cm⁻²)
Car    = 5.8;     % Carotenoid content (µg·cm⁻²)
Cbrown = 0.0;   % Brown pigment content (arbitrary units)
Cw     = 0.01;  % Equivalent water thickness, EWT (cm)
Cm     = 0.01; % Leaf mass per area, LMA (g·cm⁻²)
N      = 2.5;   % Leaf structure coefficient (dimensionless)
Cdust  = 0;    % Dust retention per unit leaf area (g·m⁻²)

data = dataSpec_DR();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Soil Reflectance Properties	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsoil1 = dry soil
% rsoil2 = wet soil
Rsoil1  = data(:,10);Rsoil2=data(:,11);
psoil	= 1;		% soil factor (psoil=0: wet soil / psoil=1: dry soil)
rsoil0  = psoil*Rsoil1+(1-psoil)*Rsoil2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	4SAIL canopy structure parm	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAI   = 5.;    % Leaf Area Index (m²·m⁻²)
hspot = 0.01;  % Hot spot parameter (dimensionless)
tts   = 30.;   % Solar zenith angle (degrees)
tto   = 0.;    % Observer zenith angle (degrees)
psi   = 90.;   % Relative azimuth angle (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          CALL PRO4SAIL       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rdot: hemispherical-directional reflectance factor in viewing direction    
% rsot: bi-directional reflectance factor
% rsdt: directional-hemispherical reflectance factor for solar incident flux
% rddt: bi-hemispherical reflectance factor
[rdot,rsot,rddt,rsdt]=PRO4SAIL(N,Cab,Car,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil0,Cdust);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	direct / diffuse light	%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the direct and diffuse light are taken into account as proposed by:
% Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo 
% measurements into total shortwave broadband albedo using a canopy 
% radiative transfer model, Agronomie
% Es = direct
% Ed = diffuse
Es  = data(:,8);
Ed  = data(:,9);
rd  = pi/180;
skyl	=	0.847- 1.61*sin((90-tts)*rd)+ 1.04*sin((90-tts)*rd)*sin((90-tts)*rd); % % diffuse radiation
PARdiro	=	(1-skyl)*Es;
PARdifo	=	(skyl)*Ed;

% resv: Directional canopy reflectance calculation
resv = (rdot .* PARdifo + rsot .* PARdiro) ./ (PARdiro + PARdifo);

% Save reflectance data to file
dlmwrite('Refl_CAN.txt', [data(:,1), resv], 'delimiter', '\t', 'precision', 5);

% Extract wavelength data
wavelength = data(:,1);

% Plot the entire canopy reflectance spectrum
figure;
plot(wavelength, resv, 'k', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Reflectance');
grid on;
axis([400 2500 0 0.6]);

% Adjust title and legend based on dust presence
if Cdust > 0
    title('Canopy Reflectance with Dust Effect');
    legend(sprintf('C_{dust} = %.2f g·m^{-2}', Cdust), 'Location', 'Best');
else
    title('Canopy Reflectance without Dust Effect');
end


