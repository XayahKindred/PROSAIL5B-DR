% _______________________________________________________________________
%
% prospect_DR.m (Extended PROSPECT model considering Dust Retention)
% Subroutines required: tav.m, dataSpec_DR.m
% _______________________________________________________________________
% for any question or request, please contact: 
% _______________________________________________________________________
% 
% This program calculates plant leaf reflectance and transmittance spectra
% from 400 nm to 2500 nm (1 nm step) using an extended PROSPECT model
% (PROSPECT-DR), which includes the new dust retention parameter (Cdust).
%
% % Model parameters:
%
%   - N       : Leaf structure parameter (dimensionless)
%   - Cab     : Chlorophyll a+b content (µg/cm²)
%   - Car     : Carotenoid content (µg/cm²)
%   - Cbrown  : Brown pigment content (arbitrary units)
%   - Cw      : Equivalent water thickness (cm or g/cm²)
%   - Cm      : Dry matter content (g/cm²)
%   - Cdust  : Dust retention per unit leaf area (g/m²)
%
% Examples of typical dust retention values (Cdust) for common urban plants, 
% based on field samples collected in Shanghai, China, are provided for reference.
% ---------------------------------------------------
%   Plant species                Cdust (g/m²)
% ---------------------------------------------------
%   Osmanthus fragrans             1.85–16.21
%   Erythrina coralloides          2.50–30.53
%   Ligustrum × vicaryi            1.29–31.33
%   Photinia × fraseri             1.84–29.81
%   Buxus sempervirens 'Variegata' 2.52–42.05
%   Cinnamomum camphora            1.07–14.94
%   Photinia serratifolia          1.14–35.06
%   Rhododendron spp.              1.73–16.85
%   Pittosporum tobira             1.48–12.92
%   Camellia japonica              2.57–6.73
% ---------------------------------------------------
% _______________________________________________________________________

function LRT=prospect_DR(N,Cab,Car,Cbrown,Cw,Cm,Cdust)

% ***********************************************************************
% Jacquemoud S., Baret F. (1990), PROSPECT: a model of leaf optical
% properties spectra, Remote Sens. Environ., 34:75-91.
% Feret et al. (2008), PROSPECT-4 and 5: Advances in the Leaf Optical 
% Properties Model Separating Photosynthetic Pigments, Remote Sensing of
% Environment, 112:3030-3043
% The specific absorption coefficient corresponding to brown pigment is
% provided by Frederic Baret (EMMAH, INRA Avignon, baret@avignon.inra.fr)
% and used with his autorization.
% ***********************************************************************

data=dataSpec_DR();
l=data(:,1);
n=data(:,2);
k_dust = data(:, end);  
k=(Cab*data(:,3)+Car*data(:,4)+Cbrown*data(:,5)+Cw*data(:,6)+Cm*data(:,7)+Cdust*k_dust)./(N);
k(find(k==0))=eps;
trans=(1-k).*exp(-k)+k.^2.*expint(k);

% ***********************************************************************
% reflectance and transmittance of one layer
% ***********************************************************************
% Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
% Interaction of isotropic ligth with a compact plant leaf, J. Opt.
% Soc. Am., 59(10):1376-1379.
% ***********************************************************************

% reflectivity and transmissivity at the interface
%-------------------------------------------------
alpha=40;
t12=tav(alpha,n);
t21=tav(90,n)./n.^2;
r12=1-t12;
r21=1-t21;
x=tav(alpha,n)./tav(90,n);
y=x.*(tav(90,n)-1)+1-tav(alpha,n);

% reflectance and transmittance of the elementary layer N = 1
%------------------------------------------------------------

ra=r12+(t12.*t21.*r21.*trans.^2)./(1-r21.^2.*trans.^2);
ta=(t12.*t21.*trans)./(1-r21.^2.*trans.^2);
r90=(ra-y)./x;
t90=ta./x;

% ***********************************************************************
% reflectance and transmittance of N layers
% ***********************************************************************
% Stokes G.G. (1862), On the intensity of the light reflected from
% or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
% 11:545-556.
% ***********************************************************************
delta=((t90.^2-r90.^2-1).^2-4*r90.^2);
beta=(1+r90.^2-t90.^2-sqrt(delta))./(2*r90);
va=(1+r90.^2-t90.^2+sqrt(delta))./(2*r90);
vb=sqrt(beta.*(va-r90)./(va.*(beta-r90)));
vbNN = vb.^(N-1);
vbNNinv = 1./vbNN;
vainv = 1./va;
s1=ta.*t90.*(vbNN-vbNNinv);
s2=ta.*(va-vainv);
s3=va.*vbNN-vainv.*vbNNinv-r90.*(vbNN-vbNNinv);

RN=ra+s1./s3;
TN=s2./s3;
LRT=[l RN TN];