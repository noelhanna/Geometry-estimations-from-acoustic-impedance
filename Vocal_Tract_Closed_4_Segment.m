function [ Z ] = ...
    Vocal_Tract_Closed_4_Segment(x, galpha, gl, gr,...
    l1,l2,l3,l4,...
    r1,r2,r3,r4,...
    vt_alpha, w1, w2)

%% Define default parameters (if none specified)

% Vocal Tract
VT.length = [l1 l2 l3 l4]';
VT.radius = [r1 r2 r3 0.013]';
for i = 1:size(VT.length)
    VT.alpha_multiplier(i)= vt_alpha;
end

% Subglottal Tract
SG.length = [0.2]';
SG.radius = [0.01]';
SG.alpha_multiplier = [vt_alpha]';


% compliance of tissue
Flex.C_tissue = 4.9*10^-8; %2.9*10^-8;
% inertance of tissue
Flex.L_tissue = w1;
% resistance of tissue
Flex.R_tissue = w2;

% glottis
glottis.alpha_multiplier = galpha;
glottis.radius = gr;
glottis.length = gl;

% General Parameters
PARAMS.Temp = 20.2;
PARAMS.humidity = 0.4;
% used for the visco-thermal loss calculations
PARAMS.DeltaT = PARAMS.Temp-26.85;
PARAMS.resolution = 5.3833; % frequency resolution in Hz
PARAMS.freq = x';
PARAMS.w = 2.*pi.*PARAMS.freq; % angular frequency
PARAMS.c = 340; % standard
PARAMS.rho = 1.14; % air density


Airway.VT = VT;
Airway.Glottis = glottis;
Airway.SG = SG;

% pass all parameters into appropriate impedance calculation

Z = Vocal_Tract_test_GENERAL(Airway, PARAMS, Flex);

end