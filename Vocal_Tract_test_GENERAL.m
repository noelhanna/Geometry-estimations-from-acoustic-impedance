
function [ Z ] = ...
    Vocal_Tract_test_GENERAL(Airway, PARAMS, Flex)
%% Output Impedance Z from an airway model
% Airway is a structure containing:
%   VT, Glottis, SG
%   Each contain: radius, length, alpha_multiplier



%% define the calibration (radiation) impedance
% seen by the source for a single microphone ACUZ measurement

SourceTube.length = 10000;
SourceTube.radius = 0.003;
SourceTube.alpha_multiplier = 1;

% characteristic impedance
SourceTube.Z0(1, :) = z0(PARAMS, SourceTube.radius, 1);

% wavenumber with adjusted attenuation coefficent
SourceTube.K(1, :) = wavenumber(PARAMS, SourceTube, 1, SourceTube.radius);

% flange load impedance
Airway.Z_rad_source(1, :) = zradiation(SourceTube, 1, SourceTube.radius);

%% define flange load impedance at lips
% (the lip end is the highest numbered Airway segment)
lip_counter = size(Airway.VT.length,1);

% characteristic impedance
Airway.VT.Z0(lip_counter, :) = z0(PARAMS, Airway.VT.radius, lip_counter);

% wavenumber with adjusted attenuation coefficent
Airway.VT.K(lip_counter, :) = wavenumber(PARAMS, Airway.VT, (lip_counter), Airway.VT.radius);

% flange load impedance
Airway.VT.Z_rad(lip_counter, :) = zradiation(Airway.VT, lip_counter, Airway.VT.radius);

Airway.VT.Z_load(lip_counter, :) = Airway.VT.Z_rad(lip_counter, :);

% input impedance with flange termination
Airway.VT.Z(lip_counter, :) ...
    = zinput(Airway.VT.length, Airway.VT, Airway.VT.Z_load(lip_counter,:), lip_counter);

%%%%%%%%%%% define outgoing wave amplitude as 1
Airway.VT.A(lip_counter, :) = ones(1,size(PARAMS.freq,1));

% use a simple function for transfer functions TF = x/y
% use TF = transfer_fn(specified_out, specified_in)


%% General recursion for all other Airway segments
for k = lip_counter-1 : -1 : 1 % count backwards from inside lips
    % characteristic impedance
    Airway.VT.Z0(k,:) = z0(PARAMS, Airway.VT.radius, k);
    
    % wavenumber with adjusted attenuation coefficent
    Airway.VT.K(k,:) = wavenumber(PARAMS, Airway.VT, k, Airway.VT.radius);
    
    Airway.VT.Z_load(k, :) = Airway.VT.Z(k+1, :);
    
    % load impedance is the next Airway, i.e. Airway.VT.Z(k+1)
    Airway.VT.Z(k,:) = zinput(Airway.VT.length, Airway.VT, Airway.VT.Z(k+1,:), k);
    
    % assume outgoing wave amplitude is 1 - temporary
    Airway.VT.A(k, :) = ones(1,size(PARAMS.freq,1));
    
end

%% include non-rigid properties

Z_tissue_comp = zcompliance(PARAMS.w, Flex.C_tissue);
Z_tissue_inert = zinert(PARAMS.w, Flex.L_tissue);

% add the inertance and compliance in series
Ztissue = Z_tissue_inert + Z_tissue_comp;
% add the damping coefficient
Ztissue = Ztissue + Flex.R_tissue;

%%%%%%%%%%%% need to scale the amount added to the radius of the Airway %%
for o = 1:lip_counter
    % add non-rigid impedance to walls in parallel to Airway
    Airway.VT.Z_flex(o,:) = zparallel(Airway.VT.Z(o, :), Ztissue);
end



%% use lungs as Z_load (impedance in the other direction)
Airway.SG.Z0(1,:) = z0(PARAMS, Airway.SG.radius, 1);
Airway.SG.K(1,:) = wavenumber(PARAMS, Airway.SG, 1, Airway.SG.radius);

% flange load impedance
Airway.Z_lung = zradiation(Airway.SG, 1, Airway.SG.radius);

% SG tract
Airway.Z_SG(1,:) = zinput(Airway.SG.length, Airway.SG, Airway.Z_lung, 1);

% if more than 1 segment of SG tract is specified
if size(Airway.SG.radius,1) > 1
    for k = 2:1:size(Airway.SG.radius,1)-1 % count from lungs to glottis
        % characteristic impedance
        Airway.SG.Z0(k,:) = z0(PARAMS, Airway.SG.radius, k);
        
        % wavenumber with adjusted attenuation coefficent
        Airway.SG.K(k,:) = wavenumber(PARAMS, Airway.SG, k, Airway.SG.radius);
        
        % load impedance is the next Airway, i.e. Airway.Z(k-1)
        Airway.Z_SG(k,:) = zinput(Airway.SG.length, Airway.SG, Airway.Z_SG(k-1,:), k);
    end
end

Airway.Glottis.Z0 = z0(PARAMS, Airway.Glottis.radius, 1);
Airway.Glottis.K = wavenumber(PARAMS, Airway.Glottis, 1, Airway.Glottis.radius);

% -----------------------------------------
% TEMP 2017 use Flex params in parallel with SG
Z_SG_comp = zcompliance(PARAMS.w, Flex.C_tissue);
Z_SG_inert = zinert(PARAMS.w, Flex.L_tissue);

% add the inertance and compliance in series
ZSGflex = Z_SG_inert + Z_SG_comp;
% add the damping coefficient
ZSGflex = ZSGflex + Flex.R_tissue;
Airway.Z_SG_flex = zparallel(Airway.Z_SG(end, :), ZSGflex);
% -----------------------------------------

% glottis
Airway.Z_glottis = zinput(Airway.Glottis.length, Airway.Glottis, Airway.Z_SG(end,:), 1);

if Airway.Glottis.radius < 0.00001
    % use ideally stopped glottis
    Airway.Z_in(1, :) = zstop( Airway.VT, 1, Airway.VT.length);
else
    % use glottal impedance loaded by SG tract
    Airway.Z_in(1, :) = zinput(Airway.VT.length, Airway.VT, Airway.Z_glottis, 1);
end

%%%%%%%%%%%%% add non-rigidity to walls in parallel to Airway %%%%%%
%%%%%%% perhaps not same parameters%%
Airway.Z_in_flex(:,1) = zparallel(Airway.Z_in(1, :), Ztissue);


%% recursion for inward impedance
for m = 2:lip_counter
    % replace with load impedance of next Airway, i.e. Airway.Z(m-1)
    Airway.Z_in(m,:) ...
        = zinput(Airway.VT.length, Airway.VT, Airway.Z_in(m-1, :), m);
    
    %%%%%%%%%%%%%% check how this is calculated - scale with radius? %%%%%
    Airway.Z_in_flex(:,m) ...
        = zparallel(Airway.Z_in(m, :), Ztissue);
    
end


%% Z and T from lips
Airway.in.Z = Airway.Z_in(lip_counter,:);

example_C = 0.01*10^-8;
Z_example_comp  = zcompliance( PARAMS.w, example_C );
Airway.in.Z_with_compliance = zparallel(Airway.in.Z, Z_example_comp);

Airway.in.Z_flex = Airway.Z_in_flex(:,lip_counter);


% %ratios
Airway.in.Z_ratio_1 = (zparallel(Airway.VT.Z_rad(end,:), Airway.in.Z))...
    ./(Airway.VT.Z_rad(end,:));
Airway.in.Z_ratio_2 = (zparallel(Airway.VT.Z_rad(end,:), Airway.in.Z))...
    ./(Airway.Z_rad_source);
Airway.in.Z_ratio_3 = (zparallel(Airway.Z_rad_source, Airway.in.Z))...
    ./(Airway.Z_rad_source(end,:));
Airway.in.Z_ratio_4 = zparallel(Airway.Z_rad_source,...
    (zparallel(Airway.VT.Z_rad(end,:), Airway.in.Z)))...
    ./(Airway.VT.Z_rad(end,:));



%% Output data and calculate bandwidths
% % % formant = formant_bandwidth(PARAMS, Airway.T_pp_flex(1,:));
% % Airway.formant.f = formant.max.f;
% % Airway.formant.z = formant.max.z;
% % Airway.formant.b = formant.max.b;


% impedance = impedance_bandwidth(PARAMS, Airway.Z(1,:));
% Airway.impedance.max_f = impedance.max.f;
% Airway.impedance.max_z = impedance.max.z;
% Airway.impedance.max_b = impedance.max.b;
%
% Airway.impedance.min_f = impedance.min.f;
% Airway.impedance.min_z = impedance.min.z;
% Airway.impedance.min_b = impedance.min.b;


% calculate Airway area and radius as a function of length
% length_ = [Airway.Airway.SG.length; Airway.Glottis.length+sum(Airway.Airway.SG.length);...
%     sum(Airway.Glottis.length)+sum(Airway.SG.length)+Airway.length];
% %
% for g = 2:length(Airway.length)
%     length_(g+2) = Airway.length(g)+length_(g+1);
%     g
% end

length_SG(1) = Airway.SG.length(1);
for g = 2:size(Airway.SG.length,1)
    length_SG(g) = Airway.SG.length(g)+length_SG(g-1);
end

length_glottis(1) = length_SG(end)+Airway.Glottis.length(1);
for h = 2:size(Airway.Glottis.length,1)
    length_glottis(h) = Airway.Glottis.length(h)+length_glottis(h-1);
end

length_VT(1) = length_glottis(end)+Airway.VT.length(1);
for j = 2:size(Airway.VT.length,1)
    length_VT(j) = Airway.VT.length(j)+length_VT(j-1);
end

length_ = [length_SG'; length_glottis'; length_VT'];

radius_ = [Airway.SG.radius; Airway.Glottis.radius; Airway.VT.radius];

area_ = (radius_.^2).*pi;

Airway.area_fn = [length_; area_];
Airway.radius_fn = [length_; radius_];

% adjust profile to glottal 0
Airway.distance = Airway.radius_fn(:,1)-sum(Airway.SG.length)...
    -sum(Airway.Glottis.length);

% plot(Airway.area_fn(:,1), Airway.area_fn(:,2), 's-'); hold all


% output.freq = PARAMS.freq;

% Airway.Flex = Flex;

Z = abs(log10(Airway.in.Z_flex));

% output.PARAMS = PARAMS;




end









%% Background functions
% The functions below are called from the main code above



function soundSpeed = calculatespeedofsound(Parameters)
% calculatespeedofsound Calculates the speed of sound
% based on Owen Cramer (1993) J. Acoust. Soc. Am. 93(5)
% p2510-2616; formula at p. 2514
% SPEED_SOUND(t,h) calculates the speed of sound given
% the temperature in Celsius t and the relative humidity h
% (expressed as a fraction).
% The speed is calculated at one atm pressure and typical mole fraction
% of CO_2.
if (Parameters.humidity < 0) || (Parameters.humidity > 1),...
        error('The relative humidity must be between 0 and 1.'); end
% Use p = 1 atm and C0_2 mole fraction given by Cramer Table 1 (p. 2511)
PRESSURE_ATMOSPHERIC = 101325;
X_C = 0.000314;
% Calculate mole fraction of water vapour using equation in
% Cramer Appendix
% (p. 2515)
moleFractionOfWater = 1.00062 + 3.14e-8*PRESSURE_ATMOSPHERIC...
    + 5.6e-7 * Parameters.temperature^2;
temperatureKelvin = Parameters.temperature + 273.15;
p_sv = exp(1.2811805e-5*temperatureKelvin^2 - 1.9509874e-2...
    * temperatureKelvin + 34.04926034 - 6.3536311e3 / temperatureKelvin);
x_w = Parameters.humidity * moleFractionOfWater ...
    * p_sv/PRESSURE_ATMOSPHERIC;
soundSpeed = calculatespeedofsoundcramer(Parameters,...
    PRESSURE_ATMOSPHERIC,x_w,X_C);
end

function soundSpeed = calculatespeedofsoundcramer(Parameters,...
    PRESSURE_ATMOSPHERIC,x_w,X_C)
% SPEED_SOUND_CRAMER Calculate the speed of sound.
% Based on Owen Cramer (1993) J. Acoust. Soc. Am. 93(5) p2510-2616;
% formula at p2514
% SPEED_SOUND_CRAMER(t,p,x_w,x_c) calculates the speed of sound given the
% temperature in Celsius t, the pressure in Pa p and the mole fraction of
% water vapour and CO_2 x_w and x_c.
if (Parameters.temperature < 0) || (Parameters.temperature > 30),...
        error('The temperature must be between 0 and 30 degrees C.'); end
if (PRESSURE_ATMOSPHERIC < 75e3) || (PRESSURE_ATMOSPHERIC > 102e3),...
        error ('The pressure must be between 75 and 102 kPa.'); end
if (x_w < 0) || (x_w > 0.06),...
        error ('The H2O mole fraction must be between 0 and 0.06.'); end
if (X_C < 0) || (X_C > 0.01),...
        error ('The CO2 mole fraction must be between 0 and 0.01.'); end
a = [
    331.5024
    0.603055
    -0.000528
    51.471935
    0.1495874
    -0.000782
    -1.82e-7
    3.73e-8
    -2.93e-10
    -85.20931
    -0.228525
    5.91e-5
    -2.835149
    -2.15e-13
    29.179762
    0.000486
    ]';
coeff = [
    1
    Parameters.temperature
    Parameters.temperature^2
    x_w
    Parameters.temperature*x_w
    Parameters.temperature^2*x_w
    PRESSURE_ATMOSPHERIC
    Parameters.temperature*PRESSURE_ATMOSPHERIC
    Parameters.temperature^2*PRESSURE_ATMOSPHERIC
    X_C
    Parameters.temperature*X_C
    Parameters.temperature^2*X_C
    x_w^2
    PRESSURE_ATMOSPHERIC^2
    X_C^2
    x_w*PRESSURE_ATMOSPHERIC*X_C
    ];
soundSpeed = a*coeff;
end


%% characteristic impedance of a cylinder of specified radius
function [ Z0 ] = z0(PARAMS, RAD, counter)
rho = PARAMS.rho;
c = PARAMS.c;
radius = RAD(counter);

Z0 = rho*c/(pi.*(radius.^2));
end

%% impedance of an inertance with length l, area A,
function [ Zinert ] = zinert( w, L )
Zinert = 1i.*w.*L;

end

%% impedance of a known compliance C = Volume/rho*c^2
function [ Zcomp ] = zcompliance( w, C )
% Zcompliance = 1/jwC,
Zcomp = 1./(1i.*w.*C);

end


%% complex wavenumber including wall losses
function [ K ] = wavenumber(PARAMS, tube, counter, RAD)

freq = PARAMS.freq;
w = PARAMS.w;
c = PARAMS.c;
DeltaT = PARAMS.DeltaT;

radius = RAD(counter);
alpha_multiplier = tube.alpha_multiplier(counter);

rv = 632.8...
    .*radius...
    .*(freq.^0.5)...
    .*(1-0.0029...
    .*(DeltaT)); % approx

rt = 532.2...
    .*radius...
    .*(freq.^0.5)...
    .*(1-0.0031...
    .*(DeltaT)); % approx

% Cp./Cv ratio of specific heats (approx 1.403 at 0, 1.401 at 100)
gamma = 1.400;

v = c...
    .*( 1 - (1./(rv.*sqrt(2)))...
    -((gamma-1)./(rt.*sqrt(2))) );

% attenuation coefficient alpha
alpha = (w./c)...
    .*( (1./(rv.*sqrt(2))) + ((gamma-1)./(rt.*sqrt(2))) );

% adjust alpha to approximate real VT measurements (alpha_multiplier 5)
new_alpha = alpha.*alpha_multiplier;

% complex propogation coefficient
K = ((w./v) - (1i.*new_alpha));

end

%% flanged opening radiation impedance
function [ Z_flange ] = zradiation(tube, counter, RAD)
% RADIATION IMPEDANCE CALCULATION FROM DALMONT ET AL. 2001
% calculates the radiation impedance with a flange

K = tube.K(counter,:);
radius = RAD(counter);
Z0 = tube.Z0(counter);

ka = K*radius;

% % % % simple approximation from Fletcher and Rossing
% % % % real part of radiation impedance (F&R eq 8.29)
% % % Z_flange_R = Z0.*( (((ka).^2)./2) - (((ka).^4)./(2.*2.*3))...
% % %     + (((ka).^6)./(2.*2.*3.*3.*4)) );
% % % % imaginary part of radiation impedance (F&R eq 8.30)
% % % Z_flange_X = (Z0./(pi.*ka.^2)).*( (((2.*ka).^3)./3)...
% % %     - (((2.*ka).^5)./(3.*3.*5)) ...
% % %     + (((2.*ka).^7)./(3.*3.*5.*5.*7)) );
% % %
% % % Z_flange = Z_flange_R + 1i.*Z_flange_X;


% alternative approximation from Dalmont et al. (2001)
% define the end correction for the low frequency limit
d_simple = 0.8216 * radius;

% determine the frequency-dependent end correction (15a)
d_simple = d_simple ./ (1 + (0.77 * ka).^2 ./ (1 + 0.77 * ka));

% determine the modulus of the reflection coefficient (15b)
modR = (1 + (0.323 .* ka) - (0.077 .* ka.^2))...
    ./ (1 + (0.323 .* ka) + ((1 - 0.077) .* ka.^2));

% calculate the imaginary part of the end correction
% since R = -e^(-2kjd(complex)) = -modR*e^(-2kjd(real))
% so, d(complex) = d*ln(modR))
di = log(modR)./(2 * K);
d = d_simple + 1i.*di;

% calculate the impedance (9)
Z_flange = 1i * Z0 * tan(K .* d); % Flanged opening

end

%% Input impedance of a cylinder with known Z0, Z_load, K and L
function [ Zin ] = zinput(LGR, tube, Z_load, counter)

Z0 = tube.Z0(counter);
K = tube.K(counter,:);
L = LGR(counter);

% from Fletcher and Rossing Chapter 8 (8.23)
Zin = Z0...
    .*( ( (Z_load.*cos(K.*L)) + (1i.*Z0.*sin(K.*L)) ) ...
    ./ ( (1i.*Z_load.*sin(K.*L)) + (Z0.*cos(K.*L)) ) );
end

%% impedance of an ideally stopped end Z = -jZocotkL
function [ Zstop ] = zstop( tube, counter, LGR)

K = tube.K(counter,:);
L = LGR(counter);
Z0 = tube.Z0(counter);

% from Fletcher and Rossing Chapter 8 (8.24)
% Zstop = -jZocotkL
Zstop = -1i.*Z0.*cot(K.*L);

end

%% Parallel addition
function [ Zparallel ] = zparallel(Z1, Z2)
% parallel addition of impedances
Zparallel = (Z1.*Z2)./(Z1 + Z2);
end


%% determine input and output pressure and flow at each segment
% function [ p_out, p_in, u_out, u_in, BonA ] = ...
%     calc_pressure_flow(tube, counter, LGR)
% % from J Smith via Fletcher and Rossing
%
% K = tube.K(counter,:);
% L = LGR(counter);
% Z_load = tube.Z_load(counter,:);
% Z0 = tube.Z0(counter);
% A = tube.A(counter);
%
% % exponential term
% eTerm = exp(1i.*K.*L);
%
% % definitions p = A/eTerm + BeTerm, u = 1/Z0 (A/eTerm - BeTerm)
% % reflection coefficient is B/A
% BonA = (Z_load - Z0)...
%     ./((Z_load + Z0).*eTerm.*eTerm);
%
% % pressure calculations
% p_out = A.*((1./eTerm) + (BonA.*eTerm));
% p_in = A.*(1 + BonA);
% % flow calculations
% u_out = A.*(1./Z0).*((1./eTerm) - (BonA.*eTerm));
% u_in = A.*(1./Z0).*(1 - BonA);

% end



%% determine amplitude adjusted pressure and flow at each segment
function [A_adjusted] = adjust_amplitude(A, BonA)
% The first segment sends an outgoing wave with amplitude A = 1;
% the reflection coefficent R = B/A means B/A is reflected at each segment
% so the next segment has an input wave amplitude A* = 1-BonA
A_adjusted = A.*BonA;

end

%% General transfer function calculation %%%%%%%%% not yet used %%%%%%%%%%
function [TF] = transfer_fn(specified_output, specified_input)
% care should be taken to choose appropriate real/abs values
TF = specified_output/specified_input;

end

%% BANDWIDTH MEASUREMENTS (Optional)
function [output] = formant_bandwidth(PARAMS, Transfer)

freq = PARAMS.freq;
resolution = PARAMS.resolution;

% need to specify upper and lower bounds of the frequency
% for each pair of max and min
max_first = [150 1050 2250 3250];
max_last = [1000 2200 3200 4000];
% create vectors of the right size to be filled below
z_max = [nan nan nan nan];
f_max = [nan nan nan nan];
max_bandwidth = [nan nan nan nan];

for n = 1:4
    % convert to index of the frequency vector
    [~, c2] =find(freq>=(max_first(n)));
    low2(n) = min(c2); % max_first index above the extreme
    [~, c2] =find(freq<=(max_last(n)));
    high2(n) = max(c2); % max_first index above the extreme
    
    [value2, index2] = max(abs(Transfer(low2(n):high2(n))));
    % zmax is the maximum impedance value
    z_max(n) = value2;
    % fmax is the frequency of the minimum impedance
    f_max(n) = max_first(n) - resolution + (index2.*resolution);
    [~, max_index] = find(abs(Transfer(low2(n):high2(n)))...
        >=(max(abs(Transfer(low2(n):high2(n))))./sqrt(2)));
    f_max_range = max_first(n) - resolution + (max_index.*resolution);
    if max(f_max_range) >= 1
        max_bandwidth(n) = max(f_max_range) - min(f_max_range);
    end
    
    % ignore max_bandwidth greater than 250
    if max_bandwidth(n) > 250;
        max_bandwidth(n) = nan;
    end
    
end

% OUTPUT
output.T = Transfer;
output.max.z(:,1) = z_max;
output.max.f(:,1) = f_max;
output.max.b(:,1) = max_bandwidth;

end



%% BANDWIDTH MEASUREMENTS (Optional)
function [output] = impedance_bandwidth(PARAMS, Transfer)

freq = PARAMS.freq;
resolution = PARAMS.resolution;

% need to specify upper and lower bounds of the frequency
% for each pair of max and min
max_first = [150 700 1800 2600 3600];
max_last = [400 1500 2500 3500 4200];
% create vectors of the right size to be filled below
z_max = [nan nan nan nan nan];
f_max = [nan nan nan nan nan];
max_bandwidth = [nan nan nan nan nan];

% min_first = [10 200 1000 2000 3100];
% min_last = [100 800 2000 3000 4000];

min_first = [10 400 1100 2000 3000 ];
min_last = [100 800 1800 2600 3600 ];

% create vectors of the right size to be filled below
z_min = [nan nan nan nan nan];
f_min = [nan nan nan nan nan];
min_bandwidth = [nan nan nan nan nan];

for n = 2:5 %%%%%%%%%%%%%% change this to look at other resonances
    % convert to index of the frequency vector
    [~, c2] =find(freq>=(max_first(n)));
    low2(n) = min(c2); % max_first index above the extreme
    [~, c2] =find(freq<=(max_last(n)));
    high2(n) = max(c2); % max_first index above the extreme
    
    [value2, index2] = max(abs(Transfer(low2(n):high2(n))));
    % zmax is the maximum impedance value
    z_max(n) = value2;
    % fmax is the frequency of the minimum impedance
    f_max(n) = max_first(n) - resolution + (index2.*resolution);
    [~, max_index] = find(abs(Transfer(low2(n):high2(n)))...
        >=(max(abs(Transfer(low2(n):high2(n))))./sqrt(2)));
    f_max_range = max_first(n) - resolution + (max_index.*resolution);
    if max(f_max_range) >= 1
        max_bandwidth(n) = max(f_max_range) - min(f_max_range);
    end
    
    % % %     % limit max_bandwidth to 250
    % % %     if max_bandwidth(n) > 250;
    % % %         max_bandwidth(n) = nan;
    % % %     end
    
    
    % convert to index of the frequency vector
    [~, c1] =find(freq>=(min_first(n)));
    low1(n) = min(c1); % min_first index above the extreme
    [~, c1] =find(freq<=(min_last(n)));
    high1(n) = max(c1); % min_first index above the extreme
    
    
    [value1, index1] = min(abs(Transfer(low1(n):high1(n))));
    % zmin is the minimum impedance value
    z_min(n) = value1;
    % fmin is the frequency of the minimum impedance
    f_min(n) = min_first(n) - resolution + (index1.*resolution);
    [~, min_index] = find(abs(Transfer(low1(n):high1(n)))...
        <=(min(abs(Transfer(low1(n):high1(n)))).*sqrt(2)));
    f_min_range = min_first(n) - resolution + (min_index.*resolution);
    if min(f_min_range) >= 1
        min_bandwidth(n) = max(f_min_range) - min(f_min_range);
    end
    
    % % %     % limit the min_bandwidth to 250
    % % %     if min_bandwidth(n) > 250;
    % % %         min_bandwidth(n) = nan;
    % % %     end
    
end

% OUTPUT
% Choose which variables to save to the output
output.T = Transfer;
output.max.z(:,1) = z_max;
output.max.f(:,1) = f_max;
output.max.b(:,1) = max_bandwidth;

output.min.z(:,1) = z_min;
output.min.f(:,1) = f_min;
output.min.b(:,1) = min_bandwidth;

end



