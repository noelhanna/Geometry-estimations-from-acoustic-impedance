function fit_general_cylinder_model_v0_3(handles, VT, Glottis, SG, Rigidity)

%% Load impedance data to fit

load(handles.LoadImpedanceFile.String);
% % load('Closed_In_Ex.mat')
Z_measured = Iteration.Z(:, handles.IterationSelection.Value);
x = Parameters.frequencyVector;


%% PARAMETERS INITIALIZATION
% % the last radius (at the mouth) is set at 0.013 in the function

radius = VT.radius';

length_VT = VT.length';

length_SG = SG.length';

radius_SG = SG.radius';

radius_glottis = Glottis.radius;
gr = radius_glottis;

length_glottis = Glottis.length;
gl = length_glottis;

alpha = Rigidity.alpha_multiplier;

L = Rigidity.inertanceInitial;
R = 5*10^5;

% attenuation coefficients
galpha = Rigidity.alpha_multiplier;
vt_alpha = Rigidity.alpha_multiplier;
subgalpha = Rigidity.alpha_multiplier;

%% -----------------------------
% set fitting function with the correct number of free parameters
if Glottis.state == 1 % closed glottis
    switch max(size(radius)) % no. of cylindrical VT segments
        case    1
            ft = fittype( 'Vocal_Tract_Closed_1_Segment(x, galpha, gl, gr, l1, r1, vt_alpha, w1, w2)' );
        case    2
            ft = fittype( 'Vocal_Tract_Closed_2_Segment(x, galpha, gl, gr, l1, l2, r1, r2, vt_alpha, w1, w2)' );
        case    3
            ft = fittype( 'Vocal_Tract_Closed_3_Segment(x, galpha, gl, gr, l1, l2, l3, r1, r2, r3, vt_alpha, w1, w2)' );
        case    4
            ft = fittype( 'Vocal_Tract_Closed_4_Segment(x, galpha, gl, gr, l1, l2, l3, l4, r1, r2, r3, r4, vt_alpha, w1, w2)' );
        case    5
            ft = fittype( 'Vocal_Tract_Closed_5_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, r1, r2, r3, r4, r5, vt_alpha, w1, w2)' );
        case    6
            ft = fittype( 'Vocal_Tract_Closed_6_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, r1, r2, r3, r4, r5, r6, vt_alpha, w1, w2)' );
        case    7
            ft = fittype( 'Vocal_Tract_Closed_7_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, r1, r2, r3, r4, r5, r6, r7, vt_alpha, w1, w2)' );
        case    8
            ft = fittype( 'Vocal_Tract_Closed_8_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, l8, r1, r2, r3, r4, r5, r6, r7, r8, vt_alpha, w1, w2)' );
        case    9
            ft = fittype( 'Vocal_Tract_Closed_9_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, l8, l9, r1, r2, r3, r4, r5, r6, r7, r8, r9, vt_alpha, w1, w2)' );
        case 10
            ft = fittype( 'Vocal_Tract_Closed_10_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, vt_alpha, w1, w2)' );
    end
    
else % open glottis
    switch  max(size(radius_SG)) % no. of cylindrical SG segments
        case 1
            ft = fittype( 'Vocal_Tract_Open_1_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, r1, r2, r3, r4, r5, r6, r7, sgl1, sgr1, subgalpha, vt_alpha, w1, w2)' );
            
        case 2
            ft = fittype( 'Vocal_Tract_Open_2_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, r1, r2, r3, r4, r5, r6, r7, sgl1, sgl2, sgr1, sgr2, subgalpha, vt_alpha, w1, w2)' );
            
        case 3
            ft = fittype( 'Vocal_Tract_Open_3_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, r1, r2, r3, r4, r5, r6, r7, sgl1, sgl2, sgl3, sgr1, sgr2, sgr3, subgalpha, vt_alpha, w1, w2)' );
            
        case 4
            ft = fittype( 'Vocal_Tract_Open_4_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, r1, r2, r3, r4, r5, r6, r7, sgl1, sgl2, sgl3, sgl4, sgr1, sgr2, sgr3, sgr4, subgalpha, vt_alpha, w1, w2)' );

%         case 10
%             %         ft = fittype( 'Vocal_Tract_test_OPENGLOTTIS_v0_2(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, sgl1, sgl2, sgl3, sgl4, sgl5, sgl6, sgl7, sgl8, sgl9, sgl10, sgr1, sgr2, sgr3, sgr4, sgr5, sgr6, sgr7, sgr8, sgr9, sgr10, subgalpha, vt_alpha, w1, w2)' );
%             ft = fittype( 'Vocal_Tract_Open_10_Segment(x, galpha, gl, gr, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, sgl1, sgl2, sgl3, sgl4, sgl5, sgl6, sgl7, sgl8, sgl9, sgl10, sgr1, sgr2, sgr3, sgr4, sgr5, sgr6, sgr7, sgr8, sgr9, sgr10, subgalpha, vt_alpha, w1, w2)' );
            
    end
end


%% FITOPTIONS SETTING
fo = fitoptions(ft);

if Glottis.state == 1 % closed glottis
    fo.StartPoint = [alpha length_glottis radius_glottis length_VT ...
        radius alpha L R];
else % open glottis
    fo.StartPoint = [alpha length_glottis radius_glottis length_VT ...
        radius length_SG radius_SG alpha alpha L R];
end

% % upper and lower limits of the variables
% glottal alpha
fo.Upper(1) = Rigidity.alphaLimits(2);
fo.Lower(1) = Rigidity.alphaLimits(1);

% glottal length
fo.Upper(2) = Glottis.lengthLimits(2);
fo.Lower(2) = Glottis.lengthLimits(1);

% glottal radius
fo.Upper(3) = Glottis.radiusLimits(2);
fo.Lower(3) = Glottis.radiusLimits(1);

% vocal tract length and radius
for i=4:3+(2*((max(size(radius)))))
    if i<=(3+max(size(radius)))
        fo.Upper(i) = VT.lengthLimits(2);
        fo.Lower(i) = VT.lengthLimits(1);
    else
        fo.Upper(i) = VT.radiusLimits(2);
        fo.Lower(i) = VT.radiusLimits(1);
    end
end

% lip radius fixed to match large impedance head
fo.Upper(3+(2*(max(size(radius))))) = 0.0131;
fo.Lower(3+(2*(max(size(radius))))) = 0.0129;

if Glottis.state == 1
    % alpha multiplier
    fo.Upper(3+(2*max(size(radius)))+1) = Rigidity.alphaLimits(2);
    fo.Lower(3+(2*max(size(radius)))+1) = Rigidity.alphaLimits(1);
    % L tissue
    fo.Upper(3+(2*max(size(radius)))+2) = Rigidity.inertanceLimits(2);
    fo.Lower(3+(2*max(size(radius)))+2) = Rigidity.inertanceLimits(1);
    % R tissue
    fo.Upper(3+(2*max(size(radius)))+3) = 10*10^5;
    fo.Lower(3+(2*max(size(radius)))+3) = 1*10^5;
    
else % include subglottal tract
    % subglottal tract length and radius
    for i = (4+(2*max(size(radius)))):...
            (4+(2*(max(size(radius)))))+((max(size(radius_SG)))-1)
        fo.Upper(i) = SG.lengthLimits(2);
        fo.Lower(i) = SG.lengthLimits(1);
    end
    for i = (4+(2*(max(size(radius)))))+((max(size(radius_SG)))) : ...
            (4+(2*(max(size(radius)))))+((max(size(radius_SG))))...
            +(max(size(radius_SG)))-1
        fo.Upper(i) = SG.radiusLimits(2);
        fo.Lower(i) = SG.radiusLimits(1);
    end
    % alpha multiplier
    fo.Upper((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+1) = Rigidity.alphaLimits(2);
    fo.Lower((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+1) = Rigidity.alphaLimits(1);
    fo.Upper((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+2) = Rigidity.alphaLimits(2);
    fo.Lower((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+2) = Rigidity.alphaLimits(1);
    % L tissue
    fo.Upper((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+3) = Rigidity.inertanceLimits(2);
    fo.Lower((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+3) = Rigidity.inertanceLimits(1);
    % R tissue
    fo.Upper((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+4) = 10*10^5;
    fo.Lower((3+(2*(max(size(radius)))))+((max(size(radius_SG))))...
        +(max(size(radius_SG)))+4) = 1*10^5;
end

%% Perfrom fit to find data arguments(f) and determine goodness
axes(handles.axes2);
[f, gof, fit_output] = fit( x, abs(log10(Z_measured)), ft, fo);
[f_isnormal_chi2, f_P_isnormal] = chi2gof(fit_output.residuals);
% f = fit( x, opg_imp3, ft, fo);

%% Update structures VT, SG, Glottis etc. with fitted data
if Glottis.state == 1 % closed glottis
    
    switch max(size(radius))
        case 1
            TotalAirway.radius = [f.gr 0.013];  % total radius vector
            VT.radius = TotalAirway.radius(2:end);
            TotalAirway.length = [f.gl f.l1 ];   % total length vector
            VT.length = TotalAirway.length(2:end);
            VT.alpha_multiplier = f.vt_alpha;
            VT.radius = TotalAirway.radius(2:end);
            SG.totalLength = 0;
            VT.totalLength = f.l1;  % vocal tract length
            len(1) = f.gl;
            for i =2:2
                len(i)=len(i-1)+TotalAirway.length(i);
                VT.alpha_multiplier(i) = f.vt_alpha;
            end
            
        case 2
            TotalAirway.radius = [f.gr f.r1 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2;  % vocal tract length
            len(1) = f.gl;
            for i =2:3
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 3
            TotalAirway.radius = [f.gr f.r1 f.r2 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3;  % vocal tract length
            len(1) = f.gl;
            for i =2:4
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 4
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4  ];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4;  % vocal tract length
            len(1) = f.gl;
            for i =2:5
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 5
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 ];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5;  % vocal tract length
            len(1) = f.gl;
            for i =2:6
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 6
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 f.r5 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 ];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6;  % vocal tract length
            len(1) = f.gl;
            for i =2:7
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 7
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7 ];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7;  % vocal tract length
            len(1) = f.gl;
            for i =2:8
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 8
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 f.r7 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7 f.l8];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7+ f.l8;  % vocal tract length
            len(1) = f.gl;
            for i =2:9
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 9
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 f.r7 f.r8 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7 f.l8 f.l9 ];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7+ f.l8+ f.l9;  % vocal tract length
            len(1) = f.gl;
            for i =2:10
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
            
        case 10
            TotalAirway.radius = [f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 f.r7 f.r8 f.r9 0.013];  % total radius vector
            TotalAirway.length = [f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7 f.l8 f.l9 f.l10];   % total length vector
            SG.totalLength = 0;
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7+ f.l8+ f.l9+f.l10;  % vocal tract length
            len(1) = f.gl;
            for i =2:11
                len(i)=len(i-1)+TotalAirway.length(i);
            end
    end
    
else % open glottis
    
    switch max(size(radius_SG))
        case 1
            TotalAirway.radius = [f.sgr1 f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 0.013];  % total radius vector
            TotalAirway.length = [f.sgl1 f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7];   % total length vector
            SG.totalLength = f.sgl1; % subglottal tract length
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7;  % vocal tract length
            len(1) = f.sgl1;
            for i =2:9
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 2
            TotalAirway.radius = [f.sgr2 f.sgr1 f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 0.013];  % total radius vector
            TotalAirway.length = [f.sgl2 f.sgl1 f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7];   % total length vector
            SG.totalLength = f.sgl1 +f.sgl2; % subglottal tract length
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7;  % vocal tract length
            len(1) = f.sgl2;
            for i =2:10
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 3
            TotalAirway.radius = [f.sgr3 f.sgr2 f.sgr1 f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 0.013];  % total radius vector
            TotalAirway.length = [f.sgl3 f.sgl2 f.sgl1 f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7];   % total length vector
            SG.totalLength = f.sgl1 + f.sgl2 + f.sgl3; % subglottal tract length
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7;  % vocal tract length
            len(1) = f.sgl3;
            for i =2:11
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
        case 4
            TotalAirway.radius = [f.sgr4 f.sgr3 f.sgr2 f.sgr1 f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 0.013];  % total radius vector
            TotalAirway.length = [f.sgl4 f.sgl3 f.sgl2 f.sgl1 f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7];   % total length vector
            SG.totalLength = f.sgl1 + f.sgl2 + f.sgl3 + f.sgl4; % subglottal tract length
            VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7;  % vocal tract length
            len(1) = f.sgl4;
            for i =2:12
                len(i)=len(i-1)+TotalAirway.length(i);
            end
            
%         case 10
%             TotalAirway.radius = [f.sgr10 f.sgr9 f.sgr8 f.sgr7 f.sgr6 f.sgr5 f.sgr4 f.sgr3 f.sgr2 ...
%                 f.sgr1 f.gr f.r1 f.r2 f.r3 f.r4 f.r5 f.r6 f.r7 f.r8 f.r9 0.013];  % total radius vector
%             TotalAirway.length = [f.sgl10 f.sgl9 f.sgl8 f.sgl7 f.sgl6 f.sgl5 f.sgl4 f.sgl3 ...
%                 f.sgl2 f.sgl1 f.gl f.l1 f.l2 f.l3 f.l4 f.l5 f.l6 f.l7 f.l8 f.l9 f.l10];   % total length vector
%             SG.totalLength = f.sgl1 + f.sgl2 + f.sgl3 + f.sgl4 + f.sgl5 + f.sgl6 + f.sgl7 ...
%                 + f.sgl8 + f.sgl9 + f.sgl10; % subglottal tract length
%             VT.totalLength = f.l1+ f.l2+ f.l3+ f.l4+ f.l5+ f.l6+ f.l7+ f.l8+ f.l9+f.l10;  % vocal tract length
%             len(1) = f.sgl10;
%             for i =2:21
%                 len(i)=len(i-1)+TotalAirway.length(i);
%             end
    end
end

Glottis.radius = f.gr;   % glottal radius
Glottis.length = f.gl;
Glottis.alpha_multiplier =  f.galpha;

Flex.C_tissue = 4.9*10^-8; %2.9*10^-8;
Flex.L_tissue = f.w1;
Flex.R_tissue = f.w2;







%% Plot best fit
axes(handles.axes2); hold off;
% plot(handles.axes2, x,  abs(log10(Z_measured))); hold all;
plot(f, x, abs(log10(Z_measured)));
% plot(f, x, opg_imp3)
set(handles.axes2, 'YScale', 'log');
xlabel(handles.axes2, 'Frequency (Hz)');
ylabel(handles.axes2, 'Z at the lips (10^)');
axis tight;
linkaxes([handles.axes2, handles.axes3], 'x')

% axes(handles.axes3)
% plot(handles.axes3, x,  angle((Z_measured)));
% % plot(f, x, opg_imp3)
% set(handles.axes3 , 'YScale', 'linear');
% xlabel(handles.axes3, 'Frequency (Hz)');
% ylabel(handles.axes3, 'Angle (rad)');
% axis tight;


%% Plot airway profile
% create a new radius vector with twice as many values
newRad = ones(max(size(TotalAirway.radius))*2,1);
newRad(1:2:end,:) = TotalAirway.radius*10^3';
newRad(2:2:end,:) = TotalAirway.radius*10^3';

newLen = ones(max(size(len))*2,1);
newLen(1:2:end,:) = len';
newLen(2:2:end,:) = len';

% create unit cylinder with correct radius
[xVal, yVal, zVal] = cylinder(newRad, 100);

% adjust cylinder lengths so that they do not slope
% zero y-axis above of glottis
zVal(1,:) = -SG.totalLength - newLen(1) - Glottis.length;
zVal(2,:) = newLen(2) - SG.totalLength - Glottis.length;
for jj = 3:2:length(zVal(:,1))
    zVal(jj,:) = newLen(jj-1,:) - SG.totalLength - Glottis.length;
end
for kk = 4:2:length(zVal(:,1))
    zVal(kk,:) = newLen(kk) - SG.totalLength - Glottis.length;
end

% plot 3D airway profile using radius and length
axes(handles.axes1)
plot3(handles.axes1, xVal, yVal, zVal, 'b');
% h=mesh(xVal,yVal,zVal,'facecolor','texturemap');

view(handles.axes1, [0 0]); % fix 2D viewpoint
xlabel(handles.axes1, 'radius (mm)');
xlim(handles.axes1, [-60 60]);
ylabel(handles.axes1, 'radius (mm)');
zlabel(handles.axes1, 'distance from glottis (m)');
axis tight;


%% Save data

VT.length = VT.length';
VT.radius = VT.radius;
VT.alpha_multiplier = VT.alpha_multiplier';

SG.length = SG.length';
SG.radius = SG.radius;
SG.alpha_multiplier = SG.alpha_multiplier';

Glottis.length = Glottis.length';
Glottis.radius = Glottis.radius;
Glottis.alpha_multiplier = Glottis.alpha_multiplier';



Airway.VT = VT;
Airway.SG = SG;
Airway.Glottis = Glottis;


PARAMS.Temp = 20.2;
PARAMS.humidity = 0.4;
% used for the visco-thermal loss calculations
PARAMS.DeltaT = PARAMS.Temp-26.85;
PARAMS.resolution = 5.3833; % frequency resolution in Hz
PARAMS.freq = x';
PARAMS.w = 2.*pi.*PARAMS.freq; % angular frequency
PARAMS.c = 340; % standard speed of sound
PARAMS.rho = 1.14; % air density

save('FittedData', 'f', 'gof', 'fit_output', 'f_isnormal_chi2',...
    'f_P_isnormal', 'TotalAirway',...
    'Airway', 'PARAMS', 'Flex');



end
