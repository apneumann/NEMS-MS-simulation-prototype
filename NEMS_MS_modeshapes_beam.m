
%% Doubly-clamped beam modeshapes (analytical)
global Mtot MAX_MODES

rhoSi = 2.33e6;       % Density of Silicon in g/cm^3
Da    = 1.66054e-24;  % Mass of 1 Da in g
kDa   = Da*1e3;       % Mass of 1 kDa in g
MDa   = Da*1e6;       % Mass of 1 MDa in g

% Jarvis best device.
H        = 50e-9; 
L        = 4e-6;
W        = 300e-9; 
Mtot     = H*W*L*rhoSi;  % mass of beam in grams
Mtot_fg  = Mtot*1E15;    % mass of beam in femtograms
Mtot_kDa = Mtot/kDa;
Mtot_MDa = Mtot/MDa;


global Ph_alpha Ph_k Ph_c Ph_max
% Modeshape constants pre-computed using Mathematica
Ph_alpha  = [.396477920160514;   .43902788581490554; .4371036440440953
             .43718662837199773; .42453216245187664; .42647860153584083];
Ph_k      = [4.730040744862708;  7.853204624095521;  10.995607838002778
             14.13716549127998;  17.278759658031507; 20.42035225872191];
Ph_c      = [1.0178094106701927; 0.9992232918372344; 1.000033551000223
             0.999998550104456;  1.0000000626556207; 0.9999999972924051];
Ph_max    = [1.6164302110500461; 1.5080526207840705; 1.5125939437711535
             1.5123974455422926; 1.5124059365329394; 1.5124055696030894];
            
global Ph2 d1_Ph2 d2_Ph2 d3_Ph2

% Modeshapes: symbolic and anonymous functions
clear x;
x = sym('x');
Ph2 = sym('Ph2',[1 MAX_MODES]);
d1_Ph2 = sym('d1_Ph2',[1 MAX_MODES]);
d2_Ph2 = sym('d2_Ph2',[1 MAX_MODES]);
d3_Ph2 = sym('d3_Ph2',[1 MAX_MODES]);

for i=1:MAX_MODES
    Ph2(i) = @(x) ( ( ( sinh( Ph_k(i)*x )-sin( Ph_k(i)*x ) )  ...
        + Ph_c(i)*( cos( Ph_k(i)*x )-cosh( Ph_k(i)*x ) ) ) / Ph_max(i) ).^2;
    d1_Ph2(i) = diff(Ph2(i),x);
    d2_Ph2(i) = diff(Ph2(i),x,2);
    d3_Ph2(i) = diff(Ph2(i),x,3);
end

