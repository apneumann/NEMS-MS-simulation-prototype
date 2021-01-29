
%% Cantilever modeshapes (analytical)
global Mtot Mtot_ng MAX_MODES

rhoSi = 2.33e6;       % Density of Silicon in g/cm^3
Da    = 1.66054e-24;  % Mass of 1 Da in g

% I2 Nat Nano micro device.
H        = 2e-6;  
% L        = 397e-6;  % micro device, doesn't match comsol exactly
L        = 400e-6;  % comsol
W        = 30e-6; 
Mtot     = H*W*L*rhoSi;  % mass of beam in grams
Mtot_ng  = Mtot*1E9;    % mass of beam in nanograms

% need unperturbed frequency shifts now to do modeshape correction



global Ph_alpha Ph_k Ph_c Ph_max
% Modeshape constants pre-computed using Mathematica
Ph_alpha  = [.25; .25; .25; .25; .25; .25];% .25; .25; .25];
Ph_k      = [ 1.87510406871196116645;  4.69409113297417457644;  7.85475743823761256486; ...
             10.99554073487546699067; 14.13716839104647058092; 17.27875953208823633354];
             %20.42035225104125099442; 23.5619449018064435; 26.7035375555182988];
Ph_c      = [1.36222055748513128791; 0.98186753917472991616; 1.00077610535497734027; ...
             0.99996644787406952847; 1.00000144989345210284; 0.99999993734437548486];
             %1.00000000270759494807; 0.99999999988299421; 1.00000000000505628];
Ph_max    = [-2.7244411149702625758; 1.96373507834945983232; -2.0015522107099546805; ...
             1.99993289574813905693; -2.0000028997869042057; 1.99999987468875096973];
             %-2.0000000054151898961; 1.99999999976598843; -2.0000000000101126];

global Ph Ph2 G_fun dPh_2

% Modeshapes: symbolic and anonymous functions
clear x
x = sym('x');
Ph = sym('Ph',[1 MAX_MODES]);
Ph2 = sym('Ph2',[1 MAX_MODES]);
d1_Ph2 = sym('d1_Ph2',[1 MAX_MODES]);
d2_Ph2 = sym('d2_Ph2',[1 MAX_MODES]);
d3_Ph2 = sym('d3_Ph2',[1 MAX_MODES]);
dPh_2 = sym('dPh_2',[1 MAX_MODES]);
d_dPh_2 = sym('d_dPh_2',[1 MAX_MODES]);
Ph4 = sym('Ph4',[1 MAX_MODES]);
d_Ph4 = sym('d_Ph4',[1 MAX_MODES]);
G_fun = sym('G_fun',[1 MAX_MODES]);

for i=1:MAX_MODES
    Ph(i) = @(x) ( ( ( sinh( Ph_k(i)*x )-sin( Ph_k(i)*x ) )  ...
        - Ph_c(i)*( cosh( Ph_k(i)*x )-cos( Ph_k(i)*x ) ) ) / Ph_max(i) );
    Ph2(i) = Ph(i)^2;
    d1_Ph2(i) = diff(Ph2(i),x);
    d2_Ph2(i) = diff(Ph2(i),x,2);
    d3_Ph2(i) = diff(Ph2(i),x,3);
    dPh_2(i) = diff(Ph(i),x)^2;
    d_dPh_2(i) = diff(dPh_2(i),x);
    Ph4(i) = Ph(i)^4;
    d_Ph4(i) = diff(Ph4(i),x);
    
    
    G_fun(i) = Ph4(i);
    
    for j=1:MAX_MODES
        if j==i, continue; end
        coeff_ij = f0(i)^2/(f0(j)^2-f0(i)^2)*Ph_alpha(j)/Ph_alpha(i);
        G_fun(i) = G_fun(i) - coeff_ij * Ph2(i) * Ph2(j);
    end
end

%%
% plot(u,PH2(1,:),'k-');  hold on
% plot(u,PH2(2,:),'r-');
% plot(u,PH2(3,:),'b-');
% plot(u,PH2(4,:),'g-');
% plot(u,PH2(5,:));
% plot(u,PH2(6,:));
% % legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');
% 
% % clear Ph Ph2 d1_Ph2 d2_Ph2 d3_Ph2 Ph4 d_Ph4 dPh_2 d_dPh_2 
% clear posn_start posn_end mass_points posn_points 
% clear rhoSi Da x i j H L W
