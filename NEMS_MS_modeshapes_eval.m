
%% Numeric evaluation of modeshapes at predetermined grid points
global MAX_MODES useGpu

% Grid points. m = mass axis (vertical). u = posn axis (horizontal)
% mda = mass axis in units of megadaltons
global m u

% mass_points = 2401;     % number of grid points for mass axis
mass_points = 401;     % number of grid points for mass axis
% posn_points = 1200*.5+1;    % number of grid points for position
posn_points = 400+1;    % number of grid points for position

m = linspace(mass_center*(1-mass_range),mass_center*(1+mass_range),mass_points)';
u = linspace(posn_center*(1-posn_range),posn_center*(1+posn_range),posn_points);
% m = linspace(mass_left,mass_right, mass_points)';
% u = linspace(posn_left,posn_right,posn_points);

%%
global PH2 D1PH2 D2PH2 D3PH2

% rows are modes 1, 2, 3, etc.
% columns are modeshapes evaluated at positions given in u vector.

if useGpu
    u_gpu = gpuArray(u);
    PH2 = zeros(MAX_MODES, length(u)); D1PH2 = zeros(MAX_MODES,length(u));
    PH2 = gpuArray(PH2); D1PH2 = gpuArray(D1PH2);

    for i=1:MAX_MODES
        PH2(i,:) = arrayfun(matlabFunction(Ph2(i)),u_gpu);
        D1PH2(i,:) = arrayfun(matlabFunction(d1_Ph2(i)),u_gpu);
        D2PH2(i,:) = arrayfun(matlabFunction(d2_Ph2(i)),u_gpu);
        D3PH2(i,:) = arrayfun(matlabFunction(d3_Ph2(i)),u_gpu);
    end
else
    x = u;
    for i=1:MAX_MODES
        PH2(i,:)   = double(subs(Ph2(i)));
        D1PH2(i,:) = double(subs(d1_Ph2(i)));
        D2PH2(i,:) = double(subs(d2_Ph2(i)));
        D3PH2(i,:) = double(subs(d3_Ph2(i)));
    end
end

%%
global DPH2 DDPH2 PH4 D1PH4 G DG m_gpu u_gpu

if useGpu
    u_gpu = gpuArray(u);
    m_gpu = gpuArray(m);
    PH2 = zeros(MAX_MODES,length(u)); D1PH2 = zeros(MAX_MODES,length(u));
    D2PH2 = zeros(MAX_MODES,length(u)); D3PH2 = zeros(MAX_MODES,length(u));
    DPH2 = zeros(MAX_MODES,length(u)); DDPH2 = zeros(MAX_MODES,length(u));
    PH4 = zeros(MAX_MODES,length(u)); D1PH4 = zeros(MAX_MODES,length(u));
    
    PH2 = gpuArray(PH2); D1PH2 = gpuArray(D1PH2); D2PH2 = gpuArray(D2PH2);
    D3PH2 = gpuArray(D3PH2); DPH2 = gpuArray(DPH2); DDPH2 = gpuArray(DDPH2);
    PH4 = gpuArray(PH4); D1PH4 = gpuArray(D1PH4);
    
    for i=1:MAX_MODES
        PH2(i,:)   = arrayfun(matlabFunction(Ph2(i)),u_gpu);
        D1PH2(i,:) = arrayfun(matlabFunction(d1_Ph2(i)),u_gpu);
        DPH2(i,:) = arrayfun(matlabFunction(dPh_2(i)),u_gpu);
        DDPH2(i,:) = arrayfun(matlabFunction(d_dPh_2(i)),u_gpu);
        PH4(i,:) = arrayfun(matlabFunction(Ph4(i)),u_gpu);
        D1PH4(i,:) = arrayfun(matlabFunction(d_Ph4(i)),u_gpu);      
    end
else
    x = u;
    for i=1:MAX_MODES
        PH2(i,:)   = double(subs(Ph2(i)));
        D1PH2(i,:) = double(subs(d1_Ph2(i)));
        DPH2(i,:) = double(subs(dPh_2(i)));
        DDPH2(i,:) = double(subs(d_dPh_2(i)));
        PH4(i,:) = double(subs(Ph4(i)));
        D1PH4(i,:) = double(subs(d_Ph4(i)));
    end
end
  
for i=1:MAX_MODES
    G(i,:) = PH4(i,:);
    DG(i,:) = D1PH4(i,:);
    G = gpuArray(G);
    DG = gpuArray(DG);
    
    for j=1:MAX_MODES
        if j==i, continue; end
        coeff_ij = f0(i)^2/(f0(j)^2-f0(i)^2)*Ph_alpha(j)/Ph_alpha(i);
        G(i,:) = G(i,:) - coeff_ij.*PH2(i,:).*PH2(j,:);
        DG(i,:) = DG(i,:) - coeff_ij.*(D1PH2(i,:).*PH2(j,:) + PH2(i,:).*D1PH2(j,:));
    end
end

clear coeff_ij

%%
% plot(u,PH2(1,:),'k-');  hold on
% plot(u,PH2(2,:),'r-');
% plot(u,PH2(3,:),'b-');
% plot(u,PH2(4,:),'g-');
% plot(u,PH2(5,:));
% plot(u,PH2(6,:));
% legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6');

% clear Ph Ph2 d1_Ph2 d2_Ph2 d3_Ph2 Ph4 d_Ph4 dPh_2 d_dPh_2 
clear posn_start posn_end mass_points posn_points 
clear rhoSi Da x i j H L W


