% Pull data from (deterministic) equation, and recover particle's mass, 
% position, and 2nd linear moment (via inertial imaging theory).

global MAX_MODES useGpu Mtot_kDa GroEL_rel SIGMA Ph_alpha Ph2 d2_Ph2
MAX_MODES = 6;
useGpu = 1;

NEMS_MS_modeshapes_beam;

% global m u

adev     = 1.5E-7;
rho      = 0;   % noise correlation coeff

for i = 1:MAX_MODES
    for j = 1:MAX_MODES
        if i==j, SIGMA(i,j) = adev^2;
        else, SIGMA(i,j) = adev^2*rho;
        end
    end
end


% Generate single event frequent shift without noise
GroEL_rel= 800/Mtot_kDa;
num_events = 1;%20;
mass_input = GroEL_rel;  % particle mass relative to beam
% posn_input = linspace(.2,.48,num_events);%.42;
% posn_input = linspace(0,.5,num_events);%.42;
posn_input = .42;
ext_input = 0.01;
digits(15);

df_p = zeros(num_events,MAX_MODES);
MU = zeros(1,MAX_MODES);

fprintf('integral:         ');
for i=1:num_events
    fprintf('\b\b\b\b\b\b\b%5.02f%%\n',i/num_events*100);
    xi = posn_input(i);
    a = ext_input;
    for j=1:MAX_MODES
        fun = matlabFunction(symfun(Ph2(j),x));
        df_p(i,j) = -integral(fun,xi-a/2,xi+a/2)*mass_input/(2*Ph_alpha(j)*a);
    end
end

df_noise  = mvnrnd(MU,SIGMA,num_events);
df_n = df_p + df_noise;



% only keep events outside 99.99% confidence
X3 = df_n(:,1:3);
X6 = df_n;
ch3 = chi2inv(.9999,3);
ch6 = chi2inv(.9999,6);

[V3,D3] = eig(SIGMA(1:3,1:3)*ch3);
[V6,D6] = eig(SIGMA*ch6);

pdf_at_confidence3 = mvnpdf(V3(:,1)'*sqrt(D3(1,1)),MU(1:3),SIGMA(1:3,1:3));
pdf_at_confidence6 = mvnpdf(V6(:,1)'*sqrt(D6(1,1)),MU,SIGMA);

pdf_eval3 = mvnpdf(X3,MU(1:3),SIGMA(1:3,1:3));
pdf_eval6 = mvnpdf(X6,MU,SIGMA);

ind3 = pdf_eval3 > pdf_at_confidence3;
ind6 = pdf_eval6 > pdf_at_confidence6;

X3(ind3,:) = [];
X6(ind6,:) = [];

df_n3 = X3;
df_n6 = X6;



NEMS_MS_modeshapes_eval;

%%
syms mm x vv
sol = vpasolve([df_p(1)==-mm/(2*Ph_alpha(1))*(Ph2(1)+.5*d2_Ph2(1)*vv), ...
                df_p(2)==-mm/(2*Ph_alpha(2))*(Ph2(2)+.5*d2_Ph2(2)*vv), ...
                df_p(3)==-mm/(2*Ph_alpha(3))*(Ph2(3)+.5*d2_Ph2(3)*vv)], ...
                [mm, x,vv],[1e-5;.4;0]);
            
[sol.mm*Mtot_kDa sol.x sol.vv]

% if( ~isempty(sol.x) && sol.x > 0.3 && sol.x < 0.5 )
%     m2_pos = [m2_pos; double(sol.x)];
%     m2_avg = [m2_avg; double(sol.mm) * Mtot_kDa];
% end

%% 2-mode mass & position spectra combined using 4 freq shifts.

m2_pos = []; m2_avg = []; m2_res = [];
m2_good_pos = []; m2_good_avg = []; m2_good_res = [];
m6_pos = []; m6_avg = []; m6_res = []; 
m6_good_pos = []; m6_good_avg = []; m6_good_res = [];
pdf_m2_all = 0; pdf_m2_goodpos = 0; pdf_m6_all = 0; pdf_m6_goodpos = 0;

pdf_a2 = [];
pdf_a6 = [];
ext_grid = linspace(0,ext_input*3,100);

jpdf_mu = ones(length(m),length(u));
if useGpu, jpdf_mu = gpuArray(jpdf_mu); end

syms x mm

% fprintf('jpdf_beam2:         ');
for i=1%1:size(df_n2,1)
%     fprintf('\b\b\b\b\b\b\b%5.02f%%\n',i/size(df_n3,1)*100);
    
    df3 = df_n3;%(i,:);
    df6 = df_n6;
    
    for ai = ext_grid
        jpdf_mu = jpdf_sqrtvar_combine(3,df3,ai,SIGMA);
        pdf_a2 = [pdf_a2; gather(sum(sum((jpdf_mu))))];
        jpdf_mu = jpdf_sqrtvar_combine(6,df6,ai,SIGMA);
        pdf_a6 = [pdf_a6; gather(sum(sum((jpdf_mu))))];
    end
    
    pdf_a2 = pdf_a2/sum(pdf_a2);
    pdf_a6 = pdf_a6/sum(pdf_a6);
end


% [std(m2_good_avg) std(m6_good_avg) mean(m2_good_res) mean(m6_good_res)]



figure;
plot(ext_grid,pdf_a2,'k');
hold on
plot(ext_grid,pdf_a6,'b');
xlabel('Normalized second moment');


% figure;
% hold off
% plot(m2_pos,m2_avg,'k.');
% hold on
% plot(m6_pos,m6_avg,'b.');
% ylim([700 900]);
% title('Mass (kDa) vs calculated position');
% xlabel('Calculated position');
% legend('2 modes', '6 modes');
% 
% figure;
% plot(m2_pos,m2_res,'k.');
% hold on
% plot(m6_pos,m6_res,'b.');
% ylim([0 50]);
% title('Mass resolution (kDa) vs calculated position');
% xlabel('Calculated position');
% legend('2 modes', '6 modes');



% figure;
% plot(m*Mtot_kDa,pdf_m2_goodpos,'k');
% hold on
% plot(m*Mtot_kDa,pdf_m6_all,'b');
% title('Ensemble mass spectrum (1000 events)');
% xlabel('Mass (kDa)');
% legend('2 modes (x>0.3)', '6 modes');

% ensemble
% [m2a,w2a] = find_mean_and_width(m*Mtot_kDa,pdf_m2_all);
% [m2g,w2g] = find_mean_and_width(m*Mtot_kDa,pdf_m2_goodpos);
% [m6a,w6a] = find_mean_and_width(m*Mtot_kDa,pdf_m6_all);
% [w2g w6a]


% figure;
% plot(u,pdf_u2,'k');
% hold on
% plot(u,pdf_u6,'b');

%%
% pdf_m2_all = pdf_m2_all / (sum(pdf_m2_all)*(m(2)-m(1))*Mtot_kDa);
% pdf_m2_goodpos = pdf_m2_goodpos / (sum(pdf_m2_goodpos)*(m(2)-m(1))*Mtot_kDa);
% pdf_m6_all = pdf_m6_all / (sum(pdf_m6_all)*(m(2)-m(1))*Mtot_kDa);
% 
% % [mean(m2_res(500:1000)) mean(m6_res) std(m2_avg(500:1000)) std(m6_avg)]
% [mean(m2_avg) std(m2_avg) mean(m2_pdf_avg) std(m2_pdf_avg)]
% 
% kernel_pdf = zeros(length(m),1);
% n = length(m2_avg);
% h = 1.06 * 1.5* std(m2_avg/Mtot_kDa) *(n)^(-0.2);
% m_gpu= gpuArray(m);
% kernel_gpu = gpuArray(kernel_pdf);
% 
% clear mm ki
% mm = sym('mm');
% ki = sym('ki');
% 
% fprintf('kernel_pdf:         ');
% for i=1:n
%     fprintf('\b\b\b\b\b\b\b%5.02f%%\n',i/n*100);
%     mi = m2_avg(i)/Mtot_kDa;
%     ki = @(mm) exp( -(mm-mi)^2/(2*h^2) ) / (2*pi*sqrt(h));
% %     ki = @(mm) .75*max( 1-(mm-mi)^2/h^2, 0 );
% %     ki = @(mm) piecewise(.75*( 1-(mm-mi)^2/h^2 ) > 0, .75*( 1-(mm-mi)^2/h^2 ),0);
%     mm = m;
%     kernel_pdf = kernel_pdf + double(subs(ki));
% %     kernel_gpu = kernel_gpu + arrayfun(ki,m_gpu);
% end
% 
% % kernel_pdf = gather(kernel_gpu) / (n*h);
% kernel_pdf = kernel_pdf / (n*h);
% 
% [mh,wh]=find_mean_and_width(m*Mtot_kDa,kernel_pdf);
% [mg,wg] = find_mean_and_width(m*Mtot_kDa,pdf_m2_goodpos);
% [mh wh mg wg]


%%
% % h=histfit(m2_avg,50,'kernel');
% % hdata = h(2).YData;
% h = hist(m2_avg);
% hmax=max(h);
% figure;
% hist(m2_avg);
% % histfit(m2_avg,50,'kernel');
% hold on;
% plot(m*Mtot_kDa,kernel_pdf/max(kernel_pdf)*hmax,'r');
% plot(m*Mtot_kDa,pdf_m2_goodpos/max(pdf_m2_goodpos)*hmax,'k');
% xlim([725 875]);
% title('Ensemble mass spectrum (1000 events)');
% xlabel('Mass (kDa)');

%%

% options = statset('Display','final');
% gm = fitgmdist(m2_avg,1);
% gmpdf = @(mm) arrayfun(@(m0) pdf(gm,m0),mm);
% h=hist(m2_avg);
% hmax=max(h);
% figure;
% hist(m2_avg);
% hold on
% gm_pdf = gmpdf(m*Mtot_kDa);
% plot(m*Mtot_kDa,gm_pdf/max(gm_pdf)*hmax,'r');
% plot(m*Mtot_kDa,pdf_m2_goodpos/max(pdf_m2_goodpos)*hmax,'k');
% xlim([725 875]);
% title('Ensemble mass spectrum (1000 events)');
% xlabel('Mass (kDa)');
% [mh,wh]=find_mean_and_width(m*Mtot_kDa,gm_pdf);
% [mg,wg] = find_mean_and_width(m*Mtot_kDa,pdf_m2_goodpos);
% [mh wh mg wg]

%%



% 
% figure;
% plot(m2_pdf_pos,(m2_pdf_avg/800-1)*100,'k.');
% hold on
% plot(m2_pos,(m2_avg/800-1)*100,'b.');
% title('Mass error (pct) vs calculated position');
% xlabel('Calculated position');
% legend('JPDF method', 'Simple solve method');





