% @author: Nichlas Vous Christensen
% @email: nvc@clin.au.dk
% @phone: +45 23464522
% @organization: Aarhus University, The MR Research Centre
% June 2025

%% Simulation
TR = 1; % s
time_max = 180; % s
t = linspace(0,time_max-TR,time_max/TR);
flip_P = 10;
flip_L = 10;
kpl = 0.02;
klp = 0;
r1p = 1/30;
r1l = 1; % mutiply of rp
k = 10; % multiply of kpl
t0_model = 0;
mu = 5;
sigma = 3;
noise_scale = 0.03; % scaled from maximum signal
noise_realizations = 10; % number of repeated noise simulations

r2p = -log(cosd(flip_P))/TR; % Relaxation due to pulsing, see Hill et al. 2013 (Model Free Approach)
r2l = -log(cosd(flip_L))/TR;

rp = r1p + r2p; % Effective relaxation rate on pyruvate
rl = r1l*r1p + r2l; % Effective relaxation rate on lactate

t0 = 10;
a = 3;
b = 1.5;
A = 1;
params_AIF = [t0, a, b, A];
AIF = AIF_model(params_AIF, t);
AIF = AIF/max(AIF); % Normalize

params_sim = [kpl, klp, r1p, r1l, k, t0_model, mu, sigma];

y_sim = HEMEX_model_plot(params_sim, t, AIF, r2p, r2l);

rng(42) % Set same random seed everytime for reproducible results

%% lower bound | start guess | upper bound
bounds = [0 0.01 1;... % kpl (pyrurvate-to-lactate conversion rate)
          0 0 0;... % klp (lactate-to-pyruvate conversion rate)
          0.01 1/30 0.05;... % rp (1/T1 pyruvate relaxation)
          0.5 1 2;... % rl (1/T1 lactate relaxation) [scaled from rp]
          0.05 1 20;... % k (pyruvate permability) [scaled from kpl]
          0 0 10;... % t0_delay (time-delay between AIF and voxel data)
          0.5 5 30;... % mu (input parameter to the hemodynamic residue function)
          0.5 3 30]; % sigma (input parameter to the hemodynamic residue function)
params_lb = bounds(:,1);
params0 = bounds(:,2);
params_ub = bounds(:,3);

%% Noisy simulation
kpl_fits = [];
k_fits = [];
mu_fits = [];
kpl_ratio = [];

% First go from Mz to Mxy
Pxy = y_sim(2,:)*sind(flip_P);
Lxy = y_sim(1,:)*sind(flip_L);

for i=1:noise_realizations
    noise_level = max(Pxy+Lxy)*noise_scale;
    
    noise_L = normrnd(0,noise_level,1,length(Pxy));
    noise_P = normrnd(0,noise_level,1,length(Lxy));
    
    Lxy_noisy = Lxy+noise_L;
    Pxy_noisy = Pxy+noise_P;
    
    % Fitting
    options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16,'Display','off');
    
    % Take into account the flip-angles (conversion from Mxy to Mz)
    Pz_noisy = Pxy_noisy/sind(flip_P);
    Lz_noisy = Lxy_noisy/sind(flip_L);
    
    params_fit = lsqcurvefit(@(params, t) HEMEX_model(params, t, AIF, r2p, r2l), params0, t, [Lz_noisy; Pz_noisy],params_lb,params_ub,options);
    
    kpl_fit = params_fit(1);
    klp_fit = params_fit(2);
    rp_fit = params_fit(3);
    rl_fit = params_fit(4);
    k_fit = params_fit(5);
    mu_fit = params_fit(7);
    sigma_fit = params_fit(8);
    
    kpl_fits(i) = kpl_fit;
    k_fits(i) = k_fit*kpl_fit;
    mu_fits(i) = mu_fit;
    
    kpl_ratio(i) = sum(Lz_noisy)/sum(Pz_noisy)*rl;
end

kpl_fit_mean = mean(kpl_fits);
kpl_fit_std = std(kpl_fits);
k_fit_mean = mean(k_fits);
k_fit_std = std(k_fits);
mu_fit_mean = mean(mu_fits);
mu_fit_std = std(mu_fits);

kpl_ratio_mean = mean(kpl_ratio);
kpl_ratio_std = std(kpl_ratio);

disp('=======================================')
fprintf('kPL (true): %.3f\n', kpl);
fprintf('kPL (HEMEX): %.3f ± %.3f\n', kpl_fit_mean, kpl_fit_std);
fprintf('kPL (ratio-metric): %.3f ± %.3f\n', kpl_ratio_mean, kpl_ratio_std);
disp('=======================================')

%% Plotting 
fig = figure('Renderer', 'painters', 'Position', [400 300 700 250]);
subplot(1,2,1)
plot(t,y_sim(1,:),'b-','linewidth',2)
hold on
plot(t,y_sim(3,:),'r-','linewidth',2)
plot(t,y_sim(4,:),'k-','linewidth',2)
plot(t,y_sim(2,:),':','linewidth',2)

title(['k_P_L = ',num2str(kpl,'%.3f'),' s^-^1, k = ',num2str(k*kpl,'%.3f'),' s^-^1',newline,'MTT = ',num2str(mu,'%.3f'),' s'],'fontsize',9)
legend('Lactate','Pyruvate (tissue)','Pyruvate (vascular)', 'Pyruvate (total)')
xlabel('Time (s)')
ylabel('Magnetization (M_Z)')
axis([0 max(t) 0-max(y_sim(:))*0.1 max(y_sim(:))*1.1])

subplot(1,2,2)
plot(t,Lxy_noisy,'b.')
hold on
plot(t,Pxy_noisy,'r.')

y_fit = HEMEX_model(params_fit, t, AIF, r2p, r2l);

legend('Lactate noisy','Pyruvate noisy')
xlabel('Time (s)')
ylabel('Signal (M_X_Y)')
axis([t(1), t(end), 0-max(Pxy)*0.1, max(Pxy)*1.1])

% Converting back to Mxy from Mz when plotting
plot(t,y_fit(1,:)*sind(flip_L),'b-')
plot(t,y_fit(2,:)*sind(flip_P),'r-')
legend('Lactate noisy','Pyruvate noisy','Lactate fit','Pyruvate fit')
title(['k_P_L = ',num2str(kpl_fit_mean,'%.3f'),' ',char(177),' ',num2str(kpl_fit_std,'%.3f'),' s^-^1, k = ',num2str(k_fit_mean,'%.3f'),' ',char(177),' ',num2str(k_fit_std,'%.3f'),' s^-^1',newline,'MTT = ',num2str(mu_fit_mean,'%.3f'),' ',char(177),' ',num2str(mu_fit_std,'%.3f'),'s'],'fontsize',9)

%% Model functions
function y = HEMEX_model(params, t, AIF, r2p, r2l)
    kpl = params(1);
    klp = params(2);
    rp = params(3) + r2p; % Take relaxation due to pulsing into account
    rp_scale = params(4); % Reparameterization of rp-rl dependency
    rl = rp_scale * rp + r2l; % Take relaxation due to pulsing into account
    kpl_scale = params(5);
    k = kpl_scale * kpl; % Reparameterization of kpl-k dependency
    t0_delay = params(6);
    mu = params(7);
    sigma = params(8);

    alpha = mu^2/sigma^2;
    beta = sigma^2/mu;

    AIF_shifted = interp1(t, AIF, t - t0_delay, 'linear', 0); % Delay between AIF and given voxel

    R = @(t)gammainc(t/beta, alpha, 'upper').*(t>=0);
    Pv = (1/mu)*conv(AIF_shifted,exp(-(k+rp)*t).*R(t),'full');
    Pv = Pv(1:length(t)); % Pyruvate vascular

    k_plus = -1/2*(klp+kpl+rl+rp) + 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));
    k_minus = -1/2*(klp+kpl+rl+rp) - 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));

    L = k*kpl/(k_plus-k_minus)*conv(Pv,exp(k_plus*t)-exp(k_minus*t),'full'); % Lactate tissue
    Pt = k/(k_plus-k_minus)*conv(Pv,(k_plus+klp+rl)*exp(k_plus*t)-(k_minus+klp+rl)*exp(k_minus*t),'full'); % Pyruvate tissue

    y(1,:) = L(1:length(t));  % Total lactate signal
    y(2,:) = Pt(1:length(t))+Pv(1:length(t)); % Total pyruvate signal
end

function y = HEMEX_model_plot(params, t, AIF, r2p, r2l)
    kpl = params(1);
    klp = params(2);
    rp = params(3) + r2p; % Take relaxation due to pulsing into account
    rp_scale = params(4); % Reparameterization of rp-rl dependency
    rl = rp_scale * rp + r2l; % Take relaxation due to pulsing into account
    kpl_scale = params(5);
    k = kpl_scale * kpl; % Reparameterization of kpl-k dependency
    t0_delay = params(6);
    mu = params(7);
    sigma = params(8);

    alpha = mu^2/sigma^2;
    beta = sigma^2/mu;
    
    AIF_shifted = interp1(t, AIF, t - t0_delay, 'linear', 0); % Delay between AIF and given voxel

    R = @(t)gammainc(t/beta, alpha, 'upper').*(t>=0);
    Pv = (1/mu)*conv(AIF_shifted,exp(-(k+rp)*t).*R(t),'full');
    Pv = Pv(1:length(t)); % Pyruvate vascular

    k_plus = -1/2*(klp+kpl+rl+rp) + 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));
    k_minus = -1/2*(klp+kpl+rl+rp) - 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));

    L = k*kpl/(k_plus-k_minus)*conv(Pv,exp(k_plus*t)-exp(k_minus*t),'full'); % Lactate tissue
    Pt = k/(k_plus-k_minus)*conv(Pv,(k_plus+klp+rl)*exp(k_plus*t)-(k_minus+klp+rl)*exp(k_minus*t),'full'); % Pyruvate tissue

    y(1,:) = L(1:length(t));  % Total lactate signal
    y(2,:) = Pt(1:length(t))+Pv(1:length(t)); % Total pyruvate signal
    y(3,:) = Pt(1:length(t)); % Pyruvate tissue
    y(4,:) = Pv(1:length(t)); % Pyruvate vascular
end

function [AIF] = AIF_model(params, t)
    a = params(1);
    b = params(2);
    t0 = params(3);
    A = params(4);

    AIF = A*(t-t0 > 0).*(t-t0).^(a-1).*exp(-(t-t0)/b);
end