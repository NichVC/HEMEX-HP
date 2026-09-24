function model_fit_results = HEMEX_main_4met_dualkidney_new_constraining(...
    AIF_raw, ...
    pyr_norm_kidney1, lac_norm_kidney1, bic_norm_kidney1, ala_norm_kidney1, ...
    pyr_norm_kidney2, lac_norm_kidney2, bic_norm_kidney2, ala_norm_kidney2, ...
    TR, flip_P, flip_L, flip_B, flip_A, varargin)
    % HEMEX_main_4met_dualkidney
    %   Fit the 4-metabolite HEMEX model (Pyr, Lac, Bic, Ala) to TWO kidney ROIs
    %   simultaneously using a shared parameter set, except for separate t0_delay
    %   terms for kidney 1 and kidney 2.
    %
    %   The residual vector is constructed by stacking the selected metabolite
    %   curves from kidney 1 and kidney 2. This means the optimizer finds one
    %   parameter set that best explains both kidneys at the same time.
    %
    %   Usage:
    %     s = HEMEX_main_4met_dualkidney(AIF_raw, ...
    %         pyr1, lac1, bic1, ala1, pyr2, lac2, bic2, ala2, ...
    %         TR, flipP, flipL, flipB, flipA);
    %
    %     s = HEMEX_main_4met_dualkidney(..., 'FitMetabolites', {'P','L'});
    %     s = HEMEX_main_4met_dualkidney(..., 'AIFMethod','gamma', 'Bounds', myBounds);
    %
    %   Parameter order:
    %     1  kpl
    %     2  klp
    %     3  rp_fit
    %     4  rl_scale
    %     5  k_scale
    %     6  t0_delay_kidney1
    %     7  t0_delay_kidney2
    %     8  mu
    %     9  sigma
    %     10 scale_kidney1
    %     11 scale_kidney2
    %     12 kpb
    %     13 kbp
    %     14 kpa
    %     15 kap
    %     16 rb_scale
    %     17 ra_scale

    % -------------------- Parse options --------------------
    opts = parseOptions4met_dualkidney(varargin{:});

    % -------------------- Time vector ----------------------
    nTime = numel(AIF_raw);
    t = linspace(0, nTime*TR - TR, nTime);

    % -------------------- Pulse contributions --------------
    r2p = -log(cosd(flip_P)) / TR;
    r2l = -log(cosd(flip_L)) / TR;
    r2b = -log(cosd(flip_B)) / TR;
    r2a = -log(cosd(flip_A)) / TR;

    % -------------------- Build AIF_fit --------------------
    [AIF_fit, aif_meta] = makeAIF(t, AIF_raw, opts);

    % -------------------- Prepare data ---------------------
    % Row order matches the model output: [L; P; B; A]
    data_all_k1 = [ ...
        lac_norm_kidney1(:)';
        pyr_norm_kidney1(:)';
        bic_norm_kidney1(:)';
        ala_norm_kidney1(:)'];

    data_all_k2 = [ ...
        lac_norm_kidney2(:)';
        pyr_norm_kidney2(:)';
        bic_norm_kidney2(:)';
        ala_norm_kidney2(:)'];

    % -------------------- Select metabolites to fit --------
    [sel_idx, fit_labels] = mapFitMetabolites(opts.FitMetabolites);

    data_sel_k1 = data_all_k1(sel_idx, :);
    data_sel_k2 = data_all_k2(sel_idx, :);

    % Shared metabolite weighting across both kidneys.
    % This keeps the inter-kidney amplitude differences intact so the
    % kidney-specific scale factors remain identifiable.
    maxVals = max([data_sel_k1, data_sel_k2], [], 2);
    maxVals(maxVals <= 0) = 1;
    weights = 1 ./ maxVals;

    target = [ ...
        data_sel_k1 .* weights; ...
        data_sel_k2 .* weights];

    % -------------------- Bounds/init ----------------------
    params_lb = opts.Bounds(:,1);
    params0   = opts.Bounds(:,2);
    params_ub = opts.Bounds(:,3);

    options = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-16, ...
        'FunctionTolerance',1e-16, ...
        'Display','off', ...
        'MaxFunctionEvaluations',10000, ...
        'MaxIterations',6000);

    % -------------------- Fit ------------------------------
    model_fun = @(p,tt) localModelSel_dual(...
        p, tt, AIF_fit, r2p, r2l, r2b, r2a, sel_idx, weights);

    [pars,~,~,exitflag] = lsqcurvefit(...
        model_fun, params0, t, target, params_lb, params_ub, options);

    % -------------------- Predicted curves -----------------
    y_pred_k1 = HEMEX_model_4met_singlekidney(...
        pars, t, AIF_fit, r2p, r2l, r2b, r2a, pars(6), pars(10));

    y_pred_k2 = HEMEX_model_4met_singlekidney(...
        pars, t, AIF_fit, r2p, r2l, r2b, r2a, pars(7), pars(11));

    lac_fits_kidney1 = y_pred_k1(1,:);
    pyr_fits_kidney1 = y_pred_k1(2,:);
    bic_fits_kidney1 = y_pred_k1(3,:);
    ala_fits_kidney1 = y_pred_k1(4,:);

    lac_fits_kidney2 = y_pred_k2(1,:);
    pyr_fits_kidney2 = y_pred_k2(2,:);
    bic_fits_kidney2 = y_pred_k2(3,:);
    ala_fits_kidney2 = y_pred_k2(4,:);

    % -------------------- Normalized R^2 -------------------
    Nt = numel(t);

    % Overall R^2 over both kidneys
    y_norm = [];
    y_pred_norm = [];

    for kk = 1:2
        if kk == 1
            data_all = data_all_k1;
            y_pred   = y_pred_k1;
        else
            data_all = data_all_k2;
            y_pred   = y_pred_k2;
        end

        for ii = 1:numel(sel_idx)
            r = sel_idx(ii);
            d = data_all(r, :);
            p = y_pred(r, :);

            npair = normalize([d, p], 'range')';
            y_norm      = [y_norm;      npair(1:Nt)]; %#ok<AGROW>
            y_pred_norm = [y_pred_norm; npair(Nt+1:end)]; %#ok<AGROW>
        end
    end

    SStot = sum((y_norm(:) - mean(y_norm(:))).^2);
    SSres = sum((y_norm(:) - y_pred_norm(:)).^2);
    R2    = 1 - SSres / SStot;

    % Per-kidney R^2 (useful for diagnosing mismatch)
    R2_kidney1 = calcNormalizedR2(data_all_k1, y_pred_k1, sel_idx);
    R2_kidney2 = calcNormalizedR2(data_all_k2, y_pred_k2, sel_idx);

    % -------------------- Reparameterization ----------------
    kpl   = pars(1);
    klp   = pars(2);
    rp0   = pars(3);
    rl_sc = pars(4);
    E  = pars(5);
    t0d_k1 = pars(6);
    t0d_k2 = pars(7);
    mu    = pars(8);
    sig_sc   = pars(9);

    scale_kidney1 = pars(10);
    scale_kidney2 = pars(11);

    kpb   = pars(12);
    kbp   = pars(13);
    kpa   = pars(14);
    kap   = pars(15);
    rb_sc = pars(16);
    ra_sc = pars(17);

    rp_eff = rp0 + r2p;
    rl_eff = rl_sc * rp_eff + r2l;
    rb_eff = rb_sc * rp_eff + r2b;
    ra_eff = ra_sc * rp_eff + r2a;

    sig_eff = sig_sc*mu;
    
    alpha = mu^2 / sig_eff^2;
    beta  = sig_eff^2 / mu;

    k_eff = ((1-E)^(-1/alpha)-1)/beta;     

    % -------------------- Package results -------------------
    model_fit_results.par_fits       = pars;
    model_fit_results.exitflag_fits  = exitflag;

    % fitted metabolite selection
    model_fit_results.fit_metabolites = fit_labels;
    model_fit_results.fit_rows        = sel_idx;

    % exchange rates
    model_fit_results.kpl_fits = kpl;
    model_fit_results.klp_fits = klp;
    model_fit_results.kpb_fits = kpb;
    model_fit_results.kbp_fits = kbp;
    model_fit_results.kpa_fits = kpa;
    model_fit_results.kap_fits = kap;

    % intrinsic relax
    model_fit_results.r1p_fits      = rp0;
    model_fit_results.r1l_fits      = rl_sc * rp0;
    model_fit_results.r1b_fits      = rb_sc * rp0;
    model_fit_results.r1a_fits      = ra_sc * rp0;
    model_fit_results.rl_scale_fits = rl_sc;
    model_fit_results.rb_scale_fits = rb_sc;
    model_fit_results.ra_scale_fits = ra_sc;

    % effective relax incl. pulsing
    model_fit_results.r1p_eff = rp_eff;
    model_fit_results.r1l_eff = rl_eff;
    model_fit_results.r1b_eff = rb_eff;
    model_fit_results.r1a_eff = ra_eff;

    % extraction / hemodynamics
    model_fit_results.k_fits  = k_eff;
    model_fit_results.E_fits  = E;
    model_fit_results.t0_delay_kidney1_fits = t0d_k1;
    model_fit_results.t0_delay_kidney2_fits = t0d_k2;
    model_fit_results.mu_fits       = mu;
    model_fit_results.sigma_sc_fits    = sig_sc;
    model_fit_results.sigma_fits    = sig_eff;

    model_fit_results.scale_kidney1_fits = scale_kidney1;
    model_fit_results.scale_kidney2_fits = scale_kidney2;

    % fit quality
    model_fit_results.R2 = R2;
    model_fit_results.R2_kidney1 = R2_kidney1;
    model_fit_results.R2_kidney2 = R2_kidney2;

    % curves & bookkeeping (raw)
    model_fit_results.pyr_norm_kidney1 = pyr_norm_kidney1;
    model_fit_results.lac_norm_kidney1 = lac_norm_kidney1;
    model_fit_results.bic_norm_kidney1 = bic_norm_kidney1;
    model_fit_results.ala_norm_kidney1 = ala_norm_kidney1;

    model_fit_results.pyr_norm_kidney2 = pyr_norm_kidney2;
    model_fit_results.lac_norm_kidney2 = lac_norm_kidney2;
    model_fit_results.bic_norm_kidney2 = bic_norm_kidney2;
    model_fit_results.ala_norm_kidney2 = ala_norm_kidney2;

    % curves & bookkeeping (fits)
    model_fit_results.pyr_fits_kidney1 = pyr_fits_kidney1;
    model_fit_results.lac_fits_kidney1 = lac_fits_kidney1;
    model_fit_results.bic_fits_kidney1 = bic_fits_kidney1;
    model_fit_results.ala_fits_kidney1 = ala_fits_kidney1;

    model_fit_results.pyr_fits_kidney2 = pyr_fits_kidney2;
    model_fit_results.lac_fits_kidney2 = lac_fits_kidney2;
    model_fit_results.bic_fits_kidney2 = bic_fits_kidney2;
    model_fit_results.ala_fits_kidney2 = ala_fits_kidney2;

    model_fit_results.t        = t;
    model_fit_results.AIF_raw  = AIF_raw;
    model_fit_results.AIF_fit  = AIF_fit;
    model_fit_results.AIF_info = aif_meta;

    % keep the weighted data used in the fit for reproducibility/debugging
    model_fit_results.weights = weights;
    model_fit_results.target_weighted = target;
end


% ======================================================
% ================== Local fit wrapper =================
% ======================================================

function ysel = localModelSel_dual(p, tt, AIF_fit, r2p, r2l, r2b, r2a, sel_idx, weights)
    y1 = HEMEX_model_4met_singlekidney(p, tt, AIF_fit, r2p, r2l, r2b, r2a, p(6), p(10));
    y2 = HEMEX_model_4met_singlekidney(p, tt, AIF_fit, r2p, r2l, r2b, r2a, p(7), p(11));

    ysel1 = y1(sel_idx, :) .* weights;
    ysel2 = y2(sel_idx, :) .* weights;

    ysel = [ysel1; ysel2];
end


% ======================================================
% =================== Option parsing ===================
% ======================================================

function opts = parseOptions4met_dualkidney(varargin)

    defaultBounds = [ ...
        0.005  0.05    0.2;   % 1) kpl
        0.005  0.1     0.5;   % 2) klp
        0.01   1/30    0.05;  % 3) rp_fit
        0.8    1       1.2;   % 4) rl_scale
        0.2    0.3     0.7;   % 5) E (scaled of k)
       -5      0       20;    % 6) t0_delay_kidney1
       -5      0.25    20;    % 7) t0_delay_kidney2
        1      5       30;    % 8) mu
        0.1    0.5     1.5;   % 9) sigma_scale of mu
        0.2    1       5;     % 10) scale_kidney1
        0.2    1       5;     % 11) scale_kidney2
        0      0.01    0.3;   % 12) kpb
        0      0       0;     % 13) kbp
        0      0.01    0.3;   % 14) kpa
        0      0       0;     % 15) kap
        0.8    1       1.2;   % 16) rb_scale
        0.8    1       1.2];  % 17) ra_scale

    defaultAIFBounds = [ ...
       -5   2   20;   % t0
        1   3    5;   % a
        1   1.5  3;   % b
        0.2 1    1.5];% A

    defaultFitMets = {'P','L','B','A'};

    p = inputParser;
    p.addParameter('Bounds',         defaultBounds,     @(x) isnumeric(x) && isequal(size(x,2),3));
    p.addParameter('AIFMethod',      'movmedian',       @(s) ischar(s) || isstring(s));
    p.addParameter('AIFWindow',      3,                 @(n) isnumeric(n) && isscalar(n) && n>0);
    p.addParameter('AIFBounds',      defaultAIFBounds,  @(x) isnumeric(x) && isequal(size(x),[4 3]));
    p.addParameter('FitMetabolites', defaultFitMets,    @(x) iscell(x) || isstring(x) || ischar(x));
    p.parse(varargin{:});
    opts = p.Results;

    if strcmpi(opts.AIFMethod,'movmedian') && mod(opts.AIFWindow,2)==0
        opts.AIFWindow = opts.AIFWindow + 1;
    end
end


% ======================================================
% ===================== AIF builder ====================
% ======================================================

function [AIF_fit, meta] = makeAIF(t, AIF_raw, opts)
    method = lower(string(opts.AIFMethod));

    switch method
        case "movmedian"
            win  = opts.AIFWindow;
            AIF_fit  = movmedian(AIF_raw, win);
            meta = struct('method','movmedian','window',win);

        case "gamma"
            % Fit gamma, return both sampled and bin-averaged versions
            bAIF = opts.AIFBounds;
            bAIF(4,:) = bAIF(4,:) * max(AIF_raw(:));

            [AIF_samp, AIF_bin, pfit] = AIF_model_fit_dual(t, AIF_raw, bAIF);

            AIF_fit  = AIF_samp;  % point-sampled for modelling
            meta = struct('method','gamma', ...
                          'bounds', bAIF, ...
                          'params', pfit, ...
                          'AIF_sampled', AIF_samp, ...
                          'AIF_binned',  AIF_bin);

        case "gamma_smoothed"
            % Fit gamma, use centered TR-bin-averaged AIF_fit for modelling
            bAIF = opts.AIFBounds;
            bAIF(4,:) = bAIF(4,:) * max(AIF_raw(:));

            [AIF_samp, AIF_bin, pfit] = AIF_model_fit_dual(t, AIF_raw, bAIF);

            AIF_fit  = AIF_bin;   % bin-averaged for modelling
            meta = struct('method','gamma_smoothed', ...
                          'bounds', bAIF, ...
                          'params', pfit, ...
                          'AIF_sampled', AIF_samp, ...
                          'AIF_binned',  AIF_bin);

        case "none"
            AIF_fit  = AIF_raw;
            meta = struct('method','none');

        otherwise
            error('AIFMethod "%s" not recognized. Use "movmedian", "gamma", "gamma_smoothed", or "none".', method);
    end
end


function [AIF_sampled, AIF_binned, params_fit] = AIF_model_fit_dual(t, AIF_raw, bounds_AIF)
    % Fit shifted gamma variate to sampled AIF_fit, then compute:
    % 1) AIF_sampled: gamma evaluated at t
    % 2) AIF_binned: centered TR-bin-average of the gamma curve

    opts = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-16, ...
        'FunctionTolerance',1e-16, ...
        'Display','off');

    params0   = bounds_AIF(:,2);
    params_lb = bounds_AIF(:,1);
    params_ub = bounds_AIF(:,3);

    AIF_smooth0 = movmedian(AIF_raw, 3);

    % Fit on sampled time grid
    params_fit = lsqcurvefit(@(p,tt) AIF_model(p, tt), ...
        params0, t, AIF_smooth0(:)', params_lb, params_ub, opts);

    % Point-sampled version
    AIF_sampled = AIF_model(params_fit, t);

    % Estimate TR from time vector
    if numel(t) > 1
        TR = median(diff(t));
    else
        TR = 1;
    end

    % Oversampled continuous gamma
    oversamp = 10;
    dtfine = TR / oversamp;
    tfine = t(1):dtfine:(t(end) + TR);
    AIF_fine = AIF_model(params_fit, tfine);

    % Centered bin-average around each sample time
    AIF_binned = zeros(size(t));
    for i = 1:numel(t)
        t_start = t(i) - TR/2;
        t_end   = t(i) + TR/2;

        idx = (tfine >= t_start) & (tfine < t_end);
        if any(idx)
            AIF_binned(i) = mean(AIF_fine(idx));
        else
            AIF_binned(i) = AIF_sampled(i);
        end
    end
end


function AIF_fit = AIF_model(params, t)
    t0 = params(1); a = params(2); b = params(3); A = params(4);
    x  = t - t0;
    AIF_fit = A * (x > 0) .* (x).^(a-1) .* exp(-(x)/b);
end


% ======================================================
% ================ Fit metabolite mapping ==============
% ======================================================

function [sel_idx, fit_labels] = mapFitMetabolites(fitMets)

    if ischar(fitMets)
        fitMets = {fitMets};
    end

    labels = upper(string(fitMets(:)'));
    sel_idx = zeros(size(labels));

    for i = 1:numel(labels)
        switch labels(i)
            case {"L","LAC","LACTATE"}
                sel_idx(i) = 1;  % row 1 in model output
                labels(i) = "L";
            case {"P","PYR","PYRUVATE"}
                sel_idx(i) = 2;  % row 2
                labels(i) = "P";
            case {"B","BIC","BICARBONATE","HCO3"}
                sel_idx(i) = 3;  % row 3
                labels(i) = "B";
            case {"A","ALA","ALANINE"}
                sel_idx(i) = 4;  % row 4
                labels(i) = "A";
            otherwise
                error('FitMetabolites entry "%s" not recognized.', labels(i));
        end
    end

    if isempty(sel_idx)
        error('FitMetabolites must include at least one metabolite.');
    end

    fit_labels = cellstr(labels);
end


% ======================================================
% ================== R2 helper =========================
% ======================================================

function R2 = calcNormalizedR2(data_all, y_pred, sel_idx)
    Nt = size(data_all, 2);
    y_norm = [];
    y_pred_norm = [];

    for ii = 1:numel(sel_idx)
        r = sel_idx(ii);
        d = data_all(r, :);
        p = y_pred(r, :);

        npair = normalize([d, p], 'range')';
        y_norm      = [y_norm;      npair(1:Nt)]; %#ok<AGROW>
        y_pred_norm = [y_pred_norm; npair(Nt+1:end)]; %#ok<AGROW>
    end

    SStot = sum((y_norm(:) - mean(y_norm(:))).^2);
    SSres = sum((y_norm(:) - y_pred_norm(:)).^2);
    R2 = 1 - SSres / SStot;
end


% ======================================================
% ================== 4-met ODE model ===================
% ======================================================

function y = HEMEX_model_4met_singlekidney(params, t, AIF_fit, r2p, r2l, r2b, r2a, t0_override, scale_factor)
% Output rows: [L; P_total; B; A]
%
% params may be the full 17-parameter vector. Only the shared parameters,
% the requested t0 delay, and the requested kidney scale factor are used here.

    % Extract parameters
    kpl = params(1);
    klp = params(2);

    rp  = params(3) + r2p;
    rl  = params(4) * rp + r2l;
    
    E = params(5); % extract fraction
    
    t0  = t0_override;
    mu  = params(8);
    sg  = params(9) * mu; % As a function of mu: 0.1≲𝜎/𝜇≲1.5

    kpb = params(12);
    kbp = params(13);
    kpa = params(14);
    kap = params(15);

    rb  = params(16) * rp + r2b;
    ra  = params(17) * rp + r2a;

    Nt = numel(t);
    t  = t(:)';
    AIF_fit = AIF_fit(:)';

    % Terms for residue function
    alpha = mu^2 / sg^2;
    beta  = sg^2 / mu;

    k = ((1-E)^(-1/alpha)-1)/beta; 

    % Take into account the delay from AIF to given voxel/compartment
    AIF_shifted = interp1(t, AIF_fit, t - t0, 'linear', 0);

    % Model residue function
    R  = @(x) gammainc(x / beta, alpha, 'upper') .* (x >= 0);
    Rt = R(t);

    % Model vascular pyruvate
    Pv = (1/mu) * conv(AIF_shifted, exp(-(k + rp) * t) .* Rt, 'full');
    Pv = Pv(1:Nt);
    Pv = Pv(:);

    % Finally the flux of pyruvate from vascular to extravascular compartment
    J = k * Pv;

    % General matrix form of rate constants
    K = [ kpl + kpb + kpa + rp,  -klp,              -kbp,              -kap;
         -kpl,                   klp + rl,          0,                 0;
         -kpb,                   0,                 kbp + rb,          0;
         -kpa,                   0,                 0,                 kap + ra];

    nStates = 4;
    M = zeros(nStates, Nt);
    e1 = [1; 0; 0; 0]; % Only pyruvate gets vascular input

    % Assuming constant temporal spacing between dynamic frames
    if Nt > 1
        TR_local = t(2) - t(1);
    else
        TR_local = 1;
    end

    % ------------------------------------------------------
    % Solve linear compartment model:
    % dM/dt = -K*M + e1*J(t)
    %
    % M = [Pt; L; B; A]
    % ------------------------------------------------------

    % Exact kinetic evolution over one frame (matrix exponential of K)
    kinetic_propagator = expm(-K * TR_local);

    % Exact contribution from vascular pyruvate input over one frame
    vascular_input_operator = K \ ((eye(nStates) - kinetic_propagator) * e1);

    % Propagate metabolite states through time
    for i = 1:Nt-1
        M(:, i+1) = kinetic_propagator * M(:, i) + vascular_input_operator * J(i);
    end

    Pt = M(1, :);
    L  = M(2, :);
    B  = M(3, :);
    A  = M(4, :);

    y = zeros(4, Nt);
    y(1, :) = scale_factor * L;
    y(2, :) = scale_factor * (Pt + Pv.');
    y(3, :) = scale_factor * B;
    y(4, :) = scale_factor * A;
end
