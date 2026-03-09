function model_fit_results = HEMEX_main_4met(AIF_raw, pyr_norm, lac_norm, bic_norm, ala_norm, TR, flip_P, flip_L, flip_B, flip_A, varargin)
    % HEMEX_main_4met
    %   Fit the 4-metabolite HEMEX model (Pyr, Lac, Bic, Ala) using a matrix-exponential ODE solver,
    %   with flexible selection of which metabolites to include in the objective function.
    %   
    %   Usage:
    %     s = HEMEX_main_4met(AIF_raw, pyr, lac, bic, ala, TR, flipP, flipL, flipB, flipA);
    %     s = HEMEX_main_4met(..., 'FitMetabolites', {'P','L'});
    %     s = HEMEX_main_4met(..., 'FitMetabolites', {'P','L','B','A'});
    %     s = HEMEX_main_4met(..., 'AIFMethod','gamma', 'Bounds', myBounds);

    % -------------------- Parse options --------------------
    opts = parseOptions4met(varargin{:});

    % -------------------- Time vector ----------------------
    nTime = numel(AIF_raw);
    t = linspace(0, nTime*TR - TR, nTime);

    % -------------------- Pulse contributions --------------
    r2p = -log(cosd(flip_P)) / TR;
    r2l = -log(cosd(flip_L)) / TR;
    r2b = -log(cosd(flip_B)) / TR;
    r2a = -log(cosd(flip_A)) / TR;

    % -------------------- Build AIF_fit ------------------------
    [AIF_fit, aif_meta] = makeAIF(t, AIF_raw, opts);

    % -------------------- Prepare data ---------------------
    Pz = pyr_norm(:)';  % row
    Lz = lac_norm(:)';  % row
    Bz = bic_norm(:)';  % row
    Az = ala_norm(:)';  % row

    data_all = [Lz; Pz; Bz; Az];  % matches model output row order [L; P; B; A]

    % -------------------- Select metabolites to fit --------
    [sel_idx, fit_labels] = mapFitMetabolites(opts.FitMetabolites);

    data_sel = data_all(sel_idx, :);

    % 1/max weighting per selected curve
    maxVals = max(data_sel, [], 2);
    maxVals(maxVals <= 0) = 1;
    weights = 1 ./ maxVals;  % Nsel x 1

    target = data_sel .* weights;

    % -------------------- Bounds/init ----------------------
    params_lb = opts.Bounds(:,1);
    params0   = opts.Bounds(:,2);
    params_ub = opts.Bounds(:,3);

    options = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-16, 'FunctionTolerance',1e-16, ...
        'Display','off', 'MaxFunctionEvaluations',10000, 'MaxIterations',6000);

    % -------------------- Fit ------------------------------
    model_fun = @(p,tt) localModelSel(p, tt, AIF_fit, r2p, r2l, r2b, r2a, sel_idx, weights);

    [pars,~,~,exitflag] = lsqcurvefit( ...
        model_fun, params0, t, target, params_lb, params_ub, options);

    % -------------------- Predicted curves -----------------
    y_pred   = HEMEX_model_4met(pars, t, AIF_fit, r2p, r2l, r2b, r2a);
    lac_fits = y_pred(1,:);
    pyr_fits = y_pred(2,:);
    bic_fits = y_pred(3,:);
    ala_fits = y_pred(4,:);

    % -------------------- Normalized R^2 (selected mets) ---
    Nt = numel(t);
    y_norm = [];
    y_pred_norm = [];

    for ii = 1:numel(sel_idx)
        r = sel_idx(ii);
        d = data_all(r, :);
        p = y_pred(r, :);

        npair = normalize([d, p], 'range')';
        y_norm      = [y_norm;      npair(1:Nt)];
        y_pred_norm = [y_pred_norm; npair(Nt+1:end)];
    end

    SStot = sum((y_norm(:) - mean(y_norm(:))).^2);
    SSres = sum((y_norm(:) - y_pred_norm(:)).^2);
    R2    = 1 - SSres / SStot;

    % -------------------- Reparameterization ----------------
    kpl   = pars(1);
    klp   = pars(2);
    rp0   = pars(3);
    rl_sc = pars(4);
    k_sc  = pars(5);
    t0d   = pars(6);
    mu    = pars(7);
    sig   = pars(8);

    kpb   = pars(9);
    kbp   = pars(10);
    kpa   = pars(11);
    kap   = pars(12);
    rb_sc = pars(13);
    ra_sc = pars(14);

    rp_eff = rp0 + r2p;
    rl_eff = rl_sc * rp_eff + r2l;
    rb_eff = rb_sc * rp_eff + r2b;
    ra_eff = ra_sc * rp_eff + r2a;

    k_eff  = k_sc * kpl;

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
    model_fit_results.k_fits        = k_eff;
    model_fit_results.k_scale_fits  = k_sc;
    model_fit_results.t0_delay_fits = t0d;
    model_fit_results.mu_fits       = mu;
    model_fit_results.sigma_fits    = sig;

    model_fit_results.R2 = R2;

    % curves & bookkeeping
    model_fit_results.pyr_norm = pyr_norm;
    model_fit_results.lac_norm = lac_norm;
    model_fit_results.bic_norm = bic_norm;
    model_fit_results.ala_norm = ala_norm;

    model_fit_results.pyr_fits = pyr_fits;
    model_fit_results.lac_fits = lac_fits;
    model_fit_results.bic_fits = bic_fits;
    model_fit_results.ala_fits = ala_fits;

    model_fit_results.t        = t;
    model_fit_results.AIF_raw  = AIF_raw;
    model_fit_results.AIF_fit  = AIF_fit;
    model_fit_results.AIF_info = aif_meta;
end


% ======================================================
% ================== Local fit wrapper =================
% ======================================================

function ysel = localModelSel(p, tt, AIF_fit, r2p, r2l, r2b, r2a, sel_idx, weights)
    yfull = HEMEX_model_4met(p, tt, AIF_fit, r2p, r2l, r2b, r2a);
    ysel  = yfull(sel_idx, :) .* weights;
end


% ======================================================
% =================== Option parsing ===================
% ======================================================

function opts = parseOptions4met(varargin)

    defaultBounds = [ ...
        0.005  0.05    0.3;   % 1) kpl
        0      0       0;     % 2) klp
        0.01   1/30    0.05;  % 3) rp_fit
        0.8    1       1.2;   % 4) rl_scale
        0.05   1       20;    % 5) k_scale
       -5      0       10;    % 6) t0_delay
        1      5       30;    % 7) mu
        0.5    3       30;    % 8) sigma
        0      0.01    0.3;   % 9) kpb
        0      0       0;     % 10) kbp
        0      0.01    0.3;   % 11) kpa
        0      0       0;     % 12) kap
        0.8    1       1.2;   % 13) rb_scale
        0.8    1       1.2];  % 14) ra_scale

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
% ===================== AIF_fit builder ====================
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
    % 2) AIF_binned:  centered TR-bin-average of the gamma curve

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
% ================== 4-met ODE model ===================
% ======================================================

function y = HEMEX_model_4met(params, t, AIF_fit, r2p, r2l, r2b, r2a)
% Output rows: [L; P_total; B; A]

    kpl = params(1);
    klp = params(2);

    rp  = params(3) + r2p;
    rl  = params(4) * rp + r2l;

    k   = params(5) * kpl;
    t0  = params(6);
    mu  = params(7);
    sg  = params(8);

    kpb = params(9);
    kbp = params(10);
    kpa = params(11);
    kap = params(12);

    rb  = params(13) * rp + r2b;
    ra  = params(14) * rp + r2a;

    Nt = numel(t);
    t  = t(:)';     
    AIF_fit = AIF_fit(:)';

    alpha = mu^2 / sg^2;
    beta  = sg^2 / mu;

    AIF_shifted = interp1(t, AIF_fit, t - t0, 'linear', 0);

    R  = @(x) gammainc(x / beta, alpha, 'upper') .* (x >= 0);
    Rt = R(t);

    Pv = (1/mu) * conv(AIF_shifted, exp(-(k + rp) * t) .* Rt, 'full');
    Pv = Pv(1:Nt);
    Pv = Pv(:);

    J = k * Pv;

    K = [ kpl + kpb + kpa + rp,  -klp,              -kbp,              -kap;
         -kpl,                   klp + rl,          0,                 0;
         -kpb,                   0,                 kbp + rb,          0;
         -kpa,                   0,                 0,                 kap + ra ];

    nStates = 4;
    M = zeros(nStates, Nt);
    e1 = [1; 0; 0; 0];

    dt_vec = diff(t);

    if isempty(dt_vec)
        % return zeros for tissue, vascular term still contributes to P_total
        Pt = zeros(1, Nt);
        L  = zeros(1, Nt);
        B  = zeros(1, Nt);
        A  = zeros(1, Nt);

        y = zeros(4, Nt);
        y(1,:) = L;
        y(2,:) = Pt + Pv.';
        y(3,:) = B;
        y(4,:) = A;
        return;
    end

    if all(abs(dt_vec - dt_vec(1)) < 1e-12)
        dt = dt_vec(1);
        A_step = expm(-K * dt);
        B_step = K \ ((eye(nStates) - A_step) * e1);

        for i = 1:Nt-1
            M(:, i+1) = A_step * M(:, i) + B_step * J(i);
        end
    else
        for i = 1:Nt-1
            dt = dt_vec(i);
            if dt <= 0
                error('Time vector t must be strictly increasing.');
            end
            A_step = expm(-K * dt);
            B_step = K \ ((eye(nStates) - A_step) * e1);
            M(:, i+1) = A_step * M(:, i) + B_step * J(i);
        end
    end

    Pt = M(1, :);
    L  = M(2, :);
    B  = M(3, :);
    A  = M(4, :);

    y = zeros(4, Nt);
    y(1, :) = L;
    y(2, :) = Pt + Pv.';
    y(3, :) = B;
    y(4, :) = A;
end
