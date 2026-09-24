function model_fit_results = HEMEX_main_4met_liver( ...
    input_functions_raw, pyr_norm, lac_norm, bic_norm, ala_norm, ...
    TR, flip_P, flip_L, flip_B, flip_A, varargin)
    % ======================================================
    % ============ 4-met Dual-input liver HEMEX ============
    % ======================================================
    %
    % INPUTS
    %   input_functions_raw : N x 2 matrix
    %       col 1 = combined pyruvate vascular input (arterial + portal)
    %       col 2 = portal lactate vascular input
    %
    %   pyr_norm, lac_norm, bic_norm, ala_norm : measured metabolite curves
    %   TR : temporal spacing of the dynamic data
    %   flip_P, flip_L, flip_B, flip_A : flip angles in degrees
    %
    % OPTIONAL
    %   'Bounds'           : bounds for the tissue fit
    %   'InputBoundsPyr'   : bounds for dual-gamma pyruvate fit
    %   'InputBoundsLac'   : bounds for single-gamma lactate fit
    %   'FitMetabolites'   : e.g. {'P','L'} or {'P','L','B','A'}
    %
    % OUTPUT
    %   model_fit_results : struct with fit parameters, curves, and metadata

    % -------------------- Parse options --------------------
    opts = parseOptions4metLiver(varargin{:});

    % -------------------- Basic checks --------------------
    if size(input_functions_raw, 2) ~= 2
        error('input_functions_raw must be an N x 2 matrix: [combined pyruvate, lactate].');
    end

    nTime = size(input_functions_raw, 1);
    if nTime < 2
        error('Need at least two time points.');
    end

    % -------------------- Time vector ----------------------
    t = (0:nTime-1) * TR;
    t = t(:)';

    % -------------------- Pulse contributions --------------
    r2p = -log(cosd(flip_P)) / TR;
    r2l = -log(cosd(flip_L)) / TR;
    r2b = -log(cosd(flip_B)) / TR;
    r2a = -log(cosd(flip_A)) / TR;

    % -------------------- Fit vascular input functions ------
    [input_functions_fit, input_fit_meta] = makeInputFunctions_liver(t, input_functions_raw, opts);

    % Split fitted vascular curves for convenience
    MPA_fit = input_functions_fit(:,1)';   % arterial pyruvate
    MPV_fit = input_functions_fit(:,2)';   % portal pyruvate
    MLV_fit = input_functions_fit(:,3)';   % portal lactate

    MP_total_raw = input_functions_raw(:,1)';   % measured combined pyruvate
    MLV_raw      = input_functions_raw(:,2)';    % measured portal lactate

    % -------------------- Prepare metabolite data -----------
    Pz = pyr_norm(:)';  % row
    Lz = lac_norm(:)';  % row
    Bz = bic_norm(:)';  % row
    Az = ala_norm(:)';  % row

    % Model output row order is [L; P; B; A]
    data_all = [Lz; Pz; Bz; Az];

    % -------------------- Select metabolites to fit ---------
    [sel_idx, fit_labels] = mapFitMetabolites(opts.FitMetabolites);
    data_sel = data_all(sel_idx, :);

    % Weight each selected curve by its max value
    maxVals = max(data_sel, [], 2);
    maxVals(maxVals <= 0) = 1;
    weights = 1 ./ maxVals;
    target  = data_sel .* weights;

    % -------------------- Bounds / initial guess ------------
    params_lb = opts.Bounds(:,1);
    params0   = opts.Bounds(:,2);
    params_ub = opts.Bounds(:,3);

    options = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-12, ...
        'FunctionTolerance',1e-12, ...
        'Display','off', ...
        'MaxFunctionEvaluations',10000, ...
        'MaxIterations',6000);

    % -------------------- Fit tissue model ------------------
    model_fun = @(p,tt) localModelSelLiver( ...
        p, tt, input_functions_fit, r2p, r2l, r2b, r2a, sel_idx, weights);

    [pars,~,~,exitflag] = lsqcurvefit( ...
        model_fun, params0, t, target, params_lb, params_ub, options);

    % -------------------- Predicted curves ------------------
    y_pred = HEMEX_model_4met_liver(pars, t, input_functions_fit, r2p, r2l, r2b, r2a);

    lac_fits = y_pred(1,:);
    pyr_fits = y_pred(2,:);
    bic_fits = y_pred(3,:);
    ala_fits = y_pred(4,:);

    % -------------------- Normalized R^2 (selected mets) ----
    y_norm = [];
    y_pred_norm = [];

    for ii = 1:numel(sel_idx)
        r = sel_idx(ii);
        d = data_all(r, :);
        p = y_pred(r, :);

        npair = normalize([d, p], 'range')';
        y_norm      = [y_norm;      npair(1:nTime)];
        y_pred_norm = [y_pred_norm; npair(nTime+1:end)];
    end

    SStot = sum((y_norm(:) - mean(y_norm(:))).^2);
    SSres = sum((y_norm(:) - y_pred_norm(:)).^2);
    if SStot > 0
        R2 = 1 - SSres / SStot;
    else
        R2 = NaN;
    end

    % -------------------- Reparameterization ----------------
    kpl   = pars(1);
    klp   = pars(2);
    rp0   = pars(3);
    rl_sc = pars(4);

    kP    = pars(5);
    kL    = pars(6);
    wA    = pars(7);
    t0A   = pars(8);
    t0V   = pars(9);
    mu    = pars(10);
    sig   = pars(11);

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

    % -------------------- Package results -------------------
    model_fit_results.par_fits       = pars;
    model_fit_results.exitflag_fits  = exitflag;

    model_fit_results.fit_metabolites = fit_labels;
    model_fit_results.fit_rows        = sel_idx;

    % exchange rates
    model_fit_results.kpl_fits = kpl;
    model_fit_results.klp_fits = klp;
    model_fit_results.kP_fits  = kP;
    model_fit_results.kL_fits  = kL;
    model_fit_results.kpb_fits  = kpb;
    model_fit_results.kbp_fits  = kbp;
    model_fit_results.kpa_fits  = kpa;
    model_fit_results.kap_fits  = kap;

    % vascular fraction
    model_fit_results.wA_fits = wA;
    model_fit_results.wV_fits = 1 - wA;

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

    % residue / hemodynamics
    model_fit_results.mu_fits    = mu;
    model_fit_results.sigma_fits = sig;

    model_fit_results.R2 = R2;

    % measured curves
    model_fit_results.pyr_norm = pyr_norm;
    model_fit_results.lac_norm = lac_norm;
    model_fit_results.bic_norm = bic_norm;
    model_fit_results.ala_norm = ala_norm;

    % fitted metabolite curves
    model_fit_results.pyr_fits = pyr_fits;
    model_fit_results.lac_fits = lac_fits;
    model_fit_results.bic_fits = bic_fits;
    model_fit_results.ala_fits = ala_fits;

    % time vector
    model_fit_results.t = t;

    % vascular input bookkeeping
    model_fit_results.input_functions_raw = input_functions_raw;
    model_fit_results.MP_total_raw = MP_total_raw;
    model_fit_results.MLV_raw = MLV_raw;

    model_fit_results.input_functions_fit = input_functions_fit;
    model_fit_results.MPA_fit = MPA_fit;
    model_fit_results.MPV_fit = MPV_fit;
    model_fit_results.MLV_fit = MLV_fit;
    model_fit_results.MP_total_fit = MPA_fit + MPV_fit;

    model_fit_results.input_fit_meta = input_fit_meta;
end


% ======================================================
% ============== Local fit wrapper =====================
% ======================================================

function ysel = localModelSelLiver(p, tt, input_functions_fit, r2p, r2l, r2b, r2a, sel_idx, weights)
    yfull = HEMEX_model_4met_liver(p, tt, input_functions_fit, r2p, r2l, r2b, r2a);
    ysel  = yfull(sel_idx, :) .* weights;
end


% ======================================================
% =================== Option parsing ===================
% ======================================================

function opts = parseOptions4metLiver(varargin)

    defaultBounds = [ ...
        0.005  0.05   0.3;   % 1)  kpl
        0      0      0;     % 2)  klp
        0.01   1/30   0.05;  % 3)  rp_fit
        0.8    1      1.2;   % 4)  rl_scale
        0.02   0.2    5;     % 5)  kP
        0.02   0.2    5;     % 6)  kL
        0      0.5    1;     % 7)  wA
        -5     10     30;    % 8)  t0A
        5      15     40;    % 9)  t0V
        1      5      30;    % 10) mu
        0.5    3      30;    % 11) sigma
        0      0.01   0.3;   % 12) kpb
        0      0      0;     % 13) kbp
        0      0.01   0.3;   % 14) kpa
        0      0      0;     % 15) kap
        0.8    1      1.2;   % 16) rb_scale
        0.8    1      1.2];  % 17) ra_scale

    % Dual gamma fit for combined pyruvate input:
    % [t0A, aA, bA, AA, t0V, aV, bV, AV]
    defaultInputBoundsPyr = [ ...
        0    5    20;   % t0A
        1.2  3    8;    % aA
        0.5  1.5  8;    % bA
        0.05 0.25 2;    % AA  (scaled by max input)
        0    12   40;   % t0V
        1.2  4    12;   % aV
        1    4    20;   % bV
        0.05 0.20 2];   % AV  (scaled by max input)

    % Single gamma fit for lactate input:
    % [t0, a, b, A]
    defaultInputBoundsLac = [ ...
        0    10   40;   % t0
        1.2  4    12;   % a
        0.5  4    20;   % b
        0.05 0.20 2];   % A  (scaled by max input)

    defaultFitMets = {'P','L','B','A'};

    p = inputParser;
    p.addParameter('Bounds',         defaultBounds,       @(x) isnumeric(x) && isequal(size(x,2),3));
    p.addParameter('InputBoundsPyr', defaultInputBoundsPyr, @(x) isnumeric(x) && isequal(size(x),[8 3]));
    p.addParameter('InputBoundsLac', defaultInputBoundsLac, @(x) isnumeric(x) && isequal(size(x),[4 3]));
    p.addParameter('FitMetabolites', defaultFitMets,      @(x) iscell(x) || isstring(x) || ischar(x));
    p.parse(varargin{:});
    opts = p.Results;
end


% ======================================================
% ============== Input function fitting =================
% ======================================================

function [input_functions_fit, meta] = makeInputFunctions_liver(t, input_functions_raw, opts)

    MP_total_raw = input_functions_raw(:,1)';
    MLV_raw      = input_functions_raw(:,2)';

    % Scale amplitude bounds by the data maxima
    bP = opts.InputBoundsPyr;
    bL = opts.InputBoundsLac;

    maxP = max(MP_total_raw);
    if maxP <= 0
        maxP = 1;
    end

    maxL = max(MLV_raw);
    if maxL <= 0
        maxL = 1;
    end

    bP([4 8], :) = bP([4 8], :) * maxP;
    bL(4, :)     = bL(4, :) * maxL;

    [MPA_fit, MPV_fit, pyr_meta] = fitDualGammaVariate(t, MP_total_raw, bP);
    [MLV_fit, lac_meta]          = fitSingleGammaVariate(t, MLV_raw, bL);

    input_functions_fit = [MPA_fit(:), MPV_fit(:), MLV_fit(:)];

    meta = struct();
    meta.pyruvate = pyr_meta;
    meta.lactate  = lac_meta;

    meta.MP_total_raw = MP_total_raw;
    meta.MLV_raw      = MLV_raw;

    meta.MPA_fit      = MPA_fit;
    meta.MPV_fit      = MPV_fit;
    meta.MLV_fit      = MLV_fit;
    meta.MP_total_fit = MPA_fit + MPV_fit;

    meta.input_functions_fit = input_functions_fit;
end


function [MPA_fit, MPV_fit, meta] = fitDualGammaVariate(t, raw, boundsPyr)

    opts = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-12, ...
        'FunctionTolerance',1e-12, ...
        'Display','off', ...
        'MaxFunctionEvaluations',10000, ...
        'MaxIterations',6000);

    params0   = boundsPyr(:,2);
    params_lb = boundsPyr(:,1);
    params_ub = boundsPyr(:,3);

    pfit = lsqcurvefit(@(p,tt) dualGammaModel(p, tt), ...
        params0, t, raw(:)', params_lb, params_ub, opts);

    pA = pfit(1:4);
    pV = pfit(5:8);

    MPA_fit = gammaVariate(pA, t);
    MPV_fit = gammaVariate(pV, t);

    % Ensure the arterial component is the earlier one
    [~, idxA] = max(MPA_fit);
    [~, idxV] = max(MPV_fit);
    if t(idxV) < t(idxA)
        tmp = pA;
        pA  = pV;
        pV  = tmp;

        MPA_fit = gammaVariate(pA, t);
        MPV_fit = gammaVariate(pV, t);
        pfit = [pA(:); pV(:)];
    end

    meta = struct();
    meta.params_fit = pfit;
    meta.params_A    = pA;
    meta.params_V    = pV;
    meta.MPA_fit     = MPA_fit;
    meta.MPV_fit     = MPV_fit;
    meta.MP_total_fit = MPA_fit + MPV_fit;
    meta.raw         = raw;
end


function [MLV_fit, meta] = fitSingleGammaVariate(t, raw, boundsLac)

    opts = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-12, ...
        'FunctionTolerance',1e-12, ...
        'Display','off', ...
        'MaxFunctionEvaluations',10000, ...
        'MaxIterations',6000);

    params0   = boundsLac(:,2);
    params_lb = boundsLac(:,1);
    params_ub = boundsLac(:,3);

    pfit = lsqcurvefit(@(p,tt) gammaVariate(p, tt), ...
        params0, t, raw(:)', params_lb, params_ub, opts);

    MLV_fit = gammaVariate(pfit, t);

    meta = struct();
    meta.params_fit = pfit;
    meta.MLV_fit    = MLV_fit;
    meta.raw        = raw;
end


function y = dualGammaModel(p, t)
    y = gammaVariate(p(1:4), t) + gammaVariate(p(5:8), t);
end


function y = gammaVariate(params, t)
    % params = [t0, a, b, A]
    t0 = params(1);
    a  = params(2);
    b  = params(3);
    A  = params(4);

    x = t - t0;
    y = A .* (x > 0) .* (x .^ (a - 1)) .* exp(-x ./ b);
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
                sel_idx(i) = 1;
                labels(i) = "L";
            case {"P","PYR","PYRUVATE"}
                sel_idx(i) = 2;
                labels(i) = "P";
            case {"B","BIC","BICARBONATE","HCO3"}
                sel_idx(i) = 3;
                labels(i) = "B";
            case {"A","ALA","ALANINE"}
                sel_idx(i) = 4;
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


function y = HEMEX_model_4met_liver(params, t, input_functions_fit, r2p, r2l, r2b, r2a)
% ======================================================
% =========== 4-met Dual-input liver HEMEX ============
% ======================================================
%
% vascular_inputs columns:
%   1 = MPA_fit : arterial pyruvate vascular input
%   2 = MPV_fit : portal pyruvate vascular input
%   3 = MLV_fit : portal lactate vascular input
%
% Output rows:
%   [L_total; P_total; B; A]

    % ---------------- Extract vascular inputs ----------------
    MPA_fit = input_functions_fit(:,1);
    MPV_fit = input_functions_fit(:,2);
    MLV_fit = input_functions_fit(:,3);

    % ---------------- Extract parameters ----------------
    % Metabolic exchange
    kpl = params(1);
    klp = params(2);

    % Tissue relaxation
    rp  = params(3) + r2p;
    rl  = params(4) * rp + r2l;

    % Vascular exchange rates
    kP  = params(5);   % pyruvate vascular exchange
    kL  = params(6);   % lactate vascular exchange

    % Arterial fraction (portal = 1-wA)
    wA  = params(7);
    wV  = 1 - wA;

    % Shared vascular residue parameters
    t0A = params(8);   % arterial delay
    t0V = params(9);   % portal delay
    mu  = params(10);
    sg  = params(11);

    % Additional pathways
    kpb = params(12);
    kbp = params(13);
    kpa = params(14);
    kap = params(15);

    % Relaxation scaling
    rb  = params(16) * rp + r2b;
    ra  = params(17) * rp + r2a;

    % ---------------- Setup ----------------
    Nt = numel(t);
    t  = t(:)';
    MPA_fit = MPA_fit(:)';
    MPV_fit = MPV_fit(:)';
    MLV_fit = MLV_fit(:)';

    alpha = mu^2 / sg^2;
    beta  = sg^2 / mu;

    % ---------------- Delayed vascular inputs ----------------
    MPA_shifted = interp1(t, MPA_fit, t - t0A, 'linear', 0);
    MPV_shifted = interp1(t, MPV_fit, t - t0V, 'linear', 0);
    MLV_shifted = interp1(t, MLV_fit, t - t0V, 'linear', 0);

    % ---------------- Shared vascular residue ----------------
    R  = @(x) gammainc(x / beta, alpha, 'upper') .* (x >= 0);
    Rt = R(t);

    % ---------------- Microvascular pyruvate ----------------
    PvA = (1/mu) * conv(MPA_shifted, exp(-(kP + rp) * t) .* Rt, 'full');
    PvA = PvA(1:Nt);

    PvV = (1/mu) * conv(MPV_shifted, exp(-(kP + rp) * t) .* Rt, 'full');
    PvV = PvV(1:Nt);

    % Weighted dual-input pyruvate vascular signal
    Pv = wA * PvA + wV * PvV;

    % ---------------- Microvascular lactate ----------------
    LvV = (1/mu) * conv(MLV_shifted, exp(-(kL + rl) * t) .* Rt, 'full');
    LvV = LvV(1:Nt);

    % ---------------- Vascular exchange currents ----------------
    JP = kP * Pv;
    JL = kL * LvV;

    % ---------------- Tissue kinetic matrix ----------------
    K = [ kpl + kpb + kpa + rp,  -klp,              -kbp,              -kap;
         -kpl,                   klp + rl,          0,                 0;
         -kpb,                   0,                 kbp + rb,          0;
         -kpa,                   0,                 0,                 kap + ra];

    nStates = 4;
    M = zeros(nStates, Nt);

    % Pyruvate enters Pt, lactate enters Lt
    e1 = [1; 0; 0; 0];
    e2 = [0; 1; 0; 0];
    U  = [e1, e2];

    % ---------------- ODE solve ----------------
    TR = t(2) - t(1);

    % Exact kinetic evolution over one frame
    kinetic_propagator = expm(-K * TR);

    % Exact vascular pyruvate + lactate input contribution
    vascular_input_operator = K \ ((eye(nStates) - kinetic_propagator) * U);

    % Propagate metabolite states through time
    for i = 1:Nt-1
        Jvec = [JP(i); JL(i)];
        M(:, i+1) = kinetic_propagator * M(:, i) + vascular_input_operator * Jvec;
    end

    % ---------------- Extract tissue states ----------------
    Pt = M(1,:);
    Lt  = M(2,:);
    Bt  = M(3,:);
    At  = M(4,:);

    % ---------------- Measured signals ----------------
    y = zeros(4, Nt);

    y(1,:) = Lt + LvV;   % total lactate = tissue + vascular portal lactate
    y(2,:) = Pt + Pv;   % total pyruvate = tissue + vascular dual-input
    y(3,:) = Bt; % total bicarbonate = tissue
    y(4,:) = At; % total alanine = tissue
end
