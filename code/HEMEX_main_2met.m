function model_fit_results = HEMEX_main_2met(AIF_raw, pyr_norm, lac_norm, TR, flip_P, flip_L, varargin)
    % HEMEX_main_2met
    %   Two-metabolite HEMEX model (Pyr, Lac) using a matrix-exponential ODE solver.
    %
    %   Usage:
    %     s = HEMEX_main_2met(AIF_raw, pyr, lac, TR, flipP, flipL);
    %     s = HEMEX_main_2met(..., 'AIFMethod','gamma', 'Bounds', myBounds);

    % -------------------- Parse options --------------------
    opts = parseOptions(varargin{:});

    % -------------------- Time vector ----------------------
    t = linspace(0, numel(AIF_raw)*TR - TR, numel(AIF_raw));

    % -------------------- Pulse contributions (Hill 2013) --
    r2p = -log(cosd(flip_P)) / TR;
    r2l = -log(cosd(flip_L)) / TR;

    % -------------------- Build AIF ------------------------
    [AIF, aif_meta] = makeAIF(t, AIF_raw, opts);

    % -------------------- Prepare fit ----------------------
    Pz = pyr_norm(:)';  % row
    Lz = lac_norm(:)';  % row

    % balance pyr/lac magnitudes; simple 1/max weighting
    AL = max(Lz); AP = max(Pz);
    if AL <= 0, AL = 1; end
    if AP <= 0, AP = 1; end
    weights = 1 ./ [AL; AP];   % 2x1

    params_lb = opts.Bounds(:,1);
    params0   = opts.Bounds(:,2);
    params_ub = opts.Bounds(:,3);

    options = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-16, 'FunctionTolerance',1e-16, ...
        'Display','off', 'MaxFunctionEvaluations',10000, 'MaxIterations',6000);

    % -------------- Fit (two curves stacked) ---------------
    model_fun = @(p,tt) HEMEX_model_2met_expm(p, tt, AIF, r2p, r2l) .* weights;
    target    = [Lz; Pz] .* weights;

    [pars,resnorm,residual,exitflag,output] = lsqcurvefit( ...
        model_fun, params0, t, target, params_lb, params_ub, options); %#ok<ASGLU>

    % -------------- Predicted curves -----------------------
    y_pred   = HEMEX_model_2met_expm(pars, t, AIF, r2p, r2l);
    lac_fits = y_pred(1,:);
    pyr_fits = y_pred(2,:);

    % -------------- Normalized R^2 (pyr+lac together) -----
    timepoints  = numel(Pz);
    norm_lac    = normalize([Lz, y_pred(1,:)], 'range')';
    norm_pyr    = normalize([Pz, y_pred(2,:)], 'range')';
    y_norm      = [norm_lac(1:timepoints); norm_pyr(1:timepoints)];
    y_pred_norm = [norm_lac(timepoints+1:end); norm_pyr(timepoints+1:end)];
    SStot = sum((y_norm(:) - mean(y_norm(:))).^2);
    SSres = sum((y_norm(:) - y_pred_norm(:)).^2);
    R2    = 1 - SSres / SStot;

    % -------------- Reparameterization ---------------------
    kpl   = pars(1);
    klp   = pars(2);
    rp    = pars(3);
    rl_sc = pars(4);     % rl = rl_sc * rp
    k_sc  = pars(5);     % k = k_sc  * kpl
    t0d   = pars(6);
    mu    = pars(7);
    sig   = pars(8);

    % -------------- Package results ------------------------
    model_fit_results.par_fits      = pars;
    model_fit_results.exitflag_fits = exitflag;
    model_fit_results.kpl_fits      = kpl;
    model_fit_results.klp_fits      = klp;
    model_fit_results.r1p_fits      = rp;
    model_fit_results.r1l_fits      = rl_sc * rp;      % rl derived from rp
    model_fit_results.k_fits        = k_sc  * kpl;     % k  derived from kpl
    model_fit_results.t0_delay_fits = t0d;
    model_fit_results.mu_fits       = mu;
    model_fit_results.sigma_fits    = sig;
    model_fit_results.R2            = R2;

    % curves & bookkeeping
    model_fit_results.pyr_norm = pyr_norm;
    model_fit_results.lac_norm = lac_norm;
    model_fit_results.pyr_fits = pyr_fits;
    model_fit_results.lac_fits = lac_fits;
    model_fit_results.t        = t;
    model_fit_results.AIF_raw  = AIF_raw;
    model_fit_results.AIF      = AIF;
    model_fit_results.AIF_info = aif_meta;   % struct: method, window/params/bounds
end


% ======================================================
% ================ Helper functions ====================
% ======================================================

function opts = parseOptions(varargin)
    % Same defaults as your original HEMEX_main.m
    defaultBounds = [ ...
        0.005  0.05    0.3;    % 1) kpl
        0      0       0;      % 2) klp
        0.01   1/30    0.05;   % 3) rp
        0.8    1       1.2;    % 4) rl scale (rl = rl_sc * rp)
        0.05   1       20;     % 5) k  scale (k  = k_sc  * kpl)
       -5      0       10;     % 6) t0_delay
        1      5       30;     % 7) mu
        0.5    3       30];    % 8) sigma

    % Default AIF gamma bounds (params = [t0; a; b; A])
    defaultAIFBounds = [ ...
       -5   2   20;   % t0
        1   3    5;   % a
        1   1.5  3;   % b
        0.2 1    1.5];% A (will be scaled by max(AIF_raw) inside makeAIF)

    p = inputParser;
    p.addParameter('Bounds',     defaultBounds,     @(x) isnumeric(x) && isequal(size(x,2),3));
    p.addParameter('AIFMethod',  'movmedian',       @(s) ischar(s) || isstring(s));
    p.addParameter('AIFWindow',  3,                 @(n) isnumeric(n) && isscalar(n) && n>0);
    p.addParameter('AIFBounds',  defaultAIFBounds,  @(x) isnumeric(x) && isequal(size(x),[4 3]));
    p.parse(varargin{:});
    opts = p.Results;

    % ensure movmedian window is odd
    if strcmpi(opts.AIFMethod,'movmedian') && mod(opts.AIFWindow,2)==0
        opts.AIFWindow = opts.AIFWindow + 1;
    end
end


function [AIF, meta] = makeAIF(t, AIF_raw, opts)
    method = lower(string(opts.AIFMethod));
    switch method
        case "movmedian"
            win  = opts.AIFWindow;
            AIF  = movmedian(AIF_raw, win);
            meta = struct('method','movmedian','window',win);

        case "gamma"
            % scale amplitude row to data scale
            bAIF = opts.AIFBounds;
            bAIF(4,:) = bAIF(4,:) * max(AIF_raw(:));   % scale A row
            [AIF, pfit] = AIF_model_fit(t, AIF_raw, bAIF);
            meta = struct('method','gamma','bounds',bAIF,'params',pfit);

        case "none"
            AIF  = AIF_raw;
            meta = struct('method','none');

        otherwise
            error('AIFMethod "%s" not recognized. Use "movmedian", "gamma", or "none".', method);
    end
end


function [AIF_fit, params_fit] = AIF_model_fit(t, AIF_raw, bounds_AIF)
    % Fit shifted gamma variate:
    % params = [t0; a; b; A]
    opts = optimoptions('lsqcurvefit', ...
        'OptimalityTolerance',1e-16, 'FunctionTolerance',1e-16, 'Display','off');

    params0   = bounds_AIF(:,2);
    params_lb = bounds_AIF(:,1);
    params_ub = bounds_AIF(:,3);

    % gentle pre-smoothing stabilizes fit
    AIF_smooth0 = movmedian(AIF_raw, 3);

    params_fit = lsqcurvefit(@(p,tt) AIF_model(p, tt), params0, t, AIF_smooth0(:)', ...
                             params_lb, params_ub, opts);

    AIF_fit = AIF_model(params_fit, t);
end


function AIF = AIF_model(params, t)
    % Shifted gamma variate: params = [t0; a; b; A]
    t0 = params(1); a = params(2); b = params(3); A = params(4);
    x  = t - t0;
    AIF = A * (x > 0) .* (x).^(a-1) .* exp(-(x)/b);
    % (x.^0 for x==0 is handled by (x>0) gate)
end


% ======================================================
% ================ 2-met ODE model =====================
% ======================================================

function y = HEMEX_model_2met_expm(params, t, AIF, r2p, r2l)
    % HEMEX_model_2met_expm
    %   Two-metabolite HEMEX model (Pyr, Lac) using the same
    %   matrix-exponential time stepping as the 4-met model.
    %
    %   y(1,:) = L(t)          (tissue lactate)
    %   y(2,:) = P_total(t)    (tissue + vascular pyruvate)

    % --- unpack parameters ---
    kpl = params(1);
    klp = params(2);

    rp  = params(3) + r2p;        % effective P relax (incl. pulsing)
    rl  = params(4) * rp + r2l;   % L relax

    k   = params(5) * kpl;        % extraction rate
    t0  = params(6);
    mu  = params(7);
    sg  = params(8);

    Nt = numel(t);
    t  = t(:)';                   % row
    AIF = AIF(:)';                % row

    % --- gamma-distributed residue function R(t) ---
    alpha = mu^2 / sg^2;
    beta  = sg^2 / mu;

    % time-shifted AIF
    AIF_shifted = interp1(t, AIF, t - t0, 'linear', 0);

    % residue
    R = @(x) gammainc(x / beta, alpha, 'upper') .* (x >= 0);
    Rt = R(t);

    % --- vascular pyruvate Pv(t) ---
    Pv = (1/mu) * conv(AIF_shifted, exp(-(k + rp) * t) .* Rt, 'full');
    Pv = Pv(1:Nt);               % truncate
    Pv = Pv(:);                  % column

    % flux into tissue pyruvate
    J = k * Pv;                  % column vector length Nt

    % --- build K matrix for [P, L]^T ---
    K = [ kpl + rp,  -klp;
         -kpl,       klp + rl ];

    nStates = 2;
    M = zeros(nStates, Nt);      % [P; L], tissue only
    e1 = [1; 0];                 % forcing only on P

    dt_vec = diff(t);

    % check if time step is (approximately) constant
    if all(abs(dt_vec - dt_vec(1)) < 1e-12)
        dt = dt_vec(1);
        A_step = expm(-K * dt);
        B_step = K \ ((eye(nStates) - A_step) * e1);

        for i = 1:Nt-1
            M(:, i+1) = A_step * M(:, i) + B_step * J(i);
        end
    else
        % general case: varying dt
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

    % unpack tissue components
    Pt = M(1, :);   % tissue pyruvate
    L  = M(2, :);   % lactate

    % --- construct observable signals ---
    y = zeros(2, Nt);
    y(1, :) = L;          % lactate (tissue)
    y(2, :) = Pt + Pv.';  % pyruvate (tissue + vascular)
end
