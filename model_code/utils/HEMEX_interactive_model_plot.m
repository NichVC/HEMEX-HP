% @author: Nichlas Vous Christensen
% @email: nvc@clin.au.dk
% @phone: +45 23464522
% @organization: Aarhus University, The MR Research Centre
% June 2025

function HEMEX_interactive_model_plot(t, AIF, pyr_data, lac_data, params, r2p, r2l, slider_min, slider_max)
    % Store the initial parameters
    % Note that these are not using the reparameterization of rl and k
    initial_params = params;

    % Create figure
    fig = figure('Name', 'Interactive Model Plot', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.6]);

    % Create plot panel
    plotPanel = uipanel('Parent', fig, 'Position', [0.05, 0.4, 0.9, 0.55]);
    ax = axes('Parent', plotPanel);
    
    [y_pred, ~] = HEMEX_model(params, t, AIF, r2p, r2l);
    
    % Data points
    plot(ax, t, pyr_data, 'ro'); % Static Pyruvate data points
    hold on;
    plotHandle2 = plot(ax, t, y_pred(2, :), 'r'); % Pyruvate
    plot(ax, t, lac_data, 'bo'); % Static Lactate data points    
    plotHandle1 = plot(ax, t, y_pred(1, :), 'b'); % Lactate
    
    hold off;
    legend(ax,'Pyruvate','Pyruvate fit','Lactate','Lactate fit')
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Signal');
    grid(ax, 'on');
    updateR2()

    % Slider properties
    slider_names = {'kpl', 'klp', 'rp', 'rl', 'k', 't0', 'mu', 'sigma'};
    num_sliders = length(slider_names);
    sliders = gobjects(num_sliders, 1);
    sliderTexts = gobjects(num_sliders, 1);

    % Create sliders panel
    slidersPanel = uipanel('Parent', fig, 'Position', [0.05, 0.05, 0.9, 0.3]);

    for i = 1:num_sliders
        row = floor((i - 1) / 5) + 1;
        col = mod((i - 1), 5) + 1;
        sliderTexts(i) = uicontrol('Parent', slidersPanel, 'Style', 'text', ...
            'Units', 'normalized', 'Position', [0.02 + (col - 1) * 0.19, 0.7 - (row - 1) * 0.5, 0.17, 0.2], ...
            'String', sprintf('%s: %.2f', slider_names{i}, params(i)));
        sliders(i) = uicontrol('Parent', slidersPanel, 'Style', 'slider', ...
            'Min', slider_min(i), 'Max', slider_max(i), 'Value', params(i), ...
            'Units', 'normalized', 'Position', [0.02 + (col - 1) * 0.19, 0.6 - (row - 1) * 0.5, 0.17, 0.2], ...
            'Callback', @updatePlot);
    end

    % Add reset button
    resetButton = uicontrol('Parent', slidersPanel, 'Style', 'pushbutton', 'String', 'Reset to Fit', ...
        'Units', 'normalized', 'Position', [0.8, 0.1, 0.15, 0.2], 'Callback', @resetParams);

    % Callback function to update the plot
    function updatePlot(~, ~)
        for i = 1:num_sliders
            params(i) = sliders(i).Value;
            sliderTexts(i).String = sprintf('%s: %.2f', slider_names{i}, params(i));
        end
        [y_pred, ~] = HEMEX_model(params, t, AIF, r2p, r2l);
        set(plotHandle1, 'XData', t, 'YData', y_pred(1, :));
        set(plotHandle2, 'XData', t, 'YData', y_pred(2, :));     
        updateR2()
    end

    function updateR2(~, ~)
        % First normalize the lac+lac_fit and pyr+pyr_fit together (to between [0,1])
        timepoints = length(pyr_data);
        norm_lac = normalize([lac_data',y_pred(1,:)],'range')';
        norm_pyr = normalize([pyr_data',y_pred(2,:)],'range')';
        y_norm = [norm_lac(1:timepoints);norm_pyr(1:timepoints)];
        y_pred_norm = [norm_lac(timepoints+1:end);norm_pyr(timepoints+1:end)];
        % Calculate R2 and display it (equally weighted).
        SStot = sum((y_norm(:)-mean(y_norm(:))).^2); % Total Sum-Of-Squares
        SSres = sum((y_norm(:)-y_pred_norm(:)).^2); % Residual Sum-Of-Squares
        R2 = 1-SSres/SStot;   
        title(ax,['R2 = ',num2str(R2),' (combined and equally weighted)']);
    end

    % Callback function to reset parameters
    function resetParams(~, ~)
        params = initial_params;
        for i = 1:num_sliders
            set(sliders(i), 'Value', params(i));
            sliderTexts(i).String = sprintf('%s: %.2f', slider_names{i}, params(i));
        end
        [y_pred, ~] = HEMEX_model(params, t, AIF, r2p, r2l);
        set(plotHandle1, 'XData', t, 'YData', y_pred(1, :));
        set(plotHandle2, 'XData', t, 'YData', y_pred(2, :));
        updateR2()
    end
end

function [y,Pv] = HEMEX_model(params, t, AIF, r2p, r2l)
    kpl = params(1);
    klp = params(2);
    rp = params(3) + r2p; % Take relaxation due to pulsing into account
    rl = params(4) + r2l; % Take relaxation due to pulsing into account
    k = params(5);
    t0_delay = params(6);
    mu = params(7);
    sigma = params(8);
    
    alpha = mu^2 ./ sigma^2;
    beta = sigma^2 ./ mu;

    AIF_shifted = interp1(t, AIF, t - t0_delay, 'linear', 0); % Delay between AIF and given voxel

    R = @(t)gammainc(t / beta, alpha, 'upper') .* (t >= 0);
    Pv = (1/(alpha*beta)) * conv(AIF_shifted, exp(-(k + rp) * t) .* R(t), 'full');
    Pv = Pv(1:length(t));

    k_plus = -1 / 2 * (klp + kpl + rl + rp) + 1 / 2 * sqrt((klp + kpl + rl + rp)^2 - 4 * (kpl * rl + klp * rp + rp * rl));
    k_minus = -1 / 2 * (klp + kpl + rl + rp) - 1 / 2 * sqrt((klp + kpl + rl + rp)^2 - 4 * (kpl * rl + klp * rp + rp * rl));

    L = k * kpl / (k_plus - k_minus) * conv(Pv, exp(k_plus * t) - exp(k_minus * t), 'full');
    Pt = k / (k_plus - k_minus) * conv(Pv, (k_plus + klp + rl) * exp(k_plus * t) - (k_minus + klp + rl) * exp(k_minus * t), 'full');

    y(1, :) = L(1:length(t));  % Lactate signal
    y(2, :) = Pt(1:length(t)) + Pv(1:length(t)); % Pyruvate signal
end
