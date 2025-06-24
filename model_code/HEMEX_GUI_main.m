% @author: Nichlas Vous Christensen
% @email: nvc@clin.au.dk
% @phone: +45 23464522
% @organization: Aarhus University, The MR Research Centre
% June 2025

function HEMEX_GUI_main(data, flip_P, flip_L, TR, bounds)
    % data: 5D dataset in the format [x, y, z, time, metabolite] with metabolite = 1 being pyruvate
    % flip_P: flip angle on pyruvate in degrees
    % flip_L: flip angle on lactate in degrees
    % TR: repetition time in seconds
    % bounds: parameter bounds used in fitting (optional)
    
    % Optional parameters
    if ~exist('bounds','var')
        % lower bound | start guess | upper bound
        bounds = [0 0.01 1;... % kpl (pyrurvate-to-lactate conversion rate)
                  0 0 0;... % klp (lactate-to-pyruvate conversion rate)
                  0.01 1/30 0.05;... % rp (1/T1 pyruvate relaxation)
                  0.8 1 1.2;... % rl (1/T1 lactate relaxation) [scaled from rp]
                  0.05 1 20;... % k (pyruvate permability) [scaled from kpl]
                  0 0 10;... % t0_delay (time-delay between AIF and voxel data)
                  0.5 5 30;... % mu (input parameter to the hemodynamic residue function)
                  0.5 3 30]; % sigma (input parameter to the hemodynamic residue function)
    end
    params_lb = bounds(:,1);
    params0 = bounds(:,2);
    params_ub = bounds(:,3);

    % Normalize data based on flip-angle
    data(:,:,:,:,1) = data(:,:,:,:,1)/sind(flip_P);
    data(:,:,:,:,2) = data(:,:,:,:,2)/sind(flip_L);
    
    % Normalize data to 1
    data = data/max(data(:));

    % Pulse contribution to relaxation rate, see Hill et al. 2013 (Model Free Approach)
    r2p = -log(cosd(flip_P))/TR;
    r2l = -log(cosd(flip_L))/TR;

    % Create time-array based on TR and data length
    timepoints = size(data,4);
    t = linspace(0,timepoints-1,timepoints)*TR;
    
    % Initialize the GUI figure
    fig = figure('Name', 'HEMEX GUI', 'NumberTitle', 'off', ...
                 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

    % Set up initial parameters
    slice_num = 1;
    time_point = 1;
    metabolites = {'Pyruvate', 'Lactate'};
    metabolite_index = 1;
    colormaps = ["Default", "Hot"];
    colors = [1 0 0; 0 0 1]; % Red, Blue
    colormap_index = 1;
    max_intensity = 0.2;
    mu_scale = 1;
    kPL_scale = 1;
    k_scale = 1;
    R2_threshold = 0.7;
    
    % Create axes for plots
    main_axes = axes('Parent', fig, 'Position', [0.1, 0.62, 0.8, 0.3]);
    aif_image_axes = axes('Parent', fig, 'Position', [0.1, 0.38, 0.17, 0.22]);
    aif_plot_axes = axes('Parent', fig, 'Position', [0.3, 0.38, 0.17, 0.22]);
    mask_image_axes = axes('Parent', fig, 'Position', [0.53, 0.38, 0.17, 0.22]);
    mask_plot_axes = axes('Parent', fig, 'Position', [0.73, 0.38, 0.17, 0.22]);
    MTT_image_axes = axes('Parent', fig, 'Position', [0.1, 0.03, 0.17, 0.2]);title('MTT map');
    kPL_image_axes = axes('Parent', fig, 'Position', [0.3, 0.03, 0.17, 0.2]);title('k_P_L map');
    k_image_axes = axes('Parent', fig, 'Position', [0.5, 0.03, 0.17, 0.2]);title('k map');

    % Create UI controls
    time_point_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.11, 0.93, 0.15, 0.03], ...
                                  'Min', 1, 'Max', size(data, 4), 'Value', time_point, ...
                                  'SliderStep', [1/(size(data,4)-1), 1/(size(data,4)-1)], ...
                                  'Callback', @update_plot);
    time_point_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.11, 0.96, 0.15, 0.03], ...
                                 'String', ['Time-point: ',num2str(time_point)]);

    if size(data,3) > 1 % If there are more than one slice, create slice-slider
        slice_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.33, 0.93, 0.15, 0.03], ...
                                 'Min', 1, 'Max', size(data, 3), 'Value', slice_num, ...
                                 'SliderStep', [1/(size(data,3)-1), 1/(size(data,3)-1)], ...
                                 'Callback', @update_plot);
        slice_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.33, 0.96, 0.15, 0.03], ...
                            'String', ['Slice no.: ',num2str(slice_num)]);
    end

    uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.51, 0.947, 0.06, 0.03], 'String', 'Metabolite');
    metabolite_menu = uicontrol('Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.57, 0.95, 0.08, 0.03], ...
                                'String', metabolites, 'Callback', @update_plot);
    
    uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.51, 0.917, 0.06, 0.03], 'String', 'Colormap');
    colormap_menu = uicontrol('Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.57, 0.92, 0.08, 0.03], ...
                                'String', colormaps, 'Callback', @update_plot);

    scale_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.73, 0.93, 0.15, 0.03], ...
                             'Min', 0.001, 'Max', 1, 'Value', max_intensity, ...
                             'SliderStep', [0.01, 0.1], 'Callback', @update_plot);
    scale_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.73, 0.96, 0.15, 0.03], ...
                            'String', ['Image scale: ',num2str(max_intensity)]);
                        
    R2_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.68, 0.22, 0.14, 0.03], ...
                             'Min', -10, 'Max', 0.99, 'Value', R2_threshold, ...
                             'SliderStep', [0.01, 0.1], 'Callback', @update_plot);
    R2_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.68, 0.25, 0.14, 0.03], ...
                            'String', ['R2 threshold: ',num2str(R2_threshold)]);
                        
    scale_mu_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.68, 0.15, 0.14, 0.03], ...
                             'Min', 0.001, 'Max', 1 , 'Value', mu_scale, ...
                             'SliderStep', [0.01, 0.1], 'Callback', @update_plot);
    scale_mu_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.68, 0.18, 0.14, 0.03], ...
                            'String',['MTT scale: ',num2str(mu_scale)]);
                        
    scale_kPL_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.68, 0.08, 0.14, 0.03], ...
                             'Min', 0.001, 'Max', 1 , 'Value', kPL_scale, ...
                             'SliderStep', [0.01, 0.1], 'Callback', @update_plot);
    scale_kPL_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.68, 0.11, 0.14, 0.03], ...
                            'String',['kPL scale: ',num2str(kPL_scale)]);
                        
    scale_k_slider = uicontrol('Style', 'slider', 'Units', 'normalized', 'Position', [0.68, 0.01, 0.14, 0.03], ...
                             'Min', 0.001, 'Max', 1 , 'Value', k_scale, ...
                             'SliderStep', [0.01, 0.1], 'Callback', @update_plot);
    scale_k_label = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [0.68, 0.04, 0.14, 0.03], ...
                            'String',['k scale: ',num2str(k_scale)]);
                        
    % Button interactions
    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.17, 0.3, 0.07, 0.03], ...
              'String', 'New AIF', 'Callback', @new_aif);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.25, 0.3, 0.07, 0.03], ...
              'String', 'Save AIF', 'Callback', @save_aif);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.33, 0.3, 0.07, 0.03], ...
              'String', 'Load AIF', 'Callback', @use_saved_aif);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.51, 0.3, 0.07, 0.03], ...
              'String', 'New Mask', 'Callback', @new_mask);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.59, 0.3, 0.07, 0.03], ...
              'String', 'Save Mask', 'Callback', @save_mask);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.67, 0.3, 0.07, 0.03], ...
              'String', 'Load Mask', 'Callback', @use_saved_mask);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.75, 0.3, 0.07, 0.03], ...
              'String', 'Start Fitting', 'Callback', @start_fitting);
          
    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.83, 0.3, 0.07, 0.03], ...
              'String', 'View Fit', 'Callback', @view_fit);
          
    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.83, 0.26, 0.07, 0.03], ...
              'String', 'Explore Fit', 'Callback', @explore_fit);

    uicontrol('Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.83, 0.22, 0.07, 0.03], ...
              'String', 'Save Fit', 'Callback', @save_fit);

    % Initialize variables for saving state
    AIF_voxel = [];
    AIF_fit = [];
    AIF_data = [];
    mask = ones(size(data,1), size(data,2));
    model_fit_results = [];
    model_fit = [];
    model_voxel = [];
    model_data = [];

    % Plot the initial data
    update_plot();
    
    % Function to update plots based on user input
    function update_plot(~, ~)      
        % Get current values from UI controls
        time_point = round(get(time_point_slider, 'Value'));
        metabolite_index = get(metabolite_menu, 'Value');
        colormap_index = get(colormap_menu, 'Value');
        max_intensity = get(scale_slider, 'Value');
        mu_scale = get(scale_mu_slider, 'Value');
        kPL_scale = get(scale_kPL_slider, 'Value');
        k_scale = get(scale_k_slider, 'Value');
        R2_threshold = get(R2_slider, 'Value');
        if size(data,3) > 1 % If there are more than one slice
            slice_num = round(get(slice_slider, 'Value'));
        end

        % Update labels        
        set(time_point_label, 'String', ['Time-point: ',num2str(time_point)]);
        set(scale_label, 'String',['Image scale: ',num2str(max_intensity)]);
        set(scale_mu_label, 'String', ['MTT scale: ',num2str(mu_scale)]);
        set(scale_kPL_label, 'String', ['kPL scale: ',num2str(kPL_scale)]);
        set(scale_k_label, 'String', ['k scale: ',num2str(k_scale)]);
        set(R2_label, 'String', ['R2 threshold: ',num2str(R2_threshold)]);
        if size(data,3) > 1 % If there are more than one slice
            set(slice_label, 'String', ['Slice no.: ',num2str(slice_num)]);
        end

        % Display selected slice and time-point in main_axes
        axes(main_axes);
        cla(main_axes); % Clear axes before plotting new data
        images = reshape(squeeze(data(:, :, slice_num, :, :)), [size(data,1), size(data,2), size(data,4) * size(data,5)]);
        montage(images, 'Size', [size(data,5) size(data,4)], 'ThumbnailSize', [size(data,1), size(data,2)], 'DisplayRange', [0 max_intensity]), colormap(colormaps(colormap_index));
        
        % Mark the current shown slice and time-point in the AIF and mask image
        rectangle('Position',[size(data,1)*time_point-size(data,1) size(data,2)*metabolite_index-size(data,2) size(data,1) size(data,2)],'EdgeColor',colors(colormap_index,:))

        % Add text labels next to rows
        text(-size(data,1)*0.25, size(data,2)/2, 'Pyr', 'Parent', main_axes, 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center');
        text(-size(data,1)*0.25, size(data,2)/2 + size(data,1), 'Lac', 'Parent', main_axes, 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center');

        % Plot AIF image
        axes(aif_image_axes);
        cla(aif_image_axes); % Clear axes before plotting new data
        imagesc(data(:, :, slice_num, time_point, metabolite_index), [0 max_intensity]);
        title(aif_image_axes, 'AIF Image');
        xlabel('x')
        ylabel('y')
        if ~isempty(AIF_voxel) & slice_num == AIF_voxel(3)
            hold(aif_image_axes, 'on');
            plot(AIF_voxel(1), AIF_voxel(2), 'o','Color',colors(colormap_index,:), 'MarkerSize', 10, 'LineWidth', 2);
            hold(aif_image_axes, 'off');
        end
        
        % Plot current AIF
        if ~isempty(AIF_fit)
            axes(aif_plot_axes);
            cla(aif_plot_axes);
            title(aif_plot_axes, ['AIF Curve (x = ',num2str(AIF_voxel(1)),', y = ',num2str(AIF_voxel(2)),', z = ',num2str(AIF_voxel(3)),')']);
            hold(aif_plot_axes, 'on');
            plot(t, AIF_data, 'ko','MarkerSize', 6);
            plot(t, AIF_fit, 'k-','linewidth',1.5);
            ylim([0 max(AIF_data)*1.5])
            xlabel('Time (s)')
            legend('AIF','Fit')
            hold(aif_plot_axes, 'off');
        end
        
        % Plot mask image
        axes(mask_image_axes);
        cla(mask_image_axes); % Clear axes before plotting new data
        data_masked = data.*mask; % Use masked data
        imagesc(data_masked(:, :, slice_num, time_point, metabolite_index), [0 max_intensity]);
        title(mask_image_axes, 'Mask for fit');
        xlabel('x')
        ylabel('y')
        if ~isempty(model_voxel) & slice_num == model_voxel(3)
            hold(mask_image_axes, 'on');
            plot(mask_image_axes, model_voxel(1), model_voxel(2), 'o','Color',colors(colormap_index,:), 'MarkerSize', 10, 'LineWidth', 2);
            hold(mask_image_axes, 'off');
        end
        
        % Plot current FIT
        if ~isempty(model_voxel)
            axes(mask_plot_axes);
            cla(mask_plot_axes);
            title(mask_plot_axes, ['Fit (x = ',num2str(model_voxel(1)),', y = ',num2str(model_voxel(2)),', z = ',num2str(model_voxel(3)),')']);
            hold(mask_plot_axes, 'on');
            plot(t, model_data(:,2), 'ro','MarkerSize', 6);
            plot(t, model_fit(2, :), 'r','linewidth',1.5);
            plot(t, model_data(:,1), 'bo','MarkerSize', 6);
            plot(t, model_fit(1, :), 'b','linewidth',1.5);
            xlabel('Time (s)')
            legend('Pyruvate','Pyruvate fit','Lactate','Lactate fit')           
        end
        
        if ~isempty(model_fit_results)
            R2_mask = model_fit_results.R2 < R2_threshold;
            % Plot MTT map
            axes(MTT_image_axes);
            cla(MTT_image_axes); % Clear axes before plotting new data
            R2_thresholded_mu = model_fit_results.mu_fits;
            mu_max = max(R2_thresholded_mu(:));
            R2_thresholded_mu(R2_mask) = 0;
            imagesc(R2_thresholded_mu(:, :, slice_num),[0 mu_scale*mu_max]); colorbar;
            title(MTT_image_axes, 'MTT map');
            if ~isempty(model_voxel) & slice_num == model_voxel(3)
                hold(MTT_image_axes, 'on');
                plot(MTT_image_axes, model_voxel(1), model_voxel(2), 'o','Color',colors(colormap_index,:), 'MarkerSize', 10, 'LineWidth', 2);
                hold(MTT_image_axes, 'off');
            end
            % Plot kPL map
            axes(kPL_image_axes);
            cla(kPL_image_axes); % Clear axes before plotting new data
            R2_thresholded_kPL = model_fit_results.kpl_fits;
            kPL_max = max(R2_thresholded_kPL(:));
            R2_thresholded_kPL(R2_mask) = 0;
            imagesc(R2_thresholded_kPL(:, :, slice_num),[0 kPL_scale*kPL_max]); colorbar;
            title(kPL_image_axes, 'kPL map');
            if ~isempty(model_voxel) & slice_num == model_voxel(3)
                hold(kPL_image_axes, 'on');
                plot(kPL_image_axes, model_voxel(1), model_voxel(2), 'o','Color',colors(colormap_index,:), 'MarkerSize', 10, 'LineWidth', 2);
                hold(kPL_image_axes, 'off');
            end
            % Plot k map
            axes(k_image_axes);
            cla(k_image_axes); % Clear axes before plotting new data
            % First redo the reparameterization of kpl-k dependency
            R2_thresholded_k = model_fit_results.k_fits;
            k_max = max(R2_thresholded_k(:));
            R2_thresholded_k(R2_mask) = 0;
            imagesc(R2_thresholded_k(:, :, slice_num),[0 k_scale*k_max+eps]); colorbar;
            title(k_image_axes, 'k map');
            if ~isempty(model_voxel) & slice_num == model_voxel(3)
                hold(k_image_axes, 'on');
                plot(k_image_axes, model_voxel(1), model_voxel(2), 'o','Color',colors(colormap_index,:), 'MarkerSize', 10, 'LineWidth', 2);
                hold(k_image_axes, 'off');
            end

        end
    end

    % Callback function to new AIF selection
    function new_aif(~, ~)
        axes(aif_image_axes);
        [x, y] = ginput(1);
        x = round(x); y = round(y); z = slice_num;
        AIF_voxel = [x, y, z];
        AIF_data = squeeze(data(y, x, slice_num, :, metabolite_index));
        AIF_fit = AIF_model_fit(t, AIF_data);
        update_plot();
    end

    % Callback function to save AIF
    function save_aif(~, ~)
        if ~isempty(AIF_fit)
            [file, path] = uiputfile('AIF.mat', 'Save AIF As');
            if isequal(file, 0) || isequal(path, 0)
                disp('User canceled the save dialog.');
            else
                fullFileName = fullfile(path, file);
                save(fullFileName, 'AIF_voxel', 'AIF_data', 'AIF_fit');
                disp(['AIF saved successfully to ', fullFileName]);
            end
        else
            disp('No AIF curve to save.');
        end
    end

    % Callback function to load AIF
    function use_saved_aif(~, ~)
        [file, path] = uigetfile('*.mat', 'Select AIF File');
        if isequal(file, 0) || isequal(path, 0)
            disp('User canceled file selection.');
        else
            fullFileName = fullfile(path, file);
            try
                load(fullFileName, 'AIF_voxel', 'AIF_data', 'AIF_fit');
                update_plot();
                disp(['Loaded AIF from ', fullFileName]);
            catch ME
                disp(['Failed to load AIF: ', ME.message]);
            end
        end
    end

    % Callback function to start mask selection
    function new_mask(~, ~)
        mask = ones(size(data,1),size(data,2)); % Reset mask before making a new
        update_plot();
        axes(mask_image_axes);
        h = impoly();
        mask = h.createMask();
        update_plot();
    end

    % Callback function to save mask
    function save_mask(~, ~)
        if ~isempty(mask)
            [file, path] = uiputfile('mask.mat', 'Save Mask As');
            if isequal(file, 0) || isequal(path, 0)
                disp('User canceled the save dialog.');
            else
                fullFileName = fullfile(path, file);
                save(fullFileName, 'mask');
                disp(['Mask saved successfully to ', fullFileName]);
            end
        else
            disp('No mask to save.');
        end
    end

    % Callback function to load mask
    function use_saved_mask(~, ~)
        [file, path] = uigetfile('*.mat', 'Select Mask File');
        if isequal(file, 0) || isequal(path, 0)
            disp('User canceled file selection.');
        else
            fullFileName = fullfile(path, file);
            try
                load(fullFileName, 'mask');
                update_plot();
                disp(['Loaded mask from ', fullFileName]);
            catch ME
                disp(['Failed to load mask: ', ME.message]);
            end
        end
    end

    % Callback function to start fitting
    function start_fitting(~, ~)
        disp('Initializing fitting...')
        model_fit_results = HEMEX_model_fit(data, AIF_fit, mask, t, r2p, r2l, params0, params_lb, params_ub);
        disp('Fitting finished.')
        update_plot();
    end

    % Callback function to viewing fit
    function view_fit(~, ~)
        axes(mask_image_axes);
        [x, y] = ginput(1);
        x = round(x); y = round(y); z = slice_num;
        model_voxel = [x, y, z];
        model_data = cat(2,squeeze(data(y,x,z,:,2)),squeeze(data(y,x,z,:,1)));
        if length(size(model_fit_results.par_fits)) > 3 % If there is a z-dimension
            model_fit = HEMEX_model(squeeze(model_fit_results.par_fits(y,x,z,:)), t, AIF_fit, r2p, r2l);
        else
            model_fit = HEMEX_model(squeeze(model_fit_results.par_fits(y,x,:)), t, AIF_fit, r2p, r2l);
        end 
        update_plot();
    end

    % Callback function to exploring fit
    function explore_fit(~, ~)
        x = model_voxel(1);
        y = model_voxel(2);
        z = model_voxel(3);
        pyr_data = squeeze(data(y,x,z,:,1));
        lac_data = squeeze(data(y,x,z,:,2));
        
        if length(size(model_fit_results.par_fits)) > 3 % If there is a z-dimension
            params = squeeze(model_fit_results.par_fits(y,x,z,:));
        else
            params = squeeze(model_fit_results.par_fits(y,x,:));
        end
        % Convert bounds to correct slider values for the co-dependent parameters (rl and k)
        slider_min = params_lb;
        slider_min(4) = slider_min(3)*slider_min(4);
        slider_min(5) = slider_min(1)*slider_min(5);
        slider_max = params_ub;
        slider_max(4) = slider_max(3)*slider_max(4);
        slider_max(5) = slider_max(1)*slider_max(5);
        % Reparameterize rl based on rp and k based on kpl
        params(4) = params(4) * params(3); % rl
        params(5) = params(5) * params(1); % k
        % Run the interactive plot
        HEMEX_interactive_model_plot(t, AIF_fit, pyr_data, lac_data, params, r2p, r2l, slider_min, slider_max)
    end

     % Callback function to save fit
    function save_fit(~, ~)
        if ~isempty(model_fit_results)
            [file, path] = uiputfile('model_fit_results.mat', 'Save Fit As');
            if isequal(file, 0) || isequal(path, 0)
                disp('User canceled the save dialog.');
            else
                fullFileName = fullfile(path, file);
                save(fullFileName, 'model_fit_results');
                disp(['Fit saved successfully to ', fullFileName]);
            end
        else
            disp('No fit data to save.');
        end
    end

end  

function model_fit_results = HEMEX_model_fit(data, AIF, mask, t, r2p, r2l, params0, params_lb, params_ub)    
    % Repeat mask for each slice
    mask = repmat(mask,[1,1,size(data,3)]);

    % Vectorize data
    pyr_norm = vectorize(data(:,:,:,:,1),mask);
    lac_norm = vectorize(data(:,:,:,:,2),mask);

    % Prepare variables for holding fitted curves
    pyr_fits = zeros(size(pyr_norm));
    lac_fits = zeros(size(lac_norm));
    
    nvox = size(pyr_norm,1);
    pars = zeros(nvox,length(params_lb));
    resnorm = zeros(nvox,1);
    residual = zeros(nvox,2,length(t));
    exitflag = zeros(nvox,1);
    output= cell(nvox,1);
    R2 = zeros(nvox,1);

    options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'Display','off','MaxFunctionEvaluations',10000,'MaxIterations',6000);
    
    disp(['Nvox = ', num2str(nvox)]);
    for j = 1:nvox
        if mod(nvox-j, 100) == 0
            disp(['Voxels left: ', num2str(nvox-j)]);
        end
        Pz = pyr_norm(j,:);
        Lz = lac_norm(j,:);
        AL = max(Lz);
        AP = max(Pz);
        weights = 1./[AL;AP];
        [pars(j,:),resnorm(j),residual(j,:,:),exitflag(j),output{j}]=lsqcurvefit(@(p,t)HEMEX_model(p, t, AIF, r2p, r2l).*weights,params0,t,[Lz;Pz].*weights,params_lb,params_ub,options);
               
        % First normalize the lac+lac_fit and pyr+pyr_fit together (to between [0,1])
        y_pred = HEMEX_model(pars(j,:), t, AIF, r2p, r2l);
        
        % Saved fitted curves
        pyr_fits(j,:) = y_pred(2,:);
        lac_fits(j,:) = y_pred(1,:);

        timepoints = length(Pz);
        norm_lac = normalize([Lz,y_pred(1,:)],'range')';
        norm_pyr = normalize([Pz,y_pred(2,:)],'range')';
        y_norm = [norm_lac(1:timepoints);norm_pyr(1:timepoints)];
        y_pred_norm = [norm_lac(timepoints+1:end);norm_pyr(timepoints+1:end)];
        % Calculate R2 for combined fit of both curves (equally weighted).
        SStot = sum((y_norm(:)-mean(y_norm(:))).^2); % Total Sum-Of-Squares
        SSres = sum((y_norm(:)-y_pred_norm(:)).^2); % Residual Sum-Of-Squares
        R2(j) = 1-SSres/SStot;         
    end
    
    %% Save results and perform reparameterization
    model_fit_results.par_fits = vectorize(pars,mask);
    model_fit_results.exitflag_fits =  vectorize(exitflag,mask);
    model_fit_results.kpl_fits = vectorize(pars(:,1),mask);
    model_fit_results.klp_fits = vectorize(pars(:,2),mask);
    model_fit_results.r1p_fits = vectorize(pars(:,3),mask); 
    model_fit_results.r1l_fits = vectorize(pars(:,4),mask) .* model_fit_results.r1p_fits;  % Reparameterize rl based on rp
    model_fit_results.k_fits = vectorize(pars(:,5),mask) .* model_fit_results.kpl_fits; % Reparameterize k based on kpl
    model_fit_results.t0_delay_fits = vectorize(pars(:,6),mask);
    model_fit_results.mu_fits = vectorize(pars(:,7),mask);
    model_fit_results.sigma_fits = vectorize(pars(:,8),mask);
    model_fit_results.R2 = vectorize(R2,mask);

    % Also save the curves and data for easier plotting
    model_fit_results.pyr_norm = vectorize(pyr_norm,mask);
    model_fit_results.lac_norm = vectorize(lac_norm,mask);
    model_fit_results.pyr_fits = vectorize(pyr_fits,mask);
    model_fit_results.lac_fits = vectorize(lac_fits,mask);
    model_fit_results.t = t;
end

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
    Pv = Pv(1:length(t));

    k_plus = -1/2*(klp+kpl+rl+rp) + 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));
    k_minus = -1/2*(klp+kpl+rl+rp) - 1/2 * sqrt((klp+kpl+rl+rp)^2-4*(kpl*rl+klp*rp+rp*rl));

    L = k*kpl/(k_plus-k_minus)*conv(Pv,exp(k_plus*t)-exp(k_minus*t),'full');
    Pt = k/(k_plus-k_minus)*conv(Pv,(k_plus+klp+rl)*exp(k_plus*t)-(k_minus+klp+rl)*exp(k_minus*t),'full');

    y(1,:) = L(1:length(t));  % Lactate signal
    y(2,:) = Pt(1:length(t))+Pv(1:length(t)); % Pyruvate signal
end

function [AIF_fit] = AIF_model_fit(t, AIF_raw)
    options = optimoptions('lsqcurvefit','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'Display','off');
    t0 = 2;
    a = 3;
    b = 1.5;
    A = max(AIF_raw);
    params0_AIF = [a, b, t0, A];
    params_fit_PA = lsqcurvefit(@(params, t) AIF_model(params, t), params0_AIF, t, AIF_raw', [], [], options);
    AIF_fit = AIF_model(params_fit_PA, t);
end

% Shifted Gamme variate
function [AIF] = AIF_model(params, t)
    a = params(1);
    b = params(2);
    t0 = params(3);
    A = params(4);

    AIF = A*(t-t0 > 0).*(t-t0).^(a-1).*exp(-(t-t0)/b);
end

function s = vectorize(S, mask)
    Ssz = size(S);
    Sdims = length(Ssz);
    masksz = size(mask);
    
    % Check if mask size matches all but the last dimension of S
    if Sdims > 2 && isequal(masksz, Ssz(1:ndims(mask)))
        % Permute S to bring the last dimension to the front
        s = permute(S, [Sdims, 1:Sdims-1]);
        % Apply the mask and transpose
        s = s(:, mask)';
        
    elseif ismatrix(S) && Ssz(1) == sum(mask(:))
        % If S is 2D and the size of its first dimension matches the sum of mask
        s = zeros([Ssz(2), masksz], "like", S);
        s(:, mask) = S';
        s = permute(s, [2:length(masksz) + 1, 1]);
    else
        error('Input dimensions are incompatible with mask.');
    end
end