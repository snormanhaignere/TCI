function [M, MAT_file] = modelfit_lagerr_cross_segdur_modular(L, varargin)

% Fit model to lagged correlation data to quantitatively estimate 
% integration window
% 
% 2019-04-05: Last edited, Sam NH

global root_directory;
root_directory = my_root_directory;
addpath(genpath([root_directory '/export_fig_v3']));

% optional parameters
I.plot_win = L.lag_t([1 end]);
I.plot_figure = true;
I.keyboard = false;
I.overwrite = false;
I.distr = 'gamma';
I.forcecausal = false;
I.bestcausal = true; % select only amongst causal solutions (not relevant for Gausssian)
I.plotcausal = true; % only plot errors for causal solutions

% other plotting parameters
I.winrange = [0, 6]; % in octaves relative to base
I.winbase = 0.03125; % in seconds
I.winhop = 1;
I.winsize = 1;

% can optionally not run analysis and just return MAT file
I.run = true;

switch I.distr
    case 'gauss'
        I.shape = 1;
    case 'gamma'
        I.shape = logspace(log10(1), log10(10), 10);
    otherwise
        error('No matching distribution');
end
[I, ~, C_value] = parse_optInputs_keyvalue(varargin, I);

if I.keyboard
%     keyboard;
end

% string identifying parameters of model fit
param_string_modelfit = [...
    struct2string(I, 'include_fields', {'distr'}), ...
    '_' struct2string(C_value, 'include_fields', {'shape','forcecausal','bestcausal'})];
if param_string_modelfit(end)=='_'; param_string_modelfit = param_string_modelfit(1:end-1); end

% file to save results to
MAT_file = [L.output_directory '/model_fit_' param_string_modelfit '-v2-fixed-NaN-handling.mat'];

% can stop and just return the MAT_file
% (useful to check if it exists)
if ~I.run  || isempty(L.channels)
    M = [];
    return;
end


%% Model fit

if ~exist(MAT_file, 'file') || I.overwrite
    
    [n_lags, n_seg, n_channels] = size(L.same_context);
    assert(n_channels == length(L.channels));
    assert(n_lags == length(L.lag_t))
    
    % measure to predict with the model
    % M.Y: lag x stim x channel
    same_context = L.same_context;
    diff_context = L.diff_context;
    M.Y = nan(size(same_context));
    for k = 1:n_seg
        M.Y(:,k,:) = diff_context(:,k,:)-same_context(:,k,:);
    end
    
    % move seg dimension to the end
    % -> lag x channel x seg
    Yshift = permute(M.Y, [1, 3, 2]);
    
    M.intper_sec = logspace(log10(min(L.segdurs)), log10(max(L.segdurs)), 20);
    M.delay_smp = 1:0.5*L.sr;
    M.shape = I.shape;
    M.err = zeros(length(M.intper_sec), length(M.delay_smp), length(M.shape), n_channels);
    M.causal_win = false(length(M.intper_sec), length(M.delay_smp), length(M.shape));
    for m = 1:length(M.shape)
        
        for i = 1:length(M.intper_sec)
            
            fprintf('shape: %.2f, intper: %.2f ms\n', M.shape(m), M.intper_sec(i)*1000);
            drawnow;
            for j = 1:length(M.delay_smp)
                
                                
                %                 % parameters to run code within just this loop
                %                 I.distr = 'gamma';
                %                 L.lag_t = [0 1];
                %                 L.sr = 100;
                %                 i = 1;
                %                 j = 1;
                %                 m = 1;
                %                 M.intper_sec = 0.1;
                %                 M.delay_smp = L.sr*0;
                %                 M.shape = 1;
                %                 i = 10;
                %                 j = 1;
                %                 m = 1;

                %                 % time vector in samples
                %                 t = -2*L.lag_t(end)*L.sr:2*L.lag_t(end)*L.sr;
                %
                %                 switch I.distr
                %                     case 'gauss'
                %
                %                         % gaussian window
                %                         % sigma when int period = central 95%
                %                         sig_sec = M.intper_sec(i) / 3.92;
                %                         sig_smp = L.sr * sig_sec;
                %                         h = normpdf(t - M.delay_smp(j), 0, sig_smp);
                %                         if I.forcecausal
                %                             h(t < 0) = 0;
                %                         end
                %                         h = h / sum(h);
                %
                %                         % plot
                %                         figure;
                %                         plot(t/L.sr, h);
                %                         xlim(3*M.intper_sec(i)*[-1, 1] + M.delay_smp(j)/L.sr);
                %
                %                         % Gaussian window is never causal
                %                         M.causal_win(i,j,m) = false;
                %
                %                     case 'gamma'
                %
                %                         % gaussian window
                %                         % sigma corresponding to 95%
                %                         a = M.shape(m);
                %                         b = 1/M.shape(m);
                %
                %                         % time vector in samples and seconds
                %                         t = -2*L.lag_t(end)*L.sr:2*L.lag_t(end)*L.sr;
                %                         t_sec = t/L.sr;
                %
                %                         % ratio which to scale stimulus
                %                         default_intper = gaminv(0.975,a,b) - gaminv(0.025,a,b);
                %                         r = M.intper_sec(i)/default_intper;
                %
                %                         % offset to adust delay
                %                         min_peak = max((a-1)*b,0)*r;
                %                         c = M.delay_smp(j)/L.sr - min_peak;
                %                         if M.delay_smp(j)/L.sr < min_peak
                %                             M.causal_win(i,j,m) = false;
                %                         else
                %                             M.causal_win(i,j,m) = true;
                %                         end
                %
                %                         % gamma distribution
                %                         h = gampdf((t_sec-c)/r, a, b);
                %                         if I.forcecausal
                %                             h(t_sec < 0) = 0;
                %                         end
                %                         h = h / sum(h);
                %
                %                         % plot
                %                         % figure;
                %                         % plot(t_sec, h);
                %                         % xlim(3*M.intper_sec(i)*[-1, 1] + M.delay_smp/L.sr)
                %
                %                     otherwise
                %                         error('No matching distribution');
                %                 end
                
                % calculate predictions for each segment
                
                used_predictors = [];
                used_segdurs = [];
                used_ipreds = [];
                for k = 1:n_seg
                    
                    %                     % rectangular segment
                    %                     seg = zeros(size(t));
                    %                     seg_dur_smp = round(L.sr * L.segdurs(k));
                    %                     seg(t <= seg_dur_smp & t >= 0) = 1;
                    %                     % plot(t, [h'/max(h(:)), seg'/max(seg(:))]);
                    %
                    %                     % convolve receptive field with segment
                    %                     overlap = myconv(seg', h', 'causal', false);
                    %                     % plot(t, area);
                    %                     % ylim([0 1]);
                    
                    predictor_found = false;
                    if ~isempty(used_segdurs)
                        xi = abs(L.segdurs(k)-used_segdurs)<1e-3;
                        if any(xi)
                            assert(sum(xi)==1);
                            predictor = used_predictors(:,xi);
                            ipred = used_ipreds(xi, :);
                            predictor_found = true;
                        end
                    end
                    
                    if ~predictor_found
                        [overlap,~,t_sec,M.causal(i,j,m)] = winoverlap(L.segdurs(k), ...
                            I.distr, M.intper_sec(i), M.delay_smp(j)/L.sr, 'shape', M.shape(m), ...
                            'win', [-L.lag_t(end), L.lag_t(end)], 'sr', L.sr, 'forcecausal', I.forcecausal);
                        % prdictions for each electrode
                        predictor = overlap(t_sec >= 0);
                        predictor = predictor(1:n_lags);
                        predictor = predictor(:);
                        predictor = 1-predictor;
                        ipred = pinv(predictor);
                        used_predictors = [used_predictors, predictor]; %#ok<AGROW>
                        used_ipreds = [used_ipreds; ipred]; %#ok<AGROW>
                        used_segdurs = [used_segdurs; L.segdurs(k)]; %#ok<AGROW>
                    end
                    
                    if any(any(isnan(Yshift(:,:,k))))
                        for q = 1:n_channels
                            xi = ~isnan(Yshift(:,q,k));
                            Yh = predictor(xi)*(pinv(predictor(xi))*Yshift(xi,q,k));
                            M.err(i,j,m,q) = M.err(i,j,m,q) + mean((Yh-Yshift(xi,q,k)).^2,1);
                        end
                    else
                        Yh = predictor*(ipred*Yshift(:,:,k));
                        M.err(i,j,m,:) = M.err(i,j,m,:) + reshape(mean((Yh-Yshift(:,:,k)).^2,1),[1,1,1,n_channels]);
                    end
                end
            end
        end
    end
    clear Yh Yshift xi;
    
    % find best prediction
    M.best_intper_sec = nan(n_channels,1);
    M.best_delay_smp = nan(n_channels,1);
    M.best_shape_smp = nan(n_channels,1);
    for q = 1:n_channels
        X = M.err(:,:,:,q);
        if any(M.causal_win(:)) && I.bestcausal
            X(~M.causal_win) = inf;
        end
        [~,xi] = min(X(:));
        [a,b,c] = ind2sub(size(X), xi);
        M.best_intper_sec(q) = M.intper_sec(a);
        M.best_delay_smp(q) = M.delay_smp(b);
        M.best_shape(q) = M.shape(c);
        clear X;
    end
    
    M.channels = L.channels;
    
    save(MAT_file, 'M', '-v7.3');
    
else
    
    load(MAT_file, 'M')
    
end

%% Plot the results

if I.plot_figure
    
    segdur_centers = I.winbase * 2.^(I.winrange(1):I.winhop:I.winrange(2));
    n_win = length(segdur_centers);
    
    figh = figure;
    
    plot_win_string = [num2str(I.plot_win(1)) '-' num2str(I.plot_win(2))];
    
    n_channels = length(L.channels);
    for q = 1:n_channels
        chan = L.channels(q);
        
        if strcmp(L.chnames{chan}, ['ch' num2str(chan)])
            chname = L.chnames{chan};
        else
            chname = ['ch' num2str(chan) '-' L.chnames{chan}];
        end
        
        % difference
        n_lags = length(L.lag_t);
        Z = nan(n_lags, n_win);
        Zp = zeros(n_lags, n_win);
        for i = 1:n_win
            segdur_range = segdur_centers(i) * 2.^(I.winsize*[-1 1]/2);
            segs_within_range = find(L.segdurs>segdur_range(1) & L.segdurs<segdur_range(2));
            Z(:,i) = nanmean(M.Y(:,segs_within_range,q),2);
            
            Yh = nan(n_lags, length(segs_within_range));
            for j = 1:length(segs_within_range)           
                [overlap,~,t_sec] = winoverlap(L.segdurs(segs_within_range(j)), ...
                    I.distr, M.best_intper_sec(q), M.best_delay_smp(q)/L.sr, 'shape', M.best_shape(q), ...
                    'win', [-L.lag_t(end), L.lag_t(end)], 'sr', L.sr, 'forcecausal', I.forcecausal, 'plot', false);
                predictor = overlap(t_sec >= 0);
                predictor = predictor(1:n_lags);
                predictor = predictor(:);
                predictor = 1-predictor;
                y = M.Y(:,segs_within_range(j),q);
                xi = ~isnan(y);
                Yh(xi,j) = predictor(xi)*pinv(predictor(xi))*y(xi);
            end
            Zp(:,i) = nanmean(Yh);
        end
        
        clf(figh);
        set(figh, 'Position', [100, 100, 900, 600]);
        n_rows = ceil(sqrt(n_win));
        n_cols = ceil(n_win/n_rows);
        ybounds = quantile(Z(:), [0, 0.99]);
        ybounds(1) = min(ybounds(1), 0);
        for i = 1:n_win
            subplot(n_rows, n_cols, i);
            plot(L.lag_t*1000, zeros(size(L.lag_t)), 'k-', 'LineWidth', 2);
            hold on;
            plot(segdur_centers(i)*[1 1]*1000, ybounds, 'r--', 'LineWidth', 2);
            plot(L.lag_t*1000, Z(:,i), 'LineWidth', 2);
            plot(L.lag_t*1000, Zp(:,i), 'LineWidth', 2);
            ylim(ybounds);
            xlim(I.plot_win*1000);
            title(sprintf('%.0f ms', segdur_centers(i)*1000));
            xlabel('Time (ms)');
            ylabel('Error');
        end
        fname = [L.figure_directory '/' chname '-model-prediction-lineplot_' param_string_modelfit '_plotwin' plot_win_string];
        print_wrapper([fname '.png']);
        print_wrapper([fname '.pdf']);
        
        % plot the error vs. parameters
        X = M.err(:,:,M.best_shape(q)==M.shape,q);
        if any(M.causal_win(:)) && I.plotcausal
            causal_string = '_only-causal';
            X(~M.causal_win(:,:,M.best_shape(q)==M.shape)) = NaN;
        else
            causal_string = '';
        end
        bounds = quantile(-X(:), [0.8, 1]);
        if ~all(isnan(X(:)))
            clf(figh);
            set(figh, 'Position', [100, 100, 600, 600]);
            
            imagesc(-X, bounds);
            colormap(parula(128));
            colorbar;
            ytick = interp1(log2(M.intper_sec), 1:length(M.intper_sec), log2(M.intper_sec));
            set(gca, 'YTick', ytick, 'YTickLabel', M.intper_sec*1000);
            xticks = get(gca, 'XTick');
            set(gca, 'XTickLabel', 1000*M.delay_smp(xticks)/L.sr);
            xlabel('Delay (ms)'); ylabel('Receptive Field Duration (ms)');
            title(sprintf('rf=%.f ms, delay=%.f ms', M.best_intper_sec(q)*1000, M.best_delay_smp(q)/L.sr*1000));
            set(gca, 'FontSize', 12);
            fname = [L.figure_directory '/' chname '-model-error_' param_string_modelfit causal_string];
            print_wrapper([fname '.png']);
            % export_fig([fname '.png'], '-png', '-transparent', '-r100');
            clear X;
        end
    end
    
end