% Computes matrix of confusions between models
%
% 2019-12-2: Created, Sam NH

clear I;
I.sr = 100;
I.distr = 'gamma';
I.intervalmass = 0.75;
I.intervaltype = 'highest';
I.shape = [1,1.5,2,3,5];
I.segdur_range = [1/32, 1/2];
I.intper_range = [1/32, 1/2];
I.nintper = 20;
I.delay_interval = 1/I.sr;
I.delay_range = [0, 1/4];
I.rampwin = 'hann';
I.rampdur = 1/32;
I.weight = true;
n_shapes = length(I.shape);

% integration period
clear M;
M.intper_sec = logspace(log10(I.intper_range(1)), log10(I.intper_range(2)), I.nintper);
n_intper_sec = length(M.intper_sec);
M.delay_sec_start = I.delay_range(1):I.delay_interval:I.delay_range(2);
M.delay_sec_start = round(M.delay_sec_start*I.sr)/I.sr;

lags = -2:1/I.sr:2;
n_lags = length(lags);

seg_durs = 2.^(log2(I.segdur_range(1)):1:log2(I.segdur_range(2)));
n_seg_durs = length(seg_durs);
weights = 1./seg_durs;
weights = weights / sum(weights);
Y_model = nan(n_lags, n_seg_durs, n_intper_sec, n_shapes);
for j = 1:n_shapes
    for i = 1:n_intper_sec
        fprintf('shape %.2f, intper %.0f ms\n', I.shape(j), M.intper_sec(i)*1000);drawnow;
        
        % calculate predictions for each segment
        for k = 1:n_seg_durs
            [winpow, ~, overlap] = win_power_ratio(seg_durs(k), ...
                I.distr, M.intper_sec(i), 0, ...
                'shape', I.shape(j), 'tsec', lags, ...
                'rampwin', I.rampwin, 'rampdur', I.rampdur, ...
                'intervalmass', I.intervalmass, 'delaypoint', 'start', ...
                'intervaltype', I.intervaltype);
            winpow(winpow<0)=0;
            Y_model(:,k,i,j) = (winpow);
        end
    end
end

Y_model_with_delays = add_delays(Y_model, checkint(M.delay_sec_start*I.sr));

%%

Ypow = squeeze_dims(sum(sum(bsxfun(@times, Y_model.^2, weights),1),2),[1,2]);
figure;
imagesc(Ypow)

%%


M.loss = nan(n_intper_sec, n_shapes, n_intper_sec, n_shapes);
M.best_delay = nan(n_intper_sec, n_shapes, n_intper_sec, n_shapes);
for j = 1:n_shapes
    for i = 1:n_intper_sec
        E = bsxfun(@minus, Y_model(:,:,i,j), Y_model_with_delays);
        mse = mean(mean(bsxfun(@times, abs(E).^2, weights),1),2);
        [min_mse, M.best_delay(:,:,i,j)] = min(mse,[],5);
        M.loss(:,:,i,j) = squeeze_dims(min_mse, [1,2]);
    end
end

%%
figure;
i = 15;
j = 3;
l = 5;
bounds = quantile(M.loss(:,j,i,j), [0 0.5]);
imagesc(flipud(M.loss(:,:,i,j)), bounds);
set(gca, 'YTick', 1:length(M.intper_sec), 'YTickLabel', fliplr(M.intper_sec)*1000);
% plot(log2(M.intper_sec), M.loss(:,10))

%%

figure;
n_rows = round(sqrt(n_seg_durs));
n_cols = ceil(n_seg_durs/n_rows);
h = nan(1,3);
for k = 1:n_seg_durs
    subplot(n_rows, n_cols, k)
    h(1) = plot(lags, Y_model(:,k,i,j));
    hold on;
    [~,min_intper] = min(M.loss(:,l,i,j));
    h(2) = plot(lags, Y_model_with_delays(:,k,i,l,M.best_delay(i,l,i,j)));
    h(3) = plot(lags, Y_model_with_delays(:,k,min_intper,l,M.best_delay(min_intper,l,i,j)));
    ylim([0,1]);
    legend(h, 'true', 'same', 'min');
end

% %%
% figure;
% i = 12;
% j = 4;
% plot(log2(M.intper_sec), M.loss(:,j,i,j))


