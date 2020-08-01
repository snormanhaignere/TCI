function [diff_context, same_context, same_context_err, same_context_twogroups] = ...
    cross_context_corr_twogroups_helper(...
    Y_seg, Y_embed_seg, ...
    samedur_order_pairs, diffdur_order_pairs, ...
    samecontext_rep_pairs, diffcontext_rep_pairs, ...
    samedur_seg_pairs, diffdur_seg_pairs, ...
    make_samedur_comparisons, make_diffdur_comparisons, ...
    interleave_samedur, interleave_diffdur, ...
    simfunc, tranweightfn, trancorrfn)

n_lags = size(Y_seg{1}, 2);
same_context_twogroups = zeros(n_lags,2);
diff_context = zeros(n_lags,1);
same_context_weight_twogroups = [0, 0];
diff_context_weight = 0;

%% Compare same duration segments

if make_samedur_comparisons
    
    % weight by number of segments in the comparison
    weight = tranweightfn(size(samedur_seg_pairs,1));
    X = cell(1,2);
    for k = 1:size(samedur_order_pairs,2)
        p_orders = samedur_order_pairs(:,k);
        % assert(p_orders(1)~=p_orders(2));
        
        % select segs to be compared
        X{1} = Y_seg{1}(samedur_seg_pairs(:,1), :, p_orders(1), :);
        X{2} = Y_seg{2}(samedur_seg_pairs(:,2), :, p_orders(2), :);
        assert(all(size(X{1})==size(X{2})));
        
        % optionally interleave values
        if interleave_samedur
            [X{1}, X{2}] = interleave_oddeven(X{1}, X{2});
        end
        
        % correlate pairs of repetitions
        for l = 1:size(diffcontext_rep_pairs,2)
            p_reps = diffcontext_rep_pairs(:,l);
            C = trancorrfn(simfunc(X{1}(:,:,1,p_reps(1)), X{2}(:,:,1,p_reps(2))));
            if all(~isnan(C))
                try
                    diff_context = diff_context + C' * weight;
                    diff_context_weight = diff_context_weight + weight;
                catch
                    keyboard
                end
            end
        end
        
        % reliability of each element
        for l = 1:size(samecontext_rep_pairs,2)
            p_reps = samecontext_rep_pairs(:,l);
            for m = 1:2
                C = trancorrfn(simfunc(X{m}(:,:,1,p_reps(1)), X{m}(:,:,1,p_reps(2))));
                if all(~isnan(C))
                    same_context_twogroups(:,m) = same_context_twogroups(:,m) + C' * weight;
                    same_context_weight_twogroups(m) = same_context_weight_twogroups(m) + weight;
                end
            end
        end
    end
    clear X;
end

%% Different duration

if make_diffdur_comparisons
    n_longer_seg_durs = length(diffdur_seg_pairs);
    X = cell(1,2);
    for j = 1:n_longer_seg_durs
        weight = tranweightfn(size(diffdur_seg_pairs{j},1)); % standard error
        for k = 1:size(diffdur_order_pairs,2)
            p_orders = diffdur_order_pairs(:,k);
            
            % select segs to be compared
            X{1} = Y_seg{1}(diffdur_seg_pairs{j}(:,1), :, p_orders(1), :);
            X{2} = Y_embed_seg{2}(diffdur_seg_pairs{j}(:,1), :, p_orders(2), :, j);
            if ~(all(size(X{1})==size(X{2})))
                keyboard;
            end
            
            % optionally interleave values
            if interleave_diffdur
                [X{1}, X{2}] = interleave_oddeven(X{1}, X{2});
            end
            
            % correlate pairs of repetitions
            for l = 1:size(diffcontext_rep_pairs,2)
                p_reps = diffcontext_rep_pairs(:,l);
                C = trancorrfn(simfunc(X{1}(:,:,1,p_reps(1)), X{2}(:,:,1,p_reps(2))));
                if all(~isnan(C))
                    diff_context = diff_context + C' * weight;
                    diff_context_weight = diff_context_weight + weight;
                end
            end
            
            % reliability of each element
            for l = 1:size(samecontext_rep_pairs,2)
                p_reps = samecontext_rep_pairs(:,l);
                for m = 1:2
                    C = trancorrfn(simfunc(X{m}(:,:,1,p_reps(1)), X{m}(:,:,1,p_reps(2))));
                    if all(~isnan(C))
                        same_context_twogroups(:,m) = same_context_twogroups(:,m) + C' * weight;
                        same_context_weight_twogroups(m) = same_context_weight_twogroups(m) + weight;
                    end
                end
            end
        end
    end
    clear X;
end

%% Divide by weights, combine groups

% divide by weights
same_context_twogroups = bsxfun(@times, same_context_twogroups, 1./same_context_weight_twogroups);
diff_context = diff_context/diff_context_weight;

% average and substract the two groups
same_context = same_context_twogroups(:,1)/2 + same_context_twogroups(:,2)/2;
same_context_err = (same_context_twogroups(:,1)/2 - same_context_twogroups(:,2)/2).^2;
