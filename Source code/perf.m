function temp_obj_eval = perf(u, u0_GT)

[xm ym] = size(u);

% TP pixels
temp_tp = [u == 0] & [u0_GT == 0];

% FP pixels
temp_fp = [u == 0] & [u0_GT ~= 0];

% FN pixels
temp_fn = [u ~= 0] & [u0_GT == 0];

% TN pixels 
temp_tn = [u ~= 0] & [u0_GT ~= 0];

count_tp = sum(sum(temp_tp));
count_fp = sum(sum(temp_fp));
count_fn = sum(sum(temp_fn));
count_tn = sum(sum(temp_tn));

% precision
temp_p = count_tp / (count_fp + count_tp);

% recall
temp_r = count_tp / (count_fn + count_tp);

% F-Measure
temp_F = 2 * ( ( temp_r * temp_p ) / ( temp_r + temp_p ) );

temp_obj_eval.precision = temp_p;
temp_obj_eval.recall = temp_r;
temp_obj_eval.Fmeasure = temp_F;
temp_obj_eval.Falsepositives = count_fp;

end