function ttmatgen2(fc1, fc2, f_vec, t_vec, ...
    duration_sweep, duration_silence, n_win_samps, overlap) % fam1, fam2,

% fc1 = sweep start freq (Hz)
% fc2 = sweep end freq (Hz)
% f_vec = frequency vector corresponding to the TFR (Hz)
% t_vec = time vector coresponding to the TFR (s)
% duration_sweep = length of the swepep (s)
% duration_silence = length of thew silent periods before and after the sweep (s)
% n_win_samps = umber of samples in the window function used in the TFR
% overlap = Amount of overlap used in the TFR (%)

rows_out = length(f_vec);
cols_out = length(t_vec);

% Convert freq and time vectors from numbers to strings for table indexing
f_labels = arrayfun(@num2str, f_vec, 'UniformOutput', 0);
t_labels = arrayfun(@num2str, t_vec, 'UniformOutput', 0);

% Create empty matrix with the desired dimensions
emptymat = zeros(rows_out, cols_out);

% Convert matrix to a table so we can index using non-integer freqs and times
table = array2table(emptymat, VariableNames=t_labels, RowNames=f_labels);

% t_res = stft_shortwin_tres*fs;
t_res = mean(diff(t_vec));
f_res = mean(diff(f_vec));
% alt method - might have to do this for cet because of log f_res
% f_res = round((fmax-fmin) / length(f_vec), 1);
% Could also try using the same method as for deriving n_cols_sweep.

% % Construct Vector of sweep frequency indices
% if fc1 > fc2 % case where the sweep decends in freq.
%     f_indices_sweep = fc1:-f_res:fc2;
% elseif fc1 < fc2 % case where the sweep ascends in freq.
%     f_indices_sweep = fc1:f_res:fc2;
% end

% Discrete sampling means sweep does not start at precisely the perscribed time.
% Find time index closest to prescribed start time
[~,sweep_tidx_start] = min(abs(t_vec - duration_silence));
sweep_tstart = t_vec(sweep_tidx_start);

% Discrete sampling means sweep does not end at precisely the perscribed time.
% Find time index closest to prescribed end time
[~,sweep_tidx_end] = min(abs(t_vec - ((duration_silence+duration_sweep)-t_res)));
sweep_tend = t_vec(sweep_tidx_end);

% Compressed frequency resolution means sweep may not start at precisely 
% the perscribed frequency. Find freq index closest to prescribed start
% frequency. 
[~,sweep_fidx_start] = min(abs(f_vec - fc1));
sweep_fstart = f_vec(sweep_fidx_start);

% Compressed frequency resolution means sweep may not end at precisely 
% the perscribed frequency. Find freq index closest to prescribed end
% frequency. 
[~,sweep_fidx_end] = min(abs(f_vec - fc2));
sweep_fend = f_vec(sweep_fidx_end);

% Construct vector of sweep time indices
t_indices_sweep = t_vec(sweep_tidx_start:sweep_tidx_end);
 
% Window overlap causes some extra time indices. Calculate number and
% remove:
extrasamps = floor((n_win_samps*(overlap/100)) / (t_res*fs));
t_indices_sweep = t_indices_sweep(1:end-extrasamps);

% Construct Vector of sweep frequency indices
if sweep_fstart > sweep_fend % case where the sweep decends in freq.
    f_indices_sweep = flip(f_vec(sweep_fidx_end:sweep_fidx_start));
elseif sweep_fstart < sweep_fend % case where the sweep ascends in freq.
    f_indices_sweep = f_vec(sweep_fidx_start:sweep_fidx_end);
end

% Convert sweep indices from numbers to strings for table indexing
f_sweep_labels = arrayfun(@num2str, f_indices_sweep, 'UniformOutput', 0);
t_sweep_labels = arrayfun(@num2str, t_indices_sweep, 'UniformOutput', 0);

% Write ones to the cells of the table at freq and time indices in
% f_indices_sweep and t_indices_sweep
for i=1:length(f_indices_sweep)
    f_index = f_sweep_labels(i);
    t_index = t_sweep_labels(i);
    table{f_index, t_index} = 1;
end

% % Construct amplitude modulation vector
% sig_am = rescale(sign(chirp(t_indices_sweep, fam1, t_indices_sweep(end), fam2, 'linear')));
% sig_am = repmat(sig_am, rows_out);

groundtruth = table2array(table);
imagesc(groundtruth)
end
