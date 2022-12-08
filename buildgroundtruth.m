function ground_truth = buildgroundtruth(fc1, fc2, fam1, fam2, f_vec, t_vec,...
    n_samps_sweep, n_freqs_total, sigma, silence, sig_am, fs)

% Calculate Start and End Frequencies for Upper Sideband Components
usb1_f1 = fc1 + fam1;
usb1_f2 = fc2 + fam2;
usb2_f1 = fc1 + (fam1 * 3);
usb2_f2 = fc2 + (fam2 * 3);
usb3_f1 = fc1 + (fam1 * 5);
usb3_f2 = fc2 + (fam2 * 5);

% Calculate Start and End Frequencies for Lower Sideband Components
lsb1_f1 = fc1 - fam1;
lsb1_f2 = fc2 - fam2;
lsb2_f1 = fc1 - (fam1 * 3);
lsb2_f2 = fc2 - (fam2 * 3);
lsb3_f1 = fc1 - (fam1 * 5);
lsb3_f2 = fc2 - (fam2 * 5);

% Error handling
sideband_freqs = [usb1_f1, usb1_f2, usb2_f1, usb2_f2, usb3_f1, usb3_f2,...
    lsb1_f1, lsb1_f2, lsb2_f1, lsb2_f2, lsb3_f1, lsb3_f2];

assert(all(sideband_freqs >= 0), ['ERROR: THE PRODUCT OF CARRIER AND ' ...
    'MODULATOR FREQUENCIES YOU HAVE ENTERED HAS RESULTED IN A FREQUENCY ' ...
    'COMPONENT WITH NEGATIVE FREQUENCY. THIS IS NOT SUPPORTED. PLEASE ' ...
    'INCREASE LOWEST CARRIER OR MODULATOR FREQUENCY'])

% Build amplitude vector
amp = ones(1, n_samps_sweep);

% Build matrices representing carrier and sideband components
carrier_MAT = tfmatgen2(f_vec, t_vec, ...
    duration_sweep, duration_silence, n_win_samps, overlap);

usb1_MAT = tfmatgen(usb1_f1, usb1_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);
usb2_MAT = tfmatgen(usb2_f1, usb2_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);
usb3_MAT = tfmatgen(usb3_f1, usb3_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);
lsb1_MAT = tfmatgen(lsb1_f1, lsb1_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);
lsb2_MAT = tfmatgen(lsb2_f1, lsb2_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);
lsb3_MAT = tfmatgen(lsb3_f1, lsb3_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence, cols_out, rows_out);

lsb1_MAT = lsb1_MAT(1:n_freqs_total, :);
lsb2_MAT = lsb2_MAT(1:n_freqs_total, :);
lsb3_MAT = lsb3_MAT(1:n_freqs_total, :);

% Combine matrices for all components
all_MAT = carrier_MAT + (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + (usb3_MAT .* 0.125) +...
    (lsb1_MAT .* 0.5) + (lsb2_MAT .* 0.25) +...
    (lsb3_MAT .* 0.125);

% Due to fmin being constant=0 in tfmatgen.m, and a variable in this
% function, the ground truth matrix is shifted up in frequency by the value
% of fmin/f_res+1. Use circshift to shift it back down.
all_MAT = circshift(flip(all_MAT, 1), -(fmin/f_res)-1, 1);

% That circshift moves some data that is below fmin to the high frequencies.
% It shouldn't be there. Zero it out.
all_MAT(end-(fmin/f_res):end, :) = 0;

% Apply amplitude modulation
modMAT = repmat(sig_am, [n_freqs_total, 1]);
ground_truth = all_MAT .* modMAT;

% Add some gaussian smoothing to the ground truth matrix'
ground_truth = imgaussfilt(ground_truth, sigma);

% Gaussian has changed range of values. Rescale to 0-1.
ground_truth = rescale(ground_truth);

% Add silence to start and end
ground_truth = [zeros(n_freqs_total, silence*fs), ground_truth, zeros(n_freqs_total, silence*fs)];

% Convert magnitude to power (W)
ground_truth = ground_truth.^2;

end