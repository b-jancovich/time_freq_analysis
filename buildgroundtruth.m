function groundtruth = buildgroundtruth(fc1, fc2, fam1, fam2, f_vec, t_vec,...
    sigma, duration_sweep, duration_silence, phi_am)
% This function generates a time-frequency representation for an amplitude 
% modulated signal with known chracteristics, with size matching the 
% sizes of arguments t_vec and f_vec. This is intended to serve as the
% ground truth time frequency representation to compare with the TFR 
% provided by a time-frequency analysis algorithm like the STFT.
% 
% It will model the first 3 upper and lower sideband components created by 
% the amplitude modulation. Argument "am_waveshape" selects sine, square 
% or saw wave for AM, and determins the order of sideband components to 
% generate. The function "tfmatgen2.m" is called once for each component. 
% The matrices returned for each call are summed together, then the result 
% is multiplied by a matrix containing the amplitude modulation signal.
% 
% Inputs:
% fc1               = Start Frequency of carrier signal sweep (Hz)
% fc2               = End frequency of carrier signal sweep (Hz)
% fam1              = Start Frequency of AM signal sweep (Hz)
% fam2              = End frequency of AM signal sweep (Hz)
% f_vec             = Frequency vector corresponding to rows of matrix (Hz)
% t_vec             = Time vector corresponding to columns of matric (s)
% sigma             = Standard deviation of gaussian image blur
% duration_sweep    = Duration of carrier and AM sweep signal (s)
% duration_silence  = Duration of silence before and after test signal (s)
% phi_am            = Initial phase of AM wave (degrees)
%
% Outputs:
% groundtruth       = The matrix containing a theoretically perfect 
%                       time-frequency representation of the test signal.
%                       The units are magnitude squared (power), normalized
%                       to a range of (0, 1).
%
%% Calculate Sweep Start and End Frequencies for Sideband Components

% Set Multiplication factors for sideband frequencies
o1 = 1;
o2 = 2;
o3 = 3;

% Calculate Start and End Frequencies for Upper Sideband Components
usb1_f1 = fc1 + (fam1 * o1); % First order start
usb1_f2 = fc2 + (fam2 * o1); % First order end
usb2_f1 = fc1 + (fam1 * o2); % Second order start
usb2_f2 = fc2 + (fam2 * o2); % Second order end
usb3_f1 = fc1 + (fam1 * o3); % Third order start
usb3_f2 = fc2 + (fam2 * o3); % Third order end

% Calculate Start and End Frequencies for Lower Sideband Components
lsb1_f1 = fc1 - (fam1 * o1); % First order start
lsb1_f2 = fc2 - (fam2 * o1); % First order end
lsb2_f1 = fc1 - (fam1 * o2); % Second order start
lsb2_f2 = fc2 - (fam2 * o2); % Second order end
lsb3_f1 = fc1 - (fam1 * o3); % Third order start
lsb3_f2 = fc2 - (fam2 * o3); % Third order end

% Error handling
sideband_freqs = [usb1_f1, usb1_f2, usb2_f1, usb2_f2, usb3_f1, usb3_f2,...
    lsb1_f1, lsb1_f2, lsb2_f1, lsb2_f2, lsb3_f1, lsb3_f2];

% Make sure no frequency is negative
assert(all(sideband_freqs >= 0), ['ERROR: THE PRODUCT OF CARRIER AND ' ...
    'MODULATOR FREQUENCIES YOU HAVE ENTERED HAS RESULTED IN A FREQUENCY ' ...
    'COMPONENT WITH NEGATIVE FREQUENCY. THIS IS NOT SUPPORTED. PLEASE ' ...
    'INCREASE LOWEST CARRIER OR MODULATOR FREQUENCY'])

% Build matrices representing carrier and sideband components
[n_sweep_samps, n_silence_samps, carrier_MAT] = tfmatgen2(fc1, fc2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, usb1_MAT] = tfmatgen2(usb1_f1, usb1_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, usb2_MAT] = tfmatgen2(usb2_f1, usb2_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, usb3_MAT] = tfmatgen2(usb3_f1, usb3_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, lsb1_MAT] = tfmatgen2(lsb1_f1, lsb1_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, lsb2_MAT] = tfmatgen2(lsb2_f1, lsb2_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);
[~, ~, lsb3_MAT] = tfmatgen2(lsb3_f1, lsb3_f2, f_vec, t_vec, ...
    duration_sweep, duration_silence);

% Combine matrices for all components & scale 
groundtruth = carrier_MAT + (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + (usb3_MAT .* 0.125) +...
    (lsb1_MAT .* 0.5) + (lsb2_MAT .* 0.25) +...
    (lsb3_MAT .* 0.125);

%% Amplitude modulation

% t_interval = duration_sweep / n_sweep_samps;

% AM sweep time vector
t_vec_sweep = linspace(duration_silence, duration_sweep+duration_silence, n_sweep_samps);
% t_vec_sweep = duration_silence:t_interval:duration_sweep+duration_silence;

% Construct amplitude modulation vector
sig_am = rescale(sign(chirp(t_vec_sweep, fam1, t_vec_sweep(end), fam2, 'linear', phi_am))); 

% Pad zeros
rem = mod(n_silence_samps, 2);
if rem ~=0
    sig_am = [zeros(1, n_silence_samps+rem), sig_am, zeros(1, n_silence_samps-rem)];
elseif rem == 0
    sig_am = [zeros(1, n_silence_samps), sig_am, zeros(1, n_silence_samps)];
end

% Repeat AM signal across all frequencies (rows)
sig_am = repmat(sig_am, size(groundtruth, 1), 1);

% Throw error if matrices are not the same size
assert(size(sig_am, 1) == size(groundtruth, 1) && size(sig_am, 2) == size(groundtruth, 2),...
    'ERROR: Size of AM signal matrix does not match size of carrier signal matrix');

% Apply amplitude modulation
groundtruth = groundtruth .* sig_am;

% Add some gaussian smoothing to the ground truth matrix'
groundtruth = imgaussfilt(groundtruth, sigma);

% Gaussian has changed range of values. Rescale to 0-1.
groundtruth = rescale(groundtruth);
% 
% % Convert magnitude to power (W)
% groundtruth = groundtruth.^2;

end