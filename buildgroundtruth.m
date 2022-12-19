function [groundtruth_t, groundtruth_f, groundtruth] = buildgroundtruth(fc1, fc2, fam1, fam2, f_vec, t_vec,...
    sigma, duration_sweep, duration_silence, phi_am, fs, f_res, freqscale, f_dir)
% This function generates a time-frequency representation for an amplitude
% modulated signal with known chracteristics, with aspect ratio matching the
% ratio of sizes of arguments t_vec and f_vec. This is intended to serve as the
% ground truth time frequency representation to compare with the TFR
% provided by a time-frequency analysis algorithm like the STFT.
%
% Groundtruth must be built at the same aspect ratio as the TFRs it is
% being compared with, however it cannot be made to match their size, as
% in some cases, the time and/or frequency resolution are too low to
% accurately render the groundtruth due to aliasing. For example, when 
% building stft_longwin_groundtruth, the amplitude modulation is aliased 
% beyond recognition due to the small number of samples (low time resolution).
% The groundtruth is therefore rendered at matching aspect ratio with
% minimum set by length(shortest_dimension) = fs.
%
% The function then models the first 3 upper and lower sideband components 
% created by amplitude modulation. Argument "am_waveshape" selects sine, 
% square or saw wave for AM, and determins the order of sideband components 
% to generate. The function "tfmatgen2.m" is called once for each component.
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
% fs                = Sampling frequency of the calling function (Hz)
% f_res             = The frequency resolution at which to generate the 
%                       groundtruth. ie. the line width in the groundtruth
%                       TFR. (Hz)
% freqscale         = Selects linear ('lin') or logarithmic ('log') 
%                       scaling for the frequency axis. Must match the
%                       scaling of f_vec. (optional, string, default: lin)
% f_dir             = Direction of the frequency vector. 'normal' is low to
%                       high, 'reversed' is high to low. (optional, string, 
%                       default: normal)
%
% Outputs:
% groundtruth       = The matrix containing a theoretically perfect
%                       time-frequency representation of the test signal.
%                       The units are magnitude squared (power), normalized
%                       to a range of (0, 1).
% groundtruth_t     = Time vector of the upscaled ground truth TFR (s)
% groundtruth_f     = Frequency vector of the upscaled ground truth TFR (Hz)
%

%% Init Vectors

% Set default values for optional input arguments if they are not present.

% If less than 12 args are entered, throw error.
if nargin < 12
    assert(nargin >= 12, 'ERROR: NOT ENOUGH INPUT ARGUMENTS')

    % If less than 13 args entered, init defaults for 13 and 14
elseif nargin < 13
    freqscale = 'lin';
    f_dir = 'normal';

    % If less than 14 are entered, first check 13:
elseif nargin < 14
    f_dir = 'normal';

    % If 14 args are entered, check they are valid
elseif nargin == 14

    % If arg 13 is not 'lin' or 'log', use 'lin'.
    if strcmp(freqscale, 'lin') ~= 1 && strcmp(freqscale, 'log') ~= 1
        warning('Warning: invalid selection for freqscale argument. Using default value "lin"')
        freqscale = 'lin';
    end

    % If arg 14 is not 'normal' or 'reverse', use 'normal'.
    if strcmp(f_dir, 'normal') ~= 1 && strcmp(f_dir, 'reverse') ~= 1
        warning('Warning: invalid selection for f_dir argument. Using default value "normal"')
        f_dir = 'normal';
    end

    % If more than 14 args are entered, throw error.
elseif nargin > 14
    assert(nargin > 14, 'ERROR: TOO MANY INPUT ARGUMENTS')
end

% Set frequnecy vector direction
switch f_dir
    case 'normal'
        f_start = min(f_vec);
        f_end = max(f_vec);
    case 'reverse'
        f_start = max(f_vec);
        f_end = min(f_vec);
end

%% Oversampling

% Resample the time and frequency vectors for groundtruth so that the
% aspect ratio is matched to the corresponding TFR, but with sufficient
% resolution in frequency and time to render the true signal.
% Chech for adequate time resolution

if length(t_vec) < (fs*2)
    p = floor((fs*2) / length(t_vec));
    switch freqscale
        case 'lin'
            groundtruth_f = linspace(f_start, f_end, p * length(f_vec))';
        case 'log'
            groundtruth_f = logspace2(f_start, f_end, p * length(f_vec))';
    end
    groundtruth_t = linspace(min(t_vec), max(t_vec), p * length(t_vec))';
elseif length(t_vec) >= (fs*2)
    groundtruth_f = f_vec;
    groundtruth_t = t_vec;
end

% Chech for adequate frequency resolution
if length(f_vec) < (max(f_vec)-min(f_vec))/f_res
    p = floor(fs / length(f_vec));
    switch freqscale
        case 'lin'
            groundtruth_f = linspace(f_start, f_end, p * length(f_vec))';
        case 'log'
            groundtruth_f = logspace2(f_start, f_end, p * length(f_vec))';
    end
    groundtruth_t = linspace(min(t_vec), max(t_vec), p * length(t_vec))';
elseif length(f_vec) >= (max(f_vec)-min(f_vec))/f_res
end

%% Calculate Sweep Start and End Frequencies for Sideband Components

% Set Multiplication factors for sideband frequencies
o1 = 1;
o2 = 3;
o3 = 5;
o4 = 7;
o5 = 9;
o6 = 11;

% Calculate Start and End Frequencies for Upper Sideband Components
usb1_f1 = fc1 + (fam1 * o1); % First order start
usb1_f2 = fc2 + (fam2 * o1); % First order end
usb2_f1 = fc1 + (fam1 * o2); % Second order start
usb2_f2 = fc2 + (fam2 * o2); % Second order end
usb3_f1 = fc1 + (fam1 * o3); % Third order start
usb3_f2 = fc2 + (fam2 * o3); % Third order end
usb4_f1 = fc1 + (fam1 * o4); % Fourth order start
usb4_f2 = fc2 + (fam2 * o4); % Fourth order end
usb5_f1 = fc1 + (fam1 * o5); % Fifth order start
usb5_f2 = fc2 + (fam2 * o5); % Fifth order end
usb6_f1 = fc1 + (fam1 * o6); % Sixth order start
usb6_f2 = fc2 + (fam2 * o6); % Sixth order end

% Calculate Start and End Frequencies for Lower Sideband Components
lsb1_f1 = fc1 - (fam1 * o1); % First order start
lsb1_f2 = fc2 - (fam2 * o1); % First order end
lsb2_f1 = fc1 - (fam1 * o2); % Second order start
lsb2_f2 = fc2 - (fam2 * o2); % Second order end
lsb3_f1 = fc1 - (fam1 * o3); % Third order start
lsb3_f2 = fc2 - (fam2 * o3); % Third order end
lsb4_f1 = fc1 - (fam1 * o4); % Fourth order start
lsb4_f2 = fc2 - (fam2 * o4); % Fourth order end
lsb5_f1 = fc1 - (fam1 * o5); % Fifth order start
lsb5_f2 = fc2 - (fam2 * o5); % Fifth order end
lsb6_f1 = fc1 - (fam1 * o6); % Sixth order start
lsb6_f2 = fc2 - (fam2 * o6); % Sixth order end

%% Build matrices representing carrier and sideband components. 

% Call tfmatgen.m once for each component:

% Carrier
[n_sweep_samps, n_silence_samps, carrier_MAT] = tfmatgen2(fc1, fc2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);

% Upper Sidebands
[~, ~, usb1_MAT] = tfmatgen2(usb1_f1, usb1_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, usb2_MAT] = tfmatgen2(usb2_f1, usb2_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, usb3_MAT] = tfmatgen2(usb3_f1, usb3_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, usb4_MAT] = tfmatgen2(usb4_f1, usb4_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, usb5_MAT] = tfmatgen2(usb5_f1, usb5_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, usb6_MAT] = tfmatgen2(usb6_f1, usb6_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);

% Lower Sidebands
[~, ~, lsb1_MAT] = tfmatgen2(lsb1_f1, lsb1_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, lsb2_MAT] = tfmatgen2(lsb2_f1, lsb2_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, lsb3_MAT] = tfmatgen2(lsb3_f1, lsb3_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, lsb4_MAT] = tfmatgen2(lsb4_f1, lsb4_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, lsb5_MAT] = tfmatgen2(lsb5_f1, lsb5_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);
[~, ~, lsb6_MAT] = tfmatgen2(lsb6_f1, lsb6_f2, groundtruth_f, groundtruth_t, ...
    duration_sweep, duration_silence);

% Combine matrices for all components & scale their magnitudes according to
% halving power for each increasing order number
groundtruth = carrier_MAT + ...
    (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + ...
    (usb3_MAT .* 0.125) + ...
    (usb4_MAT .* 0.0625) + ...
    (usb5_MAT .* 0.0312) + ...
    (usb6_MAT .* 0.0156) + ...
    (lsb1_MAT .* 0.5) + ...
    (lsb2_MAT .* 0.25) + ...
    (lsb3_MAT .* 0.125) + ...
    (lsb4_MAT .* 0.0625) + ...
    (lsb5_MAT .* 0.0312) + ...
    (lsb6_MAT .* 0.0156);

% Thicken lines in ground truth to match resolution set by f_res
approx_f_res = mean(diff(groundtruth_f));
target_f_res = f_res;
r = floor(target_f_res / approx_f_res);
for k = 1:r
    for i = 2:length(groundtruth_f)
        for j = 1:length(groundtruth_t)
            if groundtruth(i, j) ~= 0
                groundtruth(i-1, j) = groundtruth(i, j);
            end
        end
    end
end

%% Amplitude modulation

% AM sweep time vector
t_vec_sweep = linspace(0, duration_sweep, n_sweep_samps);

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

% Add some gaussian smoothing
groundtruth = imgaussfilt(groundtruth, 0.4);

% Apply amplitude modulation
groundtruth = groundtruth .* sig_am;

%% Final cleanup

% Add some gaussian smoothing
groundtruth = imgaussfilt(groundtruth, 0.4);

% Gaussian has changed range of values. Rescale to 0-1.
groundtruth = rescale(groundtruth);