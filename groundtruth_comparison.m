%% Comparison of Time-Frequency Analysis Algorithms with Ground Truth.

% This script generates a carrier signal (linearly frequency swept sine)
% and an amplitude modulation signal (linearly frequency swept square wave)
% and multiplies them together to obtain the test signal. 
% It then constructs a matrix representing the "Ground Truth" of the test 
% signal's power as a function of both time & frequency. 
% This matrix can be thought of as a theoretically perfect, tranform-free 
% spectrogram of the test signal. The script then runs a series of
% transform-based time-frequency analysis algorithms on the test signal and
% compares the results with ground truth.

clear 
close all
clc

%% User Variables

% Signal Generation Parameters
carrier_freq1 = 50;            % Sine Sweep start frequency (Hz)
carrier_freq2 = 30;            % Sine Sweep end frequency (Hz)
mod_freq1 = 2;          % Amplitude Modulation sweep start frequency (Hz)
mod_freq2 = 6;          % Amplitude Modulation sweep end frequency (Hz)
duration_sweep = 10;      % Sweep Duration
silence = 2;         % silence to pad start and end (seconds)
fs = 250;           % Sampling Frequency

% NOTE: carrier_freq1 must be > than carrier_freq2

% Signal Analysis Parameters
fmin = 0;           % Lowest frquency of interest
fmax = fs/2;        % Highest frequency of interest
f_res = 0.2;        % Frequency resolution (Hz)

%% Generate signals 

% Frequency & samp counts
n_samps_sweep = duration_sweep * fs;
n_freqs_sweep = (carrier_freq1-carrier_freq2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vector
t_vec_sweep = linspace(0, n_samps_sweep/fs, n_samps_sweep);
t_vec_final = linspace(0, n_samps_sweep/fs+2*silence, n_samps_sweep+(2*silence*fs));

% Frequency vectors
f_vec_sweep = linspace(carrier_freq1, carrier_freq2, n_freqs_sweep);
f_vec_total = linspace(fmin, fmax, n_freqs_total);

% Generate Sine Sweep
signal = chirp(t_vec_sweep, carrier_freq1, t_vec_sweep(end), carrier_freq2, 'linear'); 

% Generate Amp Mod Sweep
mod = rescale(sign(chirp(t_vec_sweep, mod_freq1, t_vec_sweep(end), mod_freq2, 'linear'))); 
mod_sil = [zeros(1,silence*fs), mod, zeros(1,silence*fs)];

% Amplitude Modulation
signal = signal .* mod;

% % Ratio of dimensions
r = n_samps_sweep / n_freqs_sweep;    

% Final Test Signal
signal = [zeros(1,silence*fs), signal, zeros(1,silence*fs)];

%% Ground Truth - Components As Vectors

% Generate Ground Truth Frequency Vector for Carrier
carrier = linspace(carrier_freq1, carrier_freq2, n_samps_sweep);

% Generate Ground Truth Frequency Vectors for Sideband Components
usb1_f1 = carrier_freq1 + mod_freq1;
usb1_f2 = carrier_freq2 + mod_freq2;
usb2_f1 = carrier_freq1 + (mod_freq1*3);
usb2_f2 = carrier_freq2 + (mod_freq2*3);
usb3_f1 = carrier_freq1 + (mod_freq1*5);
usb3_f2 = carrier_freq2 + (mod_freq2*5);

lsb1_f1 = carrier_freq1 - mod_freq1;
lsb1_f2 = carrier_freq2 - mod_freq2;
lsb2_f1 = carrier_freq1 - (mod_freq1*3);
lsb2_f2 = carrier_freq2 - (mod_freq2*3);
lsb3_f1 = carrier_freq1 - (mod_freq1*5);
lsb3_f2 = carrier_freq2 - (mod_freq2*5);

usb1 = linspace(usb1_f1, usb1_f2, n_samps_sweep);
usb2 = linspace(usb2_f1, usb2_f2, n_samps_sweep);
usb3 = linspace(usb3_f1, usb3_f2, n_samps_sweep);
lsb1 = linspace(lsb1_f1, lsb1_f2, n_samps_sweep);
lsb2 = linspace(lsb2_f1, lsb2_f2, n_samps_sweep);
lsb3 = linspace(lsb3_f1, lsb3_f2, n_samps_sweep);

% Aapply amplitude modulation
carrier_am = carrier .* mod;
usb1 = usb1 .* mod;
usb2 = usb2 .* mod;
usb3 = usb3 .* mod;
lsb1 = lsb1 .* mod;
lsb2 = lsb2 .* mod;
lsb3 = lsb3 .* mod;

% Convert zeros to nan
carrier_am(find(carrier_am==0)) = NaN;
usb1(find(usb1==0)) = NaN;
usb2(find(usb2==0)) = NaN;
usb3(find(usb3==0)) = NaN;
lsb1(find(lsb1==0)) = NaN;
lsb2(find(lsb2==0)) = NaN;
lsb3(find(lsb3==0)) = NaN;

% Add silence at start and end
carrier_true = [zeros(1,silence*fs-1), (0/0), carrier_am, (0/0), zeros(1,silence*fs-1)];
usb1_true = [zeros(1,silence*fs-1), (0/0), usb1, (0/0), zeros(1,silence*fs-1)];
usb2_true = [zeros(1,silence*fs-1), (0/0), usb2, (0/0), zeros(1,silence*fs-1)];
usb3_true = [zeros(1,silence*fs-1), (0/0), usb3, (0/0), zeros(1,silence*fs-1)];
lsb1_true = [zeros(1,silence*fs-1), (0/0), lsb1, (0/0), zeros(1,silence*fs-1)];
lsb2_true = [zeros(1,silence*fs-1), (0/0), lsb2, (0/0), zeros(1,silence*fs-1)];
lsb3_true = [zeros(1,silence*fs-1), (0/0), lsb3, (0/0), zeros(1,silence*fs-1)];

%% Ground Truth - Time-Freq Matrix

% Build matrices representing carrier and sidebands
% (f1, f2, f_res, nsamps, amp, rowsOUT)
amp = ones(1, n_samps_sweep);
carrier_MAT = sweep_time_freq_mat(carrier_freq1, carrier_freq2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb1_MAT = sweep_time_freq_mat(usb1_f1, usb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb2_MAT = sweep_time_freq_mat(usb2_f1, usb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb3_MAT = sweep_time_freq_mat(usb3_f1, usb3_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb1_MAT = sweep_time_freq_mat(lsb1_f1, lsb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb2_MAT = sweep_time_freq_mat(lsb2_f1, lsb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb3_MAT = sweep_time_freq_mat(lsb3_f1, lsb3_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);

% Combine matrices for all components
all_MAT = carrier_MAT + (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + (usb3_MAT .* 0.125) +...
    (lsb1_MAT .* 0.5) + (lsb2_MAT .* 0.25) +...
    (lsb3_MAT .* 0.125);

% Convert mag to power
all_MAT = all_MAT .^2;

% Apply amplitude modulation
modMAT = repmat(mod,[n_freqs_total, 1]);
all_MAT = all_MAT .* modMAT;

% Add some gaussian smoothing to the ground truth matrix
sigma = 1; % standard deviation of gaussian filter
all_MAT = imgaussfilt(all_MAT, sigma);

% Gaussian has changed range of values. Rescale to 0-1.
all_MAT = rescale(all_MAT);

% Add silence to start and end
all_MAT = flip([zeros(n_freqs_total, silence*fs), all_MAT, zeros(n_freqs_total, silence*fs)], 1);

%% Compute Transforms

% Compute STFT spectrogram
win1 = 200;              % spectrogram time window (samples)
win2 = 50;             % spectrogram time window (samples)
n_fft = n_freqs_total*2;                   % spectrogram FFT length (samples)
[stft_longwin, stftlong_freq, stftlong_time] = spectrogram(signal, win1, win1/2, n_fft, fs, "yaxis", 'reassigned');
[stft_shortwin, stftshort_freq, stftshort_time] = spectrogram(signal, win2, win2/2, n_fft, fs, "yaxis", 'reassigned');
stft_longwin = rescale(abs(stft_longwin)' .^2);    % compute power 
stft_shortwin = rescale(abs(stft_shortwin)' .^2);   % compute power 

% Compute CWT
[wavelet, frq] = cwt(signal, fs, FrequencyLimits=[0 fs/2]);
wavelet = rescale(abs(wavelet), 0, 1);
tms = (0:numel(signal)-1)/fs;

% Compute superlets
superlets = nfaslt(signal, fs, [10, fs/2], n_freqs_total, 3, [10 40], 1);
superlets = rescale(superlets, 0, 1);

%% Plotting 

% Axis limits
freqlim = [10 70];

figure (2)
tiledlayout('flow')

% plot time domain signal
nexttile
line(t_vec_final, signal, color='Blue', LineWidth=0.5);
title('Time Domain Signal', FontWeight='normal', fontsize=12)
axis on
grid on
ylabel('Amplitude (Normalized)');
xlabel('Time (Seconds)');
ylim([-1.5 1.5])
xlim([0 14])

% % plot vector "grund truth" time-frequency representation.
% nexttile
% line(t_vec_final, carrier_true, color='yellow', LineWidth=2);
% hold on
% line(t_vec_final, usb1_true, color='cyan', LineWidth=2);
% line(t_vec_final, lsb1_true, color='cyan', LineWidth=2);
% line(t_vec_final, usb2_true, color='cyan', LineWidth=2);
% line(t_vec_final, lsb2_true, color='cyan', LineWidth=2);
% line(t_vec_final, usb3_true, color='cyan', LineWidth=2);
% line(t_vec_final, lsb3_true, color='cyan', LineWidth=2);
% title('Vector Ground Truth, Amplitude Not Shown', FontWeight='normal', fontsize=12)
% axis on
% grid on
% ylabel('Frequency (Hz)');
% xlabel('Time (Seconds)');
% ylim(freqlim)
% xlim([0 14])
% set(gca,'color', [0 0 0.7]);
% hold off
% legend();
% l = legend('Carrier frequency', 'AM Sidebands');
% l.Color = 'none';
% l.Location = 'southwest';
% l.TextColor = 'white';
% l.EdgeColor = 'white';

% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_final, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('Analytical Ground Truth', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim([0 14])
set(gca, XDir="reverse", View=[90 90])

% Plot STFT with Long Window
nexttile
surf(stftlong_freq, stftlong_time, stft_longwin, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('STFT Spectrogram - 200pt Window', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim([0 14])
set(gca, XDir="reverse", View=[90 90])

% Plot STFT with Short Window
nexttile
surf(stftshort_freq, stftshort_time, stft_shortwin, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('STFT Spectrogram - 50pt Window', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim([0 14])
set(gca, XDir="reverse", View=[90 90])

% Plot CWT
nexttile
surf(frq, tms, wavelet', EdgeColor="none", FaceColor="texturemap")
a = colorbar;
ylabel(a,'Power (Normalized)');
title('CWT Scalogram', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim([0 14])
set(gca, XDir="reverse", View=[90 90])

% Plot Superlets
nexttile
surf(f_vec_total, t_vec_final, superlets', EdgeColor="none", FaceColor="texturemap")
a = colorbar;
ylabel(a,'Power (Normalized)');
title('SWT Scalogram', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim([0 14])
set(gca, XDir="reverse", View=[90 90])
