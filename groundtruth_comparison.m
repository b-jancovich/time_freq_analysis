clear 

% Signal Generation Parameters
stft_freq1 = 50;            % Sine Sweep start frequency (Hz)
stft_freq2 = 30;            % Sine Sweep end frequency (Hz)
modf1 = 2;          % Amplitude Modulation sweep start frequency (Hz)
modf2 = 6;          % Amplitude Modulation sweep end frequency (Hz)
duration = 10;      % Sweep Duration
silence = 2;         % silence to pad start and end (seconds)
fs = 250;           % Sampling Frequency

% f1 must be > than f2

% Signal Analysis Parameters
fmin = 0;           % Lowest frquency of interest
fmax = fs/2;        % Highest frequency of interest
f_res = 0.2;        % Frequency resolution (Hz)

% Frequency & samp counts
n_samps = duration * fs;
n_freqs_sweep = (stft_freq1-stft_freq2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vector
t_vec_sweep = linspace(0, n_samps/fs, n_samps);
t_vec_final = linspace(0, n_samps/fs+2*silence, n_samps+(2*silence*fs));

% Frequency vectors
f_vec_sweep = linspace(stft_freq1, stft_freq2, n_freqs_sweep);
f_vec_total = linspace(fmin, fmax, n_freqs_total);

% Generate Sine Sweep
signal = chirp(t_vec_sweep, stft_freq1, t_vec_sweep(end), stft_freq2, 'linear'); 

% Generate Amp Mod Sweep
mod = rescale(sign(chirp(t_vec_sweep, modf1, t_vec_sweep(end), modf2, 'linear'))); 
mod_sil = [zeros(1,silence*fs), mod, zeros(1,silence*fs)];

% Amplitude Modulation
signal = signal .* mod;

% % Ratio of dimensions
r = n_samps / n_freqs_sweep;    

signal = [zeros(1,silence*fs), signal, zeros(1,silence*fs)];

%% Ground Truth

% Generate Ground Truth Frequency Vectors for Carrier
carrier = linspace(stft_freq1, stft_freq2, n_samps);

% Generate Ground Truth Frequency Vectors for Sideband Components
usb1_f1 = stft_freq1 + modf1;
usb1_f2 = stft_freq2 + modf2;
usb2_f1 = stft_freq1 + (modf1*3);
usb2_f2 = stft_freq2 + (modf2*3);
usb3_f1 = stft_freq1 + (modf1*5);
usb3_f2 = stft_freq2 + (modf2*5);

lsb1_f1 = stft_freq1 - modf1;
lsb1_f2 = stft_freq2 - modf2;
lsb2_f1 = stft_freq1 - (modf1*3);
lsb2_f2 = stft_freq2 - (modf2*3);
lsb3_f1 = stft_freq1 - (modf1*5);
lsb3_f2 = stft_freq2 - (modf2*5);

usb1 = linspace(usb1_f1, usb1_f2, n_samps);
usb2 = linspace(usb2_f1, usb2_f2, n_samps);
usb3 = linspace(usb3_f1, usb3_f2, n_samps);
lsb1 = linspace(lsb1_f1, lsb1_f2, n_samps);
lsb2 = linspace(lsb2_f1, lsb2_f2, n_samps);
lsb3 = linspace(lsb3_f1, lsb3_f2, n_samps);

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

%% Compute Transforms

% Compute STFT spectrogram
win1 = 200;              % spectrogram time window (samples)
win2 = 50;             % spectrogram time window (samples)
n_fft = n_freqs_total*2;                   % spectrogram FFT length (samples)
[stft_longwin, stft_freq1, stft_time1] = spectrogram(signal, win1, win1/2, n_fft, fs, "yaxis", 'reassigned');
[stft_shortwin, stft_freq2, stft_time2] = spectrogram(signal, win2, win2/2, n_fft, fs, "yaxis", 'reassigned');
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

figure (1)
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

% plot "grund truth" time-frequency representation.
nexttile
line(t_vec_final, carrier_true, color='yellow', LineWidth=2);
hold on
line(t_vec_final, usb1_true, color='cyan', LineWidth=2);
line(t_vec_final, lsb1_true, color='cyan', LineWidth=2);
line(t_vec_final, usb2_true, color='cyan', LineWidth=2);
line(t_vec_final, lsb2_true, color='cyan', LineWidth=2);
line(t_vec_final, usb3_true, color='cyan', LineWidth=2);
line(t_vec_final, lsb3_true, color='cyan', LineWidth=2);
title('Ground Truth, Amplitude Not Shown', FontWeight='normal', fontsize=12)
axis on
grid on
ylabel('Frequency (Hz)');
xlabel('Time (Seconds)');
ylim(freqlim)
set(gca,'color', [0 0 0.7]);
hold off
legend();
l = legend('Carrier frequency', 'AM Sidebands');
l.Color = 'none';
l.TextColor = 'white';
l.EdgeColor = 'white';


% Plot STFT with Long Window
nexttile
surf(stft_freq1, stft_time1, stft_longwin, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('STFT Spectrogram - 200pt Window', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
set(gca, XDir="reverse", View=[90 90])

% Plot STFT with Short Window
nexttile
surf(stft_freq2, stft_time2, stft_shortwin, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('STFT Spectrogram - 50pt Window', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
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
set(gca, XDir="reverse", View=[90 90])

% Plot Superlets
nexttile
surf(f_vec_total, t_vec_final, superlets', EdgeColor="none", FaceColor="texturemap")
a = colorbar;
ylabel(a,'Power (Normalized)');
title('Fractional Adaptive Superresolution Wavelet Transform Scalogram', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
set(gca, XDir="reverse", View=[90 90])

%% Compute error

% Do something to get values of mod_sil into matrix at coords in tvec_sweep
% and fvec_sweep. Matlab central answer incoming...



surf(t_vec_sweep, f_vec_sweep, carriermat, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('test', FontWeight='normal', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
set(gca, XDir="reverse", View=[90 90])