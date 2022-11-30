clear all
close all

%% Signal Analysis Variables

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest (Hz)
f_res = 1;                 % Frequency resolution (Hz)
threshold = 0.5;            % Noise suppression threshold. Higher values result in more suppression.

% STFT Window Sizes
n_fft = 2048;               % spectrogram FFT length (samples)
win_long = 250;             % Window length for longSTFT (samples)
win_short = 50;             % Window length for shortSTFT (samples)
overlap_long = 75;          % Window overlap for longSTFT (%)
overlap_short = 75;             % Window Overlap for shortSTFT (%)

% CWT Parameters (Morse Wavelet)
gamma = 3;          % Symmetry of the wavelet. 3 = symmetric.
vpo = 48;           % Voices per Octave. Int, range [10, 48]
tbp = 120;          % Time bandwidth product of wavelet. Int, [3,120]

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 50];    % Interval of superresolution orders
mult = 1;           % Multiplicative (1) or additive (0) superresolution

%% Load & Prep audio file (UI Dialog)
[file, path] = uigetfile('.wav', 'Load WAV file to analyse');
[signal, original_fs] = audioread(fullfile(path, file));

assert(size(signal, 2) == 1, ...
    'Error: this script only supports mono audio files')

signal = signal ./ max(abs(signal));

%% Check if we need to resample again via a quick spectrogram plot

resample = questdlg(['Original Fs is ', num2str(original_fs), ...
    'Hz. Do you want to downsample to a lower Fs?'],...
    'Downsample Audio?', 'Yes', 'No', 'No');

switch resample
    case 'Yes'
        while strcmpi(resample, 'Yes')
            % Set Analysis Sampling Frequency (UI Dialog)
            prompt = ['Downsample audio from ', num2str(original_fs), ' Hz to:'];
            dlgtitle = 'Sample Rate Conversion';
            definput = {num2str(original_fs)};
            dims = [1 50];
            fs = str2double(inputdlg(prompt,dlgtitle,dims,definput));
            assert(fs < original_fs, 'Error: New Fs must be lower than original Fs');
            signal_downsamp = easySRC(signal, original_fs, fs, fs/2);

            % plot a quick spectrogram to check if sample rate is appropriate
            [s1, f1, t1] = spectrogram(signal, 250, 250/2, n_fft, original_fs, "yaxis");
            [s2, f2, t2] = spectrogram(signal_downsamp, 250, 250/2, n_fft, fs, "yaxis");
            s1 = 20*log10(abs(s1) .^ 2);
            s2 = 20*log10(abs(s2) .^ 2);
            s1 = s1 ./ max(abs(s1));
            s2 = s2 ./ max(abs(s2));
            figure(1)
            tiledlayout(2,1)
            nexttile
            surf(f1, t1, s1', EdgeColor = 'none', FaceColor='texturemap')
            xlabel('Frequency (Hz)');
            ylabel('Time (Seconds)');
            zlabel('Power (arb, normalized)');
            title(['Original Fs= ', num2str(original_fs), ' Hz'])
            set(gca, XDir="reverse", View=[90 90])
            nexttile
            surf(f2, t2, s2', EdgeColor = 'none', FaceColor='texturemap')
            xlabel('Frequency (Hz)');
            ylabel('Time (Seconds)');
            zlabel('Power (arb, normalized)');
            title(['Downsampled Fs= ', num2str(fs), ' Hz'])
            set(gca, XDir="reverse", View=[90 90])
            resample = questdlg('Do you need to downsample again?',...
                'Resample again?', 'Yes', 'No', 'No');
            close 1
        end
    case 'No'
end

if exist('signal_downsamp')
        signal = signal_downsamp;
end

clearvars si s2 f1 f2 t1 t2 prompt definput dims signal_downsamp
%% Vectors and Counts

% Highest frequency of interest (Hz)
fmax = fs/2;                

% Frequency & samp counts
n_samps = length(signal);
n_freqs = (fmax-fmin) / f_res;

% Time vector
t_vec = linspace(0, n_samps/fs, n_samps);

% Frequency vector
f_vec = linspace(fmin, fmax, n_freqs);

%% Compute Transforms

% Compute Short Window STFT
[stft_short, stftshort_freq, stftshort_time] = spectrogram(signal, ...
    win_short, ceil(win_short*(overlap_short/100)), n_fft, fs, "yaxis");

% Compute Short Window STFT
[stft_long, stftlong_freq, stftlong_time] = spectrogram(signal, ...
    win_long, ceil(win_long*(overlap_long/100)), n_fft, fs, "yaxis");

% Compute CWT
[cwlet, f_cwt] = cwt(signal, fs, WaveletParameters = [14, 200]);

% Compute superlets
swlet = nfaslt(signal, fs, [fmin, fmax], n_freqs, c1, order, mult);

%% Unit Conversaions

% Have confirmed the following by analysing function code:
% spectrogram() returns complex frequency domain data. 
%           abs(stft_raw) = magnitude. mag^2 for power. 
% cwt() returns complex frequency domain data. 
%           abs(cwt_raw) = magnitude. mag^2 for power. 
% nfaslt() returns magnitude. 
%           mag^2 for power.

% Convert to results to power from complex magnitude
stft_long = abs(stft_long) .^2;
stft_short = abs(stft_short) .^2;
cwlet = abs(cwlet) .^2;
swlet = swlet .^2;

% Convert linear power to Log power (dBW)
stft_long = 10*log10(stft_long);
stft_short = 10*log10(stft_short);
cwlet = 10*log10(cwlet);
swlet = 10*log10(swlet);

% Lots of -inf values in swlet produced by dB conversion.
% Replace -inf with the minimum value in the swlet matrix.
swlet(isinf(swlet)) = min(swlet(~isinf(swlet)), [], 'all');

%% Noise Floor Suppression

% rescale the data to a range of [-threshold : 1]
stft_long = rescale(stft_long, -threshold, 1);
stft_short = rescale(stft_short, -threshold, 1);
cwlet = rescale(cwlet, -threshold, 1);
swlet = rescale(swlet, -threshold, 1);

% Set all values below 0 to 0, effectively silencing low powered pixels.
stft_long(stft_long < 0) = 0;
stft_short(stft_short < 0) = 0;
cwlet(cwlet < 0) = 0;
swlet(swlet < 0) = 0;

%% Plotting
timelim = [0, t_vec(end)];
freqlim = [fmin, fmax];

% Dynamic STFT Names
stftshort_name = ['STFT, ', num2str(win_short), 'pt. Window'];
stftlong_name = ['STFT, ', num2str(win_long), 'pt. Window'];

figure(2)
t1 = tiledlayout(2, 2);

% Plot STFT with Short Window
nexttile
surf(stftshort_freq, stftshort_time, stft_short', EdgeColor = 'none', FaceColor='texturemap')
title(stftshort_name, FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot STFT with Long Window
nexttile
surf(stftlong_freq, stftlong_time, stft_long', EdgeColor = 'none', FaceColor='texturemap')
title(stftlong_name, FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot CWT
nexttile
surf(f_cwt, t_vec, cwlet', EdgeColor="none", FaceColor="texturemap")
title('CWT Scalogram', FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot Superlets
nexttile
surf(f_vec, t_vec, swlet', EdgeColor="none", FaceColor="texturemap")
title('SWT Scalogram', FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t1.TileSpacing = 'compact';
t1.Padding = 'compact';
set(gcf, 'Position', [50 100 1000 600])
savename = [file, 'TFR_comparison'];
saveas(gcf, savename, 'svg')

% Time domain
figure(3)
plot(t_vec, signal)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
title('Time Domain Signal', FontWeight='bold', fontsize=12)
xlabel('Time (Seconds)');
ylim([-1.5 1.5])

set(gcf, 'Position', [50 100 1000 350])
savename = [file, 'Timedomain_signal'];
saveas(gcf, savename, 'svg')