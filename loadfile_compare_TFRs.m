clear
close all
clc
%% Signal Analysis Variables

% NOTE Maximum frequency of interest is automatically set according
% to the sample rate, ie. fmax = fs/2
% To change the fmax, the script will provde an opportunity (Via a dialog)
% to downsample the audio to a lower sample rate.

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest (Hz)
f_res = 0.1;                % Frequency resolution (Hz)
threshold = 0.6;            % Noise suppression threshold. Higher values result in more suppression.
powerscaling = 'log';       % Plot power as 'lin' (W) or 'log' (dBW)

% STFT Window Sizes
n_fft = 2048;               % spectrogram FFT length (samples)
win_long = 250;             % Window length for longSTFT (samples)
win_short = 50;             % Window length for shortSTFT (samples)
overlap_long = 75;          % Window overlap for longSTFT (%)
overlap_short = 75;         % Window Overlap for shortSTFT (%)

% CWT Parameters (Morse Wavelet)
tbp = 200;              % Time bandwidth product of Morse wavelet.
gamma = 14;             % Symmetry of the wavelet. 3 = symmetric.
vpo = 48;               % Voices per Octave. Must be in range 10 : 48

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 40];    % Interval of superresolution orders
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
            prompt = 'Enter a new, lower sample rate:';
            dlgtitle = 'Sample Rate Conversion';
            definput = {num2str(original_fs)};
            dims = [1 50];
            fs = str2double(inputdlg(prompt,dlgtitle,dims,definput));
            assert(fs < original_fs, 'Error: New Fs must be lower than original Fs');
            signal_downsamp = easySRC(signal, original_fs, fs, fs/2);

            % plot a quick spectrogram to check if sample rate is appropriate
            W = 0.2;
            O = 0.175;
            [s1, f1, t1] = spectrogram(signal, ceil(W*original_fs), ceil(O*original_fs),...
                1024, original_fs, "yaxis");
            [s2, f2, t2] = spectrogram(signal_downsamp,  ceil(W*fs), ceil(O*fs),...
                1024, fs, "yaxis");
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
        fs = original_fs;
end

if exist('signal_downsamp', 'var')
    signal = signal_downsamp;
else
end

% Dialog box to ask whether to plot frequency as lin or log.
freqplotscaling = questdlg('Plot Frequency Axis on Linear or Logarithmic Scale?',...
    'Frequency Axis Scaling', 'lin', 'log', 'log');

% Save some RAM
clearvars si s2 f1 f2 t1 t2 prompt definput dims signal_downsamp

% Select directory to save figures to
here = pwd;
savepath = uigetdir(here, 'Select Figure Save Directory');

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
[cwlet, f_cwt] = cwt(signal, fs, WaveletParameters = [gamma, tbp], FrequencyLimits=[fmin fmax]);

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

% If powerscaling is set to plot power in dB, do log conversions:
switch powerscaling
    case 'log'
        % Convert linear power to Log power (dBW)
        stft_long = 10*log10(stft_long);
        stft_short = 10*log10(stft_short);
        cwlet = 10*log10(cwlet);
        swlet = 10*log10(swlet);
        % Log conversions will result in some -inf values.
        % Replace -inf with the minimum value in each matrix.
        stft_long(isinf(stft_long)) = min(stft_long(~isinf(stft_long)), [], 'all');
        stft_short(isinf(stft_short)) = min(stft_short(~isinf(stft_short)), [], 'all');
        cwlet(isinf(cwlet)) = min(cwlet(~isinf(cwlet)), [], 'all');
        swlet(isinf(swlet)) = min(swlet(~isinf(swlet)), [], 'all');
    case 'lin'
end

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

% % Below is an experimental noise suppressor that models the noise, then
% tries to subtract it from the recording. It doesn't work yet.
% % sample noise test noisesupp
% noise = rescale(rand(n_samps, 1), 0, 1);
% noiseswt = nfaslt(noise, fs, [fmin, fmax], n_freqs, c1, order, mult);
% swlet = swlet ./ max(abs(swlet), [], 'all');
% noiseswt = noiseswt ./ max(noiseswt, [], 'all');
% threshold = 0.7;
% noiseswtscaled = noiseswt .* threshold;
% clean = swlet(swlet<-threshold) - noiseswt;
%
% figure(1)
% imagesc(noiseswtscaled)
% p=colorbar
% figure(2)
% imagesc(swlet)
% p=colorbar
% figure(3)
% imagesc(clean)
% p=colorbar

%% Plotting

paramstoprint = ['fs=', num2str(fs), '_',...
    'fres=', num2str(f_res), '_', 'fmin=', num2str(fmin), '_',...
    'threshold=', num2str(threshold)];

timelim = [0, t_vec(end)];
freqlim = [fmin, fmax];

switch freqplotscaling
    case 'log'
        mid_tick =  round((fmax-(fmax/2)), -2);
        low_tick = round((fmax-fmin)*0.2, -2);
        freqticks = [fmin, low_tick, mid_tick, fmax];
    case 'lin'
        if (fmax-fmin) > 10 && (fmax-fmin) < 125
            freqticks = fmin:10:fmax;
        elseif (fmax-fmin) >125 && (fmax-fmin) < 500
            freqticks = fmin:50:fmax;
        elseif (fmax-fmin) > 500 && (fmax-fmin) < 1500
            freqticks = fmin:100:fmax;
        end
end

% Dynamic STFT Names
stftshort_name = ['STFT, ', num2str(win_short), 'pt. Window'];
stftlong_name = ['STFT, ', num2str(win_long), 'pt. Window'];

% Time domain
figure(1)
t1 = tiledlayout(9, 1);
nexttile([1,1])
plot(t_vec, signal)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('a - Time Domain Signal');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
xlabel('Time (Seconds)');
ylim([-1.2 1.2])
xlim(timelim)

% set(gcf, 'Position', [50 100 1000 350])

% savename = [file(1:end-4), '_Timedomain_Waveform', '.svg'];
% saveas(gcf, fullfile(savepath,savename), 'svg')

% Plot STFT with Short Window
% figure(2)
nexttile([2,1])
surf(stftshort_freq, stftshort_time, stft_short', EdgeColor = 'none', FaceColor='texturemap')
ttl = title(['b - ', stftshort_name, ' Spectrogram']);
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
xticks(freqticks)
set(gca, XDir="reverse", View=[90 90], xscale = freqplotscaling)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;
% set(gcf, 'Position', [50 100 1000 450])
% savename = [file(1:end-4), '_', stftshort_name, '_', paramstoprint, '.svg'];
% saveas(gcf, fullfile(savepath,savename), 'svg')


% Plot STFT with Long Window
% figure(3)
nexttile([2,1])
surf(stftlong_freq, stftlong_time, stft_long', EdgeColor = 'none', FaceColor='texturemap')
ttl = title(['b - ', stftlong_name, ' Spectrogram']);
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
xticks(freqticks)
set(gca, XDir="reverse", View=[90 90], xscale = freqplotscaling)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;
% set(gcf, 'Position', [50 100 1000 450])
% savename = [file(1:end-4), '_', stftlong_name, '_', paramstoprint, '.svg'];
% saveas(gcf, fullfile(savepath,savename), 'svg')

% Plot CWT
% figure(4)
nexttile([2,1])
surf(f_cwt, t_vec, cwlet', EdgeColor="none", FaceColor="texturemap")
ttl = title('d - CWT Scalogram');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
xticks(freqticks)
set(gca, XDir="reverse", View=[90 90], xscale = freqplotscaling)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;
% set(gcf, 'Position', [50 100 1000 450])
% savename = [file(1:end-4), '_CWT_', paramstoprint, '.svg'];
% saveas(gcf, fullfile(savepath,savename), 'svg')

% Plot Superlets
% figure(5)
nexttile([2,1])
surf(f_vec, t_vec, swlet', EdgeColor="none", FaceColor="texturemap")
ttl = title('e - SWT Scalogram');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize';
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
xticks(freqticks)
set(gca, XDir="reverse", View=[90 90], xscale = freqplotscaling)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;
% set(gcf, 'Position', [50 100 1000 450])
% savename = [file(1:end-4), '_SWT_', paramstoprint, '.svg'];
% saveas(gcf, fullfile(savepath,savename), 'svg')

t1.Padding = 'compact';
set(gcf, 'Position', [10 10 1000 2100])

savename = [file(1:end-4), '_', paramstoprint, '.svg'];
saveas(gcf, fullfile(savepath,savename), 'svg')

% clearvars -except path file
%% Optional Playback

play = questdlg('Do you want to playback original audio?',...
    'Play Audio?', 'Yes', 'No', 'No');

switch play
    case 'Yes'
        [signal, original_fs] = audioread(fullfile(path, file));
        sound(signal, original_fs);
    case 'No'
end

% clearvars