% Make Figures for Research Proposal

%% Simplest Visualisations - Waveform & Magresp - OK for Steady State Signals
close all
clear 
clc

f1 = 40;      % frequency (Hz)
duration = 2; % seconds
fs = 1000;
dt = 1 / fs;                   % seconds per sample
t = (0 : dt : duration - dt)';     % time vector
nyq = fs/2;

signal = sin(2 * pi * f1 * t); % Trigonometric sinusoid function
signal = signal + rescale(rand(size(signal)), -0.5, 0.5);

tvec = linspace(0, length(signal)/fs, length(signal));

f_signal = fft(signal);
signal_mag = abs(f_signal(1:length(f_signal)/2));
fvec = linspace(0, nyq, length(signal_mag));

figure (1)
tl = tiledlayout(1, 2);
nexttile
plot(tvec, signal);
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel('Amplitude (Arbitrary)')
xlabel('Time (Seconds)')
ylim([-1.5, 1.5])
xlim([0 0.5])
set(gca, 'fontsize', 12, FontName='Calibri')
grid on
nexttile
plot(fvec, signal_mag)
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel('Magnitude (Arbitrary)')
xlabel('Frequency (Hz)')
set(gca, 'fontsize', 12, FontName='Calibri')
grid on

tl.TileSpacing = 'compact';
tl.Padding = 'loose';
set(gcf, 'Position', [100 100 1000 300])
saveas(gcf,'waveform_magresp_steadystatesig','svg')

%% Simplest Visualisations - Waveform & Magresp - Not Good for Time Variant Signals

clear

f1 = 100;      % start freq (Hz)
f2 = 40;       % end freq (Hz)
duration = 1.6; % seconds
silence = 0.2; % seconds
fs = 1000; % sample rate
dt = 1 / fs; % seconds per sample
t = (0 : dt : duration - dt)'; % time vector
nyq = fs/2; % nyquist frequency
signal = chirp(t, f1, t(end), f2); % do sweep
signal = [zeros(silence*fs, 1); signal; zeros(silence*fs, 1)]; % insert silence
signal = signal + rescale(rand(size(signal)), -0.5, 0.5); % insert noise

tvec = linspace(0, length(signal)/fs, length(signal));

f_signal = fft(signal);
signal_mag = abs(f_signal(1:length(f_signal)/2));
fvec = linspace(0, nyq, length(signal_mag));

% Waveform & Magresp for Time Varying Signal - not good.
figure (2)
tl = tiledlayout(1, 2);
nexttile
plot(tvec, signal)
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel('Amplitude (Arbitrary)')
xlabel('Time (Seconds)')
set(gca, 'fontsize', 12, FontName='Calibri')
grid on
nexttile
plot(fvec, signal_mag)
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel('Magnitude (Arbitrary)')
xlabel('Frequency (Hz)')
set(gca, 'fontsize', 12, FontName='Calibri')
grid on
tl.TileSpacing = 'compact';
tl.Padding = 'loose';
set(gcf, 'Position', [100 100 1000 300])
saveas(gcf,'waveform_magresp_timevarysig','svg')

win = ceil(fs/20); % spectrogram time window (samples)
n_fft = 2048;           % spectrogram FFT length (samples)
[s, f, t] = spectrogram(signal, win, win/2, n_fft, fs, "yaxis", 'reassigned');

% compute power from complex STFT data
p = abs(s)' .^2;

% Magresp and STFT Spectrogram - intro TF analysis
figure (3)
t1=tiledlayout(1, 2);
nexttile
plot(signal_mag, fvec)
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Magnitude (Arbitrary)')
ylabel('Frequency (Hz)')
ylim([10 150])
set(gca, 'fontsize', 12, FontName='Calibri')
grid on
nexttile
surf(f, t, p, EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim([10 150])
set(gca, 'fontsize', 12, FontName='Calibri',  XDir="reverse",  View=[90 90])
tl.TileSpacing = 'compact';
tl.Padding = 'loose';
set(gcf, 'Position', [100 100 1000 300])
saveas(gcf,'magresp_STFT_timevarysig','svg')

%% Effect of Window Size

clear 

f1 = 100;      % start freq (Hz)
f2 = 40;       % end freq (Hz)
f3 = 10;
duration = 1.6; % seconds
silence = 0.2; % seconds
fs = 1000; % sample rate
dt = 1 / fs; % seconds per sample
t = (0 : dt : duration - dt)'; % time vector
nyq = fs/2; % nyquist frequency
signal = chirp(t, f1, t(end), f2, 'logarithmic'); % do sweep
mod = chirp(t, f3, t(end), f3); % do amp mod signal
signal = signal .* rescale(mod, 0, 1);
signal = [zeros(silence*fs, 1); signal; zeros(silence*fs, 1)]; % insert silence
signal = signal + rescale(rand(size(signal)), -0.5, 0.5); % insert noise

win1 = 200; % spectrogram time window (samples)
win2 = 50; % spectrogram time window (samples)
n_fft = 2048;           % spectrogram FFT length (samples)
tvec = linspace(0, length(signal)/fs, length(signal));

[s1, f1, t1] = spectrogram(signal, win1, win1/2, n_fft, fs, "yaxis", 'reassigned');
[s2, f2, t2] = spectrogram(signal, win2, win2/2, n_fft, fs, "yaxis", 'reassigned');

% compute power from complex STFT data
p1 = (abs(s1)' .^2);
p2 = (abs(s2)' .^2);

freqlim = [10 120];
timelim = [0 2];

figure (4)
tl = tiledlayout(2, 2);
nexttile([1 2])
plot(tvec, signal)
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylabel('Amplitude (Arbitrary)')
xlabel('Time (Seconds)')
ylim([-1.5, 1.5])
xlim([timelim])
set(gca, 'fontsize', 12, FontName='Calibri')
grid on

nexttile
surf(f1, t1, p1, EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", FontName='Calibri', View=[90 90])

nexttile
surf(f2, t2, p2, EdgeColor = 'none', FaceColor='texturemap')
ttl = title('c');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", FontName='Calibri', View=[90 90])

tl.TileSpacing = 'compact';
tl.Padding = 'loose';
set(gcf, 'Position', [100 100 1000 500])
saveas(gcf,'waveform_STFT_50pt_200pt_50OverPerc','svg')
