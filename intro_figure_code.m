clear

f1 = 100;           % start freq (Hz)
f2 = 40;            % end freq (Hz)
phi_carrier = 0;    % Initial phase of carrier waveform (degrees)

fs = 500; % sample rate
dt = 1 / fs; % seconds per sample
duration_sweep = 2; % seconds
duration_silence = 0.5; % seconds
duration_total = duration_sweep + 2* duration_silence;
samps_total = duration_total * fs;
t_vec = (0 : dt : duration_sweep - dt)'; % time vector for sweep
t_vec_final = linspace(0, duration_total, samps_total);

nyq = fs/2; % nyquist frequency

% Generate Sine Sweep (Carrier Signal)
signal = chirp(t_vec, f1, t_vec(end), f2, 'linear', phi_carrier); %, phi_carrier

% Pad signal and mod with zeros to insert silence at start & end.
signal = [zeros(duration_silence*fs, 1); signal; zeros(duration_silence*fs, 1)];

% Magnitude Spectrum
f_signal = fft(signal);
signal_mag = abs(f_signal(1:length(f_signal)/2+1));
fvec = linspace(0, nyq, length(signal_mag));

% STFT Spectrogram
win = 150; % spectrogram time window (samples)
overlap = round(win * 0.75);
n_fft = 2048;           % spectrogram FFT length (samples)
[s, f, t] = spectrogram(signal, win, overlap, n_fft, fs, "yaxis", 'reassigned');
s = abs(s)';

% Plotting
figure (1)
tl = tiledlayout(2, 2);

nexttile
plot(t_vec_final, signal)
ylim([-1.1 1.1])
xlim([0 duration_total])
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.FontSize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; 
ttl.HorizontalAlignment = 'left';  
ylabel('Amplitude (Arbitrary)')
xlabel('Time (Seconds)')
set(gca, 'fontsize', 12)
grid on

nexttile
plot(fvec, signal_mag)
ylim([0 60])
xlim([0 150])
tt2 = title('b');
tt2.FontWeight = 'bold';
tt2.FontSize = 12;
tt2.Units = 'Normalize'; 
tt2.Position(1) = 0; 
tt2.HorizontalAlignment = 'left';  ylabel('Magnitude (Arbitrary)')
xlabel('Frequency (Hz)')
set(gca, 'fontsize', 12)
grid on

nexttile ([1,2])
surf(f, t, s, EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Magnitude (arbitrary)');
tt3 = title('c');
tt3.FontWeight = 'bold';
tt3.FontSize = 12;
tt3.Units = 'Normalize'; 
tt3.Position(1) = 0; 
tt3.HorizontalAlignment = 'left';  axis on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Magnitude (arbitrary)');
xlim([20 125])
set(gca, 'fontsize', 12, XDir="reverse",  View=[90 90])

set(gcf, 'Position', [100 100 800 450])
saveas(gcf,'waveform_magspect_stft','svg')
