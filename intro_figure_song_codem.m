clear

% Globals
fs = 250; % sample rate
dt = 1 / fs; % seconds per sample
nyq = fs/2; % nyquist frequency
n_fft = 2048;           % spectrogram FFT length (samples)


% unit 1
f1_1 = 60;           % start freq (Hz)
f1_2 = 30;            % end freq (Hz)
d1 = 2;           % duration seconds
phi_1 = 0;    % Initial phase of carrier waveform (degrees)
t_s1 = (0 : dt : d1- dt)';
unit1 = chirp(t_s1, f1_1, t_s1(end), f1_2, 'logarithmic', phi_1); %, phi_carrier

% unit 2
f2_1 = 33;           % start freq (Hz)
f2_2 = 39;            % end freq (Hz)
f2am1 = 3;              % am start freq
f2am2 = 5;              % am end freq
d2 = 2.2;           % duration seconds
phi_1 = 0;      % Initial phase of carrier waveform (degrees)
phi_am2 = -90;  % Initial phase of AM waveform (degrees)
t_s2 = (0 : dt : d2- dt)';
unit2 = chirp(t_s2, f2_1, t_s2(end), f2_2, 'linear', phi_1); %, phi_carrier
am_unit2 = rescale(sign(chirp(t_s2, f2am1, t_s2(end), f2am2, 'logarithmic', phi_am2))); %, phi_carrier))
unit2 = unit2 .* am_unit2;

% unit 3
f3_1 = 90;           % start freq (Hz)
f3_2 = 30;            % end freq (Hz)
d3 = 0.6;           % duration seconds
phi_3 = 0;    % Initial phase of carrier waveform (degrees)
t_s3 = (0 : dt : d3- dt)';
unit3 = 1.2 .* chirp(t_s3, f3_1, t_s3(end), f3_2, 'logarithmic', phi_1); %, phi_carrier
unit3 = repmat([unit3; zeros(0.4*fs, 1)], 6, 1);

% Silent bits
s1d = 1.2; % seconds
s1 = zeros(s1d*fs, 1);
s2d = 0.6;
s2 = zeros(s2d*fs, 1);
s3d = 1.4;
s3 = zeros(s3d*fs, 1);
s4d = 1;
s4 = zeros(s4d*fs, 1);

% Construct full song
signal = [s1; unit1; s2; unit2; s3; unit3; s4];

% Normalize
signal = signal ./ max(abs(signal));

% Vectors & Counts
samps_total = length(signal);
duration_total = samps_total / fs;
t_vec_final = linspace(0, duration_total, samps_total);

% Magnitude Spectrum
f_signal = fft(signal);
signal_mag = abs(f_signal(1:length(f_signal)/2+1));
fvec = linspace(0, nyq, length(signal_mag));

% STFT Spectrogram 1
win1 = 200; % spectrogram time window (samples)
overlap1 = round(win1 * 0.9);
[s1, f1, t1] = spectrogram(signal, win1, overlap1, n_fft, fs, "yaxis", 'reassigned'); % 
s1 = abs(s1);
s1 = s1 ./ max(s1, [], 'all');

% STFT Spectrogram 2
win2 = 50; % spectrogram time window (samples)
overlap2 = round(win2 * 0.5);
[s2, f2, t2] = spectrogram(signal, win2, overlap2, n_fft, fs, "yaxis", 'reassigned'); % 
s2 = abs(s2);
s2 = s2 ./ max(s2, [], 'all');

% Plotting
freqlim = [20 125];

figure (1)
tl = tiledlayout(2, 2);

nexttile ([1, 2])
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
grid on

nexttile
TFRplot(t1, f1, s1, freqlim, [t1(1) t1(end)])
tt3 = title('b');
tt3.FontWeight = 'bold';
tt3.FontSize = 12;
tt3.Units = 'Normalize'; 
tt3.Position(1) = 0; 
tt3.Position(2) = 1; 
tt3.HorizontalAlignment = 'left';  axis on

nexttile
TFRplot(t2, f2, s2, freqlim, [t2(1) t2(end)])
tt3 = title('c');
tt3.FontWeight = 'bold';
tt3.FontSize = 12;
tt3.Units = 'Normalize'; 
tt3.Position(1) = 0; 
tt3.Position(2) = 1; 
tt3.HorizontalAlignment = 'left';  axis on


set(gcf, 'Position', [100 100 800 450])
saveas(gcf,'struct_sig_example_stft','svg')
