%% Comparison of Time-Frequency Analysis Algorithms with Ground Truth.

% This script generates a carrier signal (linearly frequency swept sine)
% and an amplitude modulation signal (linearly frequency swept square)
% and multiplies them together, creating an amplitude modulated test signal. 
% It then constructs a matrix representing the "Ground Truth" of the test 
% signal's power as a function of both time & frequency. This matrix is 
% constructed analytically, without the use of transforms, through 
% knowledge of the signal's generating parameters. The matrix, when viewed 
% as a surface plot serves as a theoretically perfect spectrogram of the 
% test signal. The script then runs a series of transform-based 
% time-frequency analysis algorithms on the test signal and compares the 
% results with ground truth.

% Note: if carrier and 

% Methods compared:
% Short Time Fourier Transform (short window)
% Short Time Fourier Transform (long window)
% Continuous Wavelet Transform
% Fractional Adaptive Superresolution Wavelet Transform

% Ben Jancovich, 2022
% Centre for Marine Science and Innovation
% School of Biological, Earth and Environmental Sciences
% University of New South Wales, Sydney, Australia

clear 
clc

%% User Variables

% NOTE: STFT Overlap is locked at 50%
% STFT FFT size is locked at n_freqs_total*2 = 1250. This results in equal
% number of frequency points between fmin and fmax for all methods. This is set
% using fmin, fmax and fres for other methods.

% Test Signal Parameters:
carrier_freq1 = 50;         % Sine Sweep start frequency (Hz)
carrier_freq2 = 35;         % Sine Sweep end frequency (Hz)
mod_freq1 = 2;              % Amplitude Modulation sweep start frequency (Hz)
mod_freq2 = 7;              % Amplitude Modulation sweep end frequency (Hz)
duration_sweep = 5;         % Sweep Duration
silence = 1;                % Silence to pad start and end (seconds)
fs = 250;                   % Sampling Frequency (Hz)

% Ground Truth Parameters:
sigma = 1;                  % Standard deviation of gaussian filter

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest
fmax = fs/2;                % Highest frequency of interest
f_res = 0.2;                % Frequency resolution (Hz)

% STFT Window Sizes
win1 = 200;                 % Spectrogram time window (samples)
win2 = 50;                  % Spectrogram time window (samples)

% CWT Parameters
tbp = 100;      % Time bandwidth product of Morse wavelet. Scalar, >= 3, <= 120.
vpo = 20;       % Voices per Octave. Int, >= 1, <= 48.

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 50];    % Interval of superresolution orders
mult = 1;           % Multiplicative ('1') or additive ('0') superresolution 
fmin_swt = 10;      % Cutoff of DC filter (Hz)

%% Generate time & Frequency Vectors

% Frequency & samp counts
n_samps_sweep = duration_sweep * fs;
n_samps_total = n_samps_sweep + (2*(silence*fs));
n_freqs_sweep = (carrier_freq1-carrier_freq2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vectors
t_vec_sweep = linspace(0, n_samps_sweep/fs, n_samps_sweep);
t_vec_total = linspace(0, n_samps_total/fs, n_samps_total);

% Frequency vectors
f_vec_sweep = linspace(carrier_freq1, carrier_freq2, n_freqs_sweep);
f_vec_total = linspace(fmin, fmax, n_freqs_total);

%% Generate signals 

% Generate Sine Sweep (Carrier Signal)
signal = chirp(t_vec_sweep, carrier_freq1, t_vec_sweep(end), carrier_freq2, 'linear'); 

% Generate Square Sweep (Amplitude Modulator)
mod = rescale(sign(chirp(t_vec_sweep, mod_freq1, t_vec_sweep(end), mod_freq2, 'linear'))); 

% Pad signal and mod with zeros to insert silence at start & end.
signal_sil = [zeros(1, silence * fs), signal, zeros(1, silence * fs)];
mod_sil = [zeros(1, silence * fs), mod, zeros(1, silence * fs)];

% Amplitude Modulation
signal = signal_sil .* mod_sil; 

%% Ground Truth - Time-Freq Matrix

% Calculate Start and End Frequencies for Upper Sideband Components
usb1_f1 = carrier_freq1 + mod_freq1;
usb1_f2 = carrier_freq2 + mod_freq2;
usb2_f1 = carrier_freq1 + (mod_freq1 * 3);
usb2_f2 = carrier_freq2 + (mod_freq2 * 3);
usb3_f1 = carrier_freq1 + (mod_freq1 * 5);
usb3_f2 = carrier_freq2 + (mod_freq2 * 5);

% Calculate Start and End Frequencies for Lower Sideband Components
lsb1_f1 = carrier_freq1 - mod_freq1;
lsb1_f2 = carrier_freq2 - mod_freq2;
lsb2_f1 = carrier_freq1 - (mod_freq1 * 3);
lsb2_f2 = carrier_freq2 - (mod_freq2 * 3);
lsb3_f1 = carrier_freq1 - (mod_freq1 * 5);
lsb3_f2 = carrier_freq2 - (mod_freq2 * 5);

%% Error handling
sideband_freqs = [usb1_f1, usb1_f2, usb2_f1, usb2_f2, usb3_f1, usb3_f2,...
    lsb1_f1, lsb1_f2, lsb2_f1, lsb2_f2, lsb3_f1, lsb3_f2];
assert(all(sideband_freqs >= 0), ['ERROR: THE PRODUCT OF CARRIER AND ' ...
    'MODULATOR FREQUENCIES YOU HAVE ENTERED HAS RESULTED IN A FREQUENCY ' ...
    'COMPONENT WITH NEGATIVE FREQUENCY. THIS IS NOT SUPPORTED. PLEASE ' ...
    'INCREASE LOWEST CARRIER OR MODULATOR FREQUENCY'])

% Build matrices representing carrier and sidebands
% (f1, f2, f_res, nsamps, amp, rowsOUT)
amp = ones(1, n_samps_sweep);

carrier_MAT = tfmatgen(carrier_freq1, carrier_freq2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb1_MAT = tfmatgen(usb1_f1, usb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb2_MAT = tfmatgen(usb2_f1, usb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
usb3_MAT = tfmatgen(usb3_f1, usb3_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb1_MAT = tfmatgen(lsb1_f1, lsb1_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb2_MAT = tfmatgen(lsb2_f1, lsb2_f2,...
    f_res, n_samps_sweep, amp, n_freqs_total);
lsb3_MAT = tfmatgen(lsb3_f1, lsb3_f2,... 
    f_res, n_samps_sweep, amp, n_freqs_total);

% Combine matrices for all components
all_MAT = carrier_MAT + (usb1_MAT .* 0.5) + ...
    (usb2_MAT .* 0.25) + (usb3_MAT .* 0.125) +...
    (lsb1_MAT .* 0.5) + (lsb2_MAT .* 0.25) +...
    (lsb3_MAT .* 0.125);

% Due to fmin being constant=0 in tfmatgen.m, and a variable in this
% script, the ground truth matrix is shifted up in frequency by the value 
% of fmin/f_res+1. Use circshift to shift it back down.
all_MAT = circshift(all_MAT, (fmin/f_res)+1, 1);

% Apply amplitude modulation
modMAT = repmat(mod,[n_freqs_total, 1]);
all_MAT = all_MAT .* modMAT;

% Add some gaussian smoothing to the ground truth matrix
all_MAT = imgaussfilt(all_MAT, sigma);

% Gaussian has changed range of values. Rescale to 0-1.
all_MAT = rescale(all_MAT);

% Add silence to start and end % REMOVE FLIP - This was a hack to fix a bug
% in tfmatgen - the bug is fexed properly
all_MAT = flip([zeros(n_freqs_total, silence*fs), all_MAT, zeros(n_freqs_total, silence*fs)], 1);

%% Compute Transforms

% spectrogram FFT length (samples)
n_fft = n_freqs_total*2;        

% Compute STFT
[stft_longwin, stftlong_freq, stftlong_time] = spectrogram(signal, ...
    win1, win1/2, n_fft, fs, "yaxis");
[stft_shortwin, stftshort_freq, stftshort_time] = spectrogram(signal, ...
    win2, win2/2, n_fft, fs, "yaxis");

% Convert magnitude to power 
stft_longwin = rescale(abs(stft_longwin) .^2);    
stft_shortwin = rescale(abs(stft_shortwin) .^2); 

% Compute CWT
[wavelet, frq] = cwt(signal, fs, FrequencyLimits=[fmin fmax], ...
    TimeBandwidth=tbp, VoicesPerOctave=vpo);
wavelet = rescale(abs(wavelet), 0, 1);
tms = (0:numel(signal)-1)/fs;

% Compute superlets
superlets = nfaslt(signal, fs, [fmin_swt, fmax], n_freqs_total, c1, order, mult);
superlets = rescale(superlets);

%% Compute Error 

% All matrices must be the same size as ground truth
stft_shortwin_resz = imresize(stft_shortwin,[n_freqs_total, n_samps_total]);
stft_longwin_resz = imresize(stft_longwin,[n_freqs_total, n_samps_total]);
wavelet_resz = imresize(wavelet,[n_freqs_total, n_samps_total]);
superlets_resz = imresize(superlets,[n_freqs_total, n_samps_total]);

% Compare similarity to ground truth via Root Mean Square Error:
stft_shortwin_freqerror = mean(rmse(stft_shortwin_resz, all_MAT, 1));
stft_shortwin_timeerror = mean(rmse(stft_shortwin_resz, all_MAT, 2));
stft_shortwin_totalerror = rmse(stft_shortwin_resz, all_MAT, 'all');

stft_longwin_freqerror = mean(rmse(stft_longwin_resz, all_MAT, 1));
stft_longwin_timeerror = mean(rmse(stft_longwin_resz, all_MAT, 2));
stft_longwin_totalerror = rmse(stft_longwin_resz, all_MAT, 'all');

wavelet_freqerror = mean(rmse(wavelet_resz, all_MAT, 1));
wavelet_timeerror = mean(rmse(wavelet_resz, all_MAT, 2));
wavelet_totalerror = mean(rmse(wavelet_resz, all_MAT, 'all'));

superlet_freqerror = mean(rmse(superlets_resz, all_MAT, 1));
superlet_timeerror = mean(rmse(superlets_resz, all_MAT, 2));
superlet_totalerror = rmse(superlets_resz, all_MAT, 'all');

% Compare similarity to ground truth via Structural Similarity Index:
stft_shortwin_SSIM = ssim(stft_shortwin_resz, all_MAT);
stft_longwin_SSIM = ssim(stft_longwin_resz, all_MAT);
wavelet_SSIM = ssim(wavelet_resz, all_MAT);
superlet_SSIM = ssim(superlets_resz, all_MAT);

% Compare similarity to ground truth via Peak Signal To Noise Ratio:
stft_shortwin_PSNR = psnr(stft_shortwin_resz, all_MAT);
stft_longwin_PSNR = psnr(stft_longwin_resz, all_MAT);
wavelet_PSNR = psnr(wavelet_resz, all_MAT);
superlet_PSNR = psnr(superlets_resz, all_MAT);

% Compare similarity to ground truth via IMMSE:
stft_shortwin_immse = immse(stft_shortwin_resz, all_MAT);
stft_longwin_immse = immse(stft_longwin_resz, all_MAT);
wavelet_immse = immse(wavelet_resz, all_MAT);
superlet_immse = immse(superlets_resz, all_MAT);

%% Plot Figure 1 - Error


% Dynamic STFT titles 
stftlong_name = ['STFT, ', num2str(win1), 'pt. Window'];
stftshort_name = ['STFT, ', num2str(win2), 'pt. Window'];

% Collate Error data for plotting
xlabels1 = categorical({'RMSE - Frequency', 'RMSE - Time', 'Total RMSE'});
xlabels1 = reordercats(xlabels1,{'RMSE - Frequency', 'RMSE - Time', 'Total RMSE'});
ydata1 = [stft_shortwin_freqerror, stft_longwin_freqerror, wavelet_freqerror, superlet_freqerror;
    stft_shortwin_timeerror, stft_longwin_timeerror, wavelet_timeerror, superlet_timeerror;
    stft_shortwin_totalerror, stft_longwin_totalerror, wavelet_totalerror, superlet_totalerror];
xlabels2 = categorical({'20pt STFT', '200pt STFT', 'CWT', 'SWT'});
xlabels2 = reordercats(xlabels2,{'20pt STFT', '200pt STFT', 'CWT', 'SWT'});
ydata2 = [stft_shortwin_SSIM, stft_longwin_SSIM, wavelet_SSIM, superlet_SSIM];
ydata3 = [stft_shortwin_PSNR, stft_longwin_PSNR, wavelet_PSNR, superlet_PSNR];
ydata4 = [stft_shortwin_immse, stft_longwin_immse, wavelet_immse, superlet_immse];

% Plot label rounding
np1 = 3;
np2 = 3;
np3 = 2;
np4 = 4;

% Init figure
figure(1)
t1 = tiledlayout('flow');

% Plot RMSE
nexttile
b1 = bar(xlabels1, ydata1);
xtips1 = b1(1).XEndPoints;
ytips1 = b1(1).YEndPoints;
labels1 = string(round(b1(1).YData, np1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
xtips2 = b1(2).XEndPoints;
ytips2 = b1(2).YEndPoints;
labels2 = string(round(b1(2).YData, np1));
text(xtips2,ytips2,labels2,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
xtips3 = b1(3).XEndPoints;
ytips3 = b1(3).YEndPoints;
labels3 = string(round(b1(3).YData, np1));
text(xtips3,ytips3,labels3,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom')
xtips4 = b1(4).XEndPoints;
ytips4 = b1(4).YEndPoints;
labels4 = string(round(b1(4).YData, np1));
text(xtips4,ytips4,labels4,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom')
ylim([0 0.15])
lg = legend(stftshort_name, stftlong_name, 'CWT', 'SWT');
lg.Location = 'Northwest';
ylabel 'RMSE re. Ground Truth'
title('Root Mean Squared Error in Time, Freq & Total', FontWeight='bold', fontsize=12)
set(gca, 'fontsize', 12)
grid on

% Plot Structural Similarity Index
nexttile
b2 = bar(xlabels2, ydata2);
xtips1_2 = b2(1).XEndPoints;
ytips1_2 = b2(1).YEndPoints;
labels1_2 = string(round(b2(1).YData, np2));
text(xtips1_2,ytips1_2,labels1_2,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
ylabel 'SSI re. Ground Truth'
title('Structural Similarity Index', FontWeight='bold', fontsize=12)
ylim([0 1])
set(gca, 'fontsize', 12)
grid on

% Plot PSNR
nexttile
b3 = bar(xlabels2, ydata3);
xtips1_3 = b3(1).XEndPoints;
ytips1_3 = b3(1).YEndPoints;
labels1_3 = string(round(b3(1).YData, np3));
text(xtips1_3,ytips1_3,labels1_3,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
ylabel 'PSNR re. Ground Truth'
title('Peak Signal to Noise Ratio', FontWeight='bold', fontsize=12)
ylim([0 30])
set(gca, 'fontsize', 12)
grid on

% Plot MSE
nexttile
b4 = bar(xlabels2, ydata4);
xtips1_4 = b4(1).XEndPoints;
ytips1_4 = b4(1).YEndPoints;
labels1_4 = string(round(b4(1).YData, np4));
text(xtips1_4,ytips1_4,labels1_4,'HorizontalAlignment','right',...
    'VerticalAlignment','bottom')
ylabel 'MSE re. Ground Truth'
title('Mean Squared Error', FontWeight='bold', fontsize=12)
ylim([0 0.025])
set(gca, 'fontsize', 12)
grid on

t1.TileSpacing = 'compact';
t1.Padding = 'compact';
set(gcf, 'Position', [500 400 800 500])
% saveas(gcf,'time_vs_freq_simplesig','svg')


%% Plot Figure 2 - Time-Freq Representations

% Common Axis limits
freqlim = [10 70];
timelim = [0 (n_samps_total/fs)];

% Init figure
figure (2)
t2 = tiledlayout('flow');

% plot time domain signal
nexttile
line(t_vec_total, signal, color='Blue', LineWidth=0.5);
title('Time Domain Signal', FontWeight='bold', fontsize=12)
axis on
grid on
ylabel('Amplitude (Normalized)');
xlabel('Time (Seconds)');
ylim([-1.5 1.5])
xlim(timelim)
set(gca, 'fontsize', 12)
ax = gca;
ax.Layer = 'bottom';
ax.GridColor = [0 0 0];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [0 0 0];
ax.MinorGridAlpha = 0.15;

% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title('Analytical Ground Truth', FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
set(gca, 'fontsize', 12)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

% Plot STFT with Short Window
nexttile
surf(stftshort_freq, stftshort_time, stft_shortwin', EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title(stftshort_name, FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
set(gca, 'fontsize', 12)
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
surf(stftlong_freq, stftlong_time, stft_longwin', EdgeColor = 'none', FaceColor='texturemap')
a = colorbar;
ylabel(a,'Power (Normalized)');
title(stftlong_name, FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
set(gca, 'fontsize', 12)
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
surf(frq, tms, wavelet', EdgeColor="none", FaceColor="texturemap")
a = colorbar;
ylabel(a,'Power (Normalized)');
title('CWT Scalogram', FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
set(gca, 'fontsize', 12)
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
surf(f_vec_total, t_vec_total, superlets', EdgeColor="none", FaceColor="texturemap")
a = colorbar;
ylabel(a,'Power (Normalized)');
title('SWT Scalogram', FontWeight='bold', fontsize=12)
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
zlabel('Power (arbitrary)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90])
set(gca, 'fontsize', 12)
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t2.TileSpacing = 'compact';
t2.Padding = 'compact';
set(gcf, 'Position', [500 400 800 500])
% saveas(gcf,'time_vs_freq_simplesig','svg')