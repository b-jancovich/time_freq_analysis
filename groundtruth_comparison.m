%% Comparison of Time-Frequency Analysis Algorithms with Ground Truth.

% This script generates a carrier signal (linearly frequency swept sine)
% and an amplitude modulation signal (linearly frequency swept square)
% and multiplies them together, creating an amplitude modulated test signal.
% It then constructs a matrix that can be thought of as a synthetic
% spectrogram representing the "Ground Truth" of the test signal's
% power as a function of both time & frequency. This spectrogram is
% constructed analytically using knowledge of the test signal's TF
% charactertistics, rather than by transforming or otherwise analysing
% the signal. This ground truth matrix serves as a theoretically perfect
% spectrogram of the test signal. Several transform based TF analyses are
% then performed and compared against the ground truth.

% Note: if "carrier_freq1" x "mod_freq1", "carrier_freq2" x "mod_freq2",
% or any of the products of their first 4 subharmonics reach a negative
% value, this script will break.

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
close
clc

%% User Variables

% STFT FFT size is locked at n_freqs_total*2 = 1250. This results in equal
% number of frequency points between fmin and fmax for all methods. This is set
% using fmin, fmax and fres for other methods.

% NOTE: fc1 and fc2 must NOT be the same frequency.

% Test Signal Parameters:
fc1 = 50;               % Sine Sweep start frequency (Hz) Must be ~= fc2
fc2 = 30;               % Sine Sweep end frequency (Hz) Must be ~= fc1
fam1 = 2;               % Amplitude Modulation sweep start frequency (Hz)
fam2 = 5;               % Amplitude Modulation sweep end frequency (Hz)
duration_sweep = 5;     % Sweep Duration
duration_silence = 1;   % Silence to pad start and end (seconds)
fs = 250;               % Sampling Frequency (Hz)

% Ground Truth Parameters:
sigma = 1;                  % Standard deviation of gaussian filter

% Signal Analysis Parameters
fmin = 10;                  % Lowest frquency of interest
fmax = fs/2;                % Highest frequency of interest
f_res = 0.2;                % Frequency resolution (Hz)
powerscaling = 'lin';       % Plot power as 'lin' (W) or 'log' (dBW)  

% STFT Window Sizes
n_fft = 2*((fs/2)/f_res);   % spectrogram FFT length (samples)
win_long = n_fft/5;             % Window length for long STFT (samples)
win_short = 50;                  % Window length for short STFT (samples)
overlap_long = 75;              % Window overlap % for long STFT
overlap_short = 75;              % Window Overlap % for short STFT

% CWT Parameters
tbp = 120;      % Time bandwidth product of Morse wavelet.
gamma = 3;      % Symmetry of the wavelet. 3 = symmetric.
vpo = 48;       % Voices per Octave. Must be in range 10 : 48

% Superlet Parameters
c1 = 3;             % Initial number of cycles in superlet.
order = [10 50];    % Interval of superresolution orders
mult = 1;           % Multiplicative ('1') or additive ('0') superresolution
dcfilt = 10;        % Cutoff of DC filter (Hz)

%% Generate time & Frequency Vectors

% Frequency & samp counts
n_samps_sweep = duration_sweep * fs;
n_samps_total = n_samps_sweep + (2*(duration_silence*fs));
n_freqs_sweep = (fc1-fc2) / f_res;
n_freqs_total = (fmax-fmin) / f_res;

% Time vectors
t_vec_sweep = linspace(0, n_samps_sweep/fs, n_samps_sweep);
t_vec_total = linspace(0, n_samps_total/fs, n_samps_total);

% Frequency vectors
f_vec_sweep = linspace(fc1, fc2, n_freqs_sweep);
f_vec_total = linspace(fmin, fmax, n_freqs_total);

%% Generate signals

% Generate Sine Sweep (Carrier Signal)
sig_c = chirp(t_vec_sweep, fc1, t_vec_sweep(end), fc2, 'linear');

% Generate Square Sweep (Amplitude Modulator)
sig_am = rescale(sign(chirp(t_vec_sweep, fam1, t_vec_sweep(end), fam2, 'linear')));

% Pad signal and mod with zeros to insert silence at start & end.
sig_c_sil = [zeros(1, duration_silence * fs), sig_c, zeros(1, duration_silence * fs)];
sig_am_sil = [zeros(1, duration_silence * fs), sig_am, zeros(1, duration_silence * fs)];

% Amplitude Modulation
signal = sig_c_sil .* sig_am_sil;

%% Compute Transforms

% Compute short windowed STFT
[stft_shortwin, stft_shortwin_f, stft_shortwin_t] = spectrogram(signal, ...
    win_short, ceil(win_short*(overlap_short/100)), n_fft, fs, "yaxis");

% Compute long windowed STFT
[stft_longwin, stft_longwin_f, stft_longwin_t] = spectrogram(signal, ...
    win_long, ceil(win_long*(overlap_long/100)), n_fft, fs, "yaxis");

% Compute CWT
[cwlet, cwlet_f] = cwt(signal, fs, VoicesPerOctave=vpo, WaveletParameters = [gamma, tbp], FrequencyLimits=[fmin fmax]);

% Compute SLT
slt = nfaslt(signal, fs, [dcfilt, fmax], n_freqs_total, c1, order, mult);

% Matrix sizes
stft_shortwin_size = size(stft_shortwin);
stft_longwin_size = size(stft_longwin);
cwlet_size = size(cwlet);
slt_size = size(slt);



%%
% % STFT Resolutions
% stft_shortwin_fres = (fmax-fmin) / win_short;
% stft_longwin_fres = (fmax-fmin) / win_long;
stft_shortwin_tres = (win_short-(ceil(win_short*(overlap_short/100))))*(1/fs);
stft_longwin_tres = (win_long-(ceil(win_long*(overlap_long/100))))*(1/fs);

[test_s, test_f, test_t] = spectrogram(sig_c, ...
    win_short, ceil(win_short*(overlap_short/100)), n_fft, fs, "yaxis");

% 
% % CWT Resolutions
% idx1 = 1;
% idx2 = 2;
% cwletfrqs = flip(cwlet_f);
% for i = 1:length(cwletfrqs)
%     cwlet_fres_vec(i) = cwletfrqs(idx2) - cwletfrqs(idx1);
%     if idx2 < length(cwletfrqs)
%         idx1 = idx1 +1;
%         idx2 = idx2 +1;
%     end
% end


%% Ground Truth - Time-Freq Matrix

% Because each analysis algorithm will return a matrix with a different
% number of rows (frequencies), and resizing them will corrupt the results,
% each analysis algorithm must be compared with its own ground truth 
% matrix, having the same number of rows.



% Construct groundtruth matrices for each analysis method.
stft_shortwin_groundtruth = buildgroundtruth(fc1, fc2, fam1, fam2, fmin, f_res,...
    n_samps_sweep, stft_shortwin_size(1), sigma, duration_silence, sig_am, fs); % runs ok

stft_longwin_groundtruth = buildgroundtruth(fc1, fc2, fam1, fam2, fmin, f_res,...
    n_samps_sweep, stft_longwin_size(1), sigma, duration_silence, sig_am, fs); % runs ok

cwlet_groundtruth = buildgroundtruth(fc1, fc2, fam1, fam2, fmin, f_res,...
    n_samps_sweep, cwlet_size(1), sigma, duration_silence, sig_am, fs);

slt_groundtruth = buildgroundtruth(fc1, fc2, fam1, fam2, fmin, f_res,...
    n_samps_sweep, slt_size(1), sigma, duration_silence, sig_am, fs);

% test plot
figure(1)
tiledlayout(1,2)
nexttile
imagesc(stft_shortwin_groundtruth)
nexttile
imagesc(abs(stft_shortwin).^2)

figure(2)
tiledlayout(1,2)
nexttile
imagesc(stft_longwin_groundtruth)
nexttile
imagesc(abs(stft_longwin).^2)

figure(3)
tiledlayout(1,2)
nexttile
imagesc(cwlet_groundtruth)
nexttile
imagesc(abs(cwlet).^2)

figure(4)
tiledlayout(1,2)
nexttile
imagesc(slt_groundtruth)
nexttile
imagesc(slt.^2)


%% Unit Convertsions & Normalizations

% Data to be plotted
% Convert algorithm outputs to power (W)
stft_longwin_plot = (abs(stft_longwin).^2);     % spectrogram() returns complex data.
stft_shortwin_plot = (abs(stft_shortwin).^2);   % spectrogram() returns complex data.
cwlet_plot = (abs(cwlet) .^2);                  % cwt() returns complex data.
slt_plot = (slt .^2);                           % nfaslt returns magnitude.

% If powerscaling is set to plot power in dB, do log conversions:
switch powerscaling
    case 'log'
        % Convert linear power (W) to log power (dBW)
        stft_longwin_plot = 10*log10(stft_longwin_plot / 1);
        stft_shortwin_plot = 10*log10(stft_shortwin_plot / 1);
        cwlet_plot = 10*log10(cwlet_plot / 1);
        slt_plot = 10*log10(slt_plot / 1);

        % Log conversions will result in some -inf values.
        % Replace -inf with the minimum value in each matrix.
        stft_longwin_plot(isinf(stft_longwin_plot)) = min(stft_longwin_plot(~isinf(stft_longwin_plot)), [], 'all');
        stft_shortwin_plot(isinf(stft_shortwin_plot)) = min(stft_shortwin_plot(~isinf(stft_shortwin_plot)), [], 'all');
        cwlet_plot(isinf(cwlet_plot)) = min(cwlet_plot(~isinf(cwlet_plot)), [], 'all');
        slt_plot(isinf(slt_plot)) = min(slt_plot(~isinf(slt_plot)), [], 'all');
    case 'lin'
end

% % Normalize to max=1
stft_shortwin_plot = stft_shortwin_plot ./ max((stft_shortwin_plot), [], 'all');
stft_longwin_plot = stft_longwin_plot ./ max((stft_longwin_plot), [], 'all');
cwlet_plot = cwlet_plot ./ max((cwlet_plot), [], 'all');
slt_plot = slt_plot ./ max((slt_plot), [], 'all');

% Data for use in measures of error
% Convert algorithm outputs to power (W)
stft_longwin = abs(stft_longwin).^2;    % spectrogram returns complex data.
stft_shortwin = abs(stft_shortwin).^2;  % spectrogram returns complex data.
cwlet = abs(cwlet) .^2;                 % cwt returns complex data.
slt = slt .^2;              % nfaslt returns magnitude.

% Normalize to max = 1
stft_shortwin = stft_shortwin ./ max((stft_shortwin), [], 'all');
stft_longwin = stft_longwin ./ max((stft_longwin), [], 'all');
cwlet = cwlet ./ max((cwlet), [], 'all');
slt = slt ./ max((slt), [], 'all');

%% Compute Error



% Compare similarity to ground truth via Root Mean Square Error:
stft_shortwin_freqerror = mean(rmse(stft_shortwin_resz, all_MAT, 1));
stft_shortwin_timeerror = mean(rmse(stft_shortwin_resz, all_MAT, 2));
stft_shortwin_totalerror = rmse(stft_shortwin_resz, all_MAT, 'all');

stft_longwin_freqerror = mean(rmse(stft_longwin_resz, all_MAT, 2));
stft_longwin_timeerror = mean(rmse(stft_longwin_resz, all_MAT, 1));
stft_longwin_totalerror = rmse(stft_longwin_resz, all_MAT, 'all');

cwlet_freqerror = mean(rmse(cwlet_resz, all_MAT, 2));
cwlet_timeerror = mean(rmse(cwlet_resz, all_MAT, 1));
cwlet_totalerror = mean(rmse(cwlet_resz, all_MAT, 'all'));

slt_freqerror = mean(rmse(slt_resz, all_MAT, 2));
slt_timeerror = mean(rmse(slt_resz, all_MAT, 1));
slt_totalerror = rmse(slt_resz, all_MAT, 'all');

% Compare similarity to ground truth via Structural Similarity Index:
stft_shortwin_SSIM = ssim(stft_shortwin_resz, all_MAT);
stft_longwin_SSIM = ssim(stft_longwin_resz, all_MAT);
cwlet_SSIM = ssim(cwlet_resz, all_MAT);
slt_SSIM = ssim(slt_resz, all_MAT);

% Dynamic STFT Names
stftshort_name = ['STFT, ', num2str(win_short), 'pt. Window']; %, num2str(overlap2), ' % Overlap'
stftlong_name = ['STFT, ', num2str(win_long), 'pt. Window']; % , num2str(overlap1), ' % Overlap'



%% Plot Figure 1 - Time Domain Signal

% Test signal
figure(1)
t1 = tiledlayout(3, 1);
nexttile
plot(t_vec_total, sig_c_sil)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])
nexttile
plot(t_vec_total, sig_am_sil)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-.05 1.5])
nexttile
plot(t_vec_total, signal)
ylabel('Amplitude (arbitrary)', FontWeight='normal');
ttl = title('c');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
xlabel('Time (Seconds)');
set(gca, FontSize=12, FontName='Calibri')
ylim([-1.5 1.5])

set(gcf, 'Position', [50 100 1000 500])
saveas(gcf,'Timedomain_Test_signal','svg')

%% Plot Figure 2 & 3 - Error

% Collate Error data for plotting
xlabels1 = categorical({'Freq Axis', 'Time Axis', 'Total RMSE'});
xlabels1 = reordercats(xlabels1,{'Freq Axis', 'Time Axis', 'Total RMSE'});
ydata1 = [stft_shortwin_freqerror, stft_longwin_freqerror, cwlet_freqerror, slt_freqerror;
    stft_shortwin_timeerror, stft_longwin_timeerror, cwlet_timeerror, slt_timeerror;
    stft_shortwin_totalerror, stft_longwin_totalerror, cwlet_totalerror, slt_totalerror];
xlabels2 = categorical({[num2str(win_short), 'pt ', 'STFT'], [num2str(win_long), 'pt ', 'STFT'], 'CWT', 'SWT'});
xlabels2 = reordercats(xlabels2,{[num2str(win_short), 'pt ', 'STFT'], [num2str(win_long), 'pt ', 'STFT'], 'CWT', 'SWT'});
ydata2 = [stft_shortwin_SSIM, stft_longwin_SSIM, cwlet_SSIM, slt_SSIM];

% Plot label rounding
np1 = 3;
np2 = 3;
np3 = 2;
np4 = 4;
np5 = 3;
np6 = 3;
np7 = 3;

% Plot RMSE
figure(2)
b1 = bar(xlabels1, ydata1);
xtips1 = b1(1).XEndPoints;
ytips1 = b1(1).YEndPoints;
labels1 = string(round(b1(1).YData, np1));
text(xtips1,ytips1,labels1,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips2 = b1(2).XEndPoints;
ytips2 = b1(2).YEndPoints;
labels2 = string(round(b1(2).YData, np1));
text(xtips2,ytips2,labels2,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips3 = b1(3).XEndPoints;
ytips3 = b1(3).YEndPoints;
labels3 = string(round(b1(3).YData, np1));
text(xtips3,ytips3,labels3,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
xtips4 = b1(4).XEndPoints;
ytips4 = b1(4).YEndPoints;
labels4 = string(round(b1(4).YData, np1));
text(xtips4,ytips4,labels4,'HorizontalAlignment','left',...
    'VerticalAlignment','middle','Rotation',90,'FontSize',12, FontName='Calibri')
ylim([0 0.14])
lg = legend(stftshort_name, stftlong_name, 'CWT', 'SWT');
lg.Location = 'Northwest';
ylabel 'RMSE re. Ground Truth'
grid on
set(gca, FontSize=12, FontName='Calibri')
set(gcf, 'Position', [50 50 1000 300])
saveas(gcf,'Final_methods_analytical_RMSE_TF','svg')

% Plot SSI
figure(3)
b2 = bar(xlabels2, ydata2);
xtips1_2 = b2(1).XEndPoints;
ytips1_2 = b2(1).YEndPoints;
labels1_2 = string(round(b2(1).YData, np2));
text(xtips1_2,ytips1_2,labels1_2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
ylabel 'SSI re. Ground Truth'
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
ylim([0 (round(max(ydata2)*1.3, np2-2))])
grid on
set(gca, FontSize=12, FontName='Calibri')
xtickangle(90)

set(gcf, 'Position', [50 100 350 300])
saveas(gcf,'Final_methods_ERROR_analytical_SSI','svg')


%% Plot More Figures - Time-Freq Representations

% Common Axis limits
freqlim = [10 70];
timelim = [0 (n_samps_total/fs)];

% Init figure
figure (6)
t1 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], FontSize=12, FontName='Calibri')
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
surf(stft_shortwin_f, stft_shortwin_t, stft_shortwin_plot', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], FontSize=12, FontName='Calibri')
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
t1.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_shortSTFT','svg')

figure (7)
t2 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90],  FontSize=12, FontName='Calibri')
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
surf(stft_longwin_f, stft_longwin_t, stft_longwin_plot', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
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
t2.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_longSTFT','svg')

figure (8)
t3 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
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
surf(cwlet_f, t_vec_total, cwlet_plot', EdgeColor="none", FaceColor="texturemap")
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t3.TileSpacing = 'compact';
t3.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_CWT','svg')


figure (9)
t4 = tiledlayout(1,2);
% Plot matrix "grund truth" time-frequency representation.
nexttile
surf(f_vec_total, t_vec_total, all_MAT', EdgeColor = 'none', FaceColor='texturemap')
ttl = title('a');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
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
surf(f_vec_total, t_vec_total, slt_plot', EdgeColor="none", FaceColor="texturemap")
ttl = title('b');
tt1.FontWeight = 'bold';
tt1.fontsize = 12;
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';  
axis on
grid on
xlabel('Frequency (Hz)');
ylabel('Time (Seconds)');
xlim(freqlim)
ylim(timelim)
set(gca, XDir="reverse", View=[90 90], fontsize=12, FontName='Calibri')
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;

t4.TileSpacing = 'compact';
t4.Padding = 'loose';
set(gcf, 'Position', [50 100 1000 450])
saveas(gcf,'Groundtruth_vs_SWT','svg')


t2.TileSpacing = 'compact';
t2.Padding = 'loose';