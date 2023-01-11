fs = 250;
duration_silence = 1;
duration_sweep = 5;
duration_total = duration_silence+duration_sweep;
fc1 = 50;
fc2 = 30;
fam1 = 2;
fam2 = 7;
fmin = 10;
fmax = fs/2;
fres = 0.2;
phi_am = -90;
sigma = 0.3;
f_vec_total = (fmin:fres:fmax)';
t_vec_total = (0:1/fs:duration_total)';

% Generate groundtruth TFR that is used for plotting only
[groundtruth_t, groundtruth_f, groundtruth] = buildgroundtruth(fc1, fc2, fam1, fam2, ...
    f_vec_total, t_vec_total, sigma, duration_sweep,...
    duration_silence, phi_am, fs, fres);

row_labels = 1:1:length(f_vec_total);
column_labels = 1:1:length(t_vec_total);

imagesc(column_labels, row_labels, groundtruth)
ylabel('Row Number', fontsize=12, FontName='Calibri', FontWeight='bold')
xlabel('Column Number', fontsize=12, FontName='Calibri', FontWeight='bold')
axis on
grid on
set(gca, fontsize=12, FontName='Calibri');
ax = gca;
ax.Layer = 'top';
ax.GridColor = [1 1 1];
ax.GridAlpha = 0.15;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.MinorGridLineStyle = ':';
ax.MinorGridColor = [1 1 1];
ax.MinorGridAlpha = 0.15;
set(gca,'YDir','normal');
c = colorbar;
c.Label.String = ('Magnitude, Arbitrary Units, Normalized');
c.FontName = 'Calibri';
c.FontSize = 12;
saveas(gcf,'Groundtruth_fs250_dursil1_durswp_5_fres0.2_fmin10.svg','svg')

