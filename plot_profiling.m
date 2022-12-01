figure (10)
tiledlayout('flow')
nexttile
F = specplot(swlet, t_vec, f_vec, 'superlet', freqlim, timelim);
nexttile
F = specplot(cwlet, t_vec, f_vec, 'contwavlet', freqlim, timelim);