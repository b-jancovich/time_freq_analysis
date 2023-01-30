# time_freq_analysis

This repository contains MATLAB functions and scripts used to evaluate algorithms for the time-frequency analysis of animal vocalisations.
The code was written as part of the following study:

Jancovich, B.A., New Methods for Visualising Animal Sounds, [DOI:preprint], (2023)
<br>Centre for Marine Science and Innovation
<br>School of Biological, Earth and Environmental Sciences
<br>University of New South Wales, Sydney, Australia

To evaluate the code used to conduct this study, please download the entire repository. 
<br>The following MATLAB products must be installed:
<br>-	MATLAB 2022b
<br>-	MATLAB Signal Processing Toolbox
<br>-	MATLAB Image Processing Toolbox
<br>-	MATLAB Wavelet Toolbox

Additionally, the following Github repository must be downloaded and added to the MATLAB Path:
<br>-	Superlets (https://github.com/TransylvanianInstituteOfNeuroscience/Superlets)

The following algorithms are compared:
<br>- Short-time Fourier transform (50pt and 250pt window sizes)
<br>- Continuous wavelet transform
<br>- Fractional adaptive superresolution wavelet transform (aka. Superlets transform, aka. SLT)

"groundtruth_comparison.m"
<br>This script generates a synthetic animal call, runs it through the four time-frequency analysis algorithms, constructs a "ground truth" time-frequency representation, and measures agreement between ground truth and the algorithm outputs. This script was the basis for the "Qualitative Evaluation - Synthetic Animal Call" and the "Qualitative Evaluation - Synthetic Animal Call" sections of this study. To run, open script, configure the "User Variables" section, hit "Run".

"loadfile_compare_TFRs.m"
<br>This script loads a .wav audio file and performs time-frequency anayses by the four different methods. Analysis parameters are configured by GUI elements. Results are plotted and saved in a user selected location. To run, open script, hit "Run". All user interactions are via GUI dialog boxes.

References:
<br>Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A. & Mureșan, R. C. Time-frequency super-resolution with superlets. Nat. Commun. 12, 337 (2021).
