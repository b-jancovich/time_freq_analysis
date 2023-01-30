# time_freq_analysis

This repository contains MATLAB functions and scripts used to evaluate algorithms for the time-frequency analysis of animal vocalisations.
The code was written as part of the following study:

Jancovich, B.A., New Methods for Visualising Animal Sounds, [DOI:preprint], (2023)
Centre for Marine Science and Innovation
School of Biological, Earth and Environmental Sciences
University of New South Wales, Sydney, Australia

To evaluate the code used to conduct this study, please download the entire repository. 
The following MATLAB products must be installed:
-	MATLAB 2022b
-	MATLAB Signal Processing Toolbox
-	MATLAB Image Processing Toolbox
-	MATLAB Wavelet Toolbox

Additionally, the following Github repository must be downloaded and added to the MATLAB Path:
-	Superlets (https://github.com/TransylvanianInstituteOfNeuroscience/Superlets)

The following algorithms are compared:
- Short-time Fourier transform (50pt and 250pt window sizes)
- Continuous wavelet transform
- Fractional adaptive superresolution wavelet transform (aka. Superlets transform, aka. SLT)

"groundtruth_comparison.m"
This script generates a synthetic animal call, runs it through the four time-frequency analysis algorithms, constructs a "ground truth" time-frequency representation, and measures agreement between ground truth and the algorithm outputs. This script was the basis for the "Qualitative Evaluation - Synthetic Animal Call" and the "Qualitative Evaluation - Synthetic Animal Call" sections of this study.

"loadfile_compare_TFRs.m"
This script loads a .wav audio file and performs time-frequency anayses by the four different methods. Analysis parameters are configured by GUI elements. Results are plotted and saved in a user selected location.

References:
Moca, V. V., Bârzan, H., Nagy-Dăbâcan, A. & Mureșan, R. C. Time-frequency super-resolution with superlets. Nat. Commun. 12, 337 (2021).
