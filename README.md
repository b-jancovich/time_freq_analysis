# time_freq_analysis

This repository contains software for the demonstration and evaluation of the Superlets Transform, aka. SLT, (Moca et al. 2021, Transylvanian Institute Of Neuroscience), and comparison of this new time-frequency analysis method with conventional methods such as the short-time Fourier transform. The intended application is in the analysis and identification of animal vocalisations. 
The repository contains MATLAB source code and a compiled MATLAB Runtime application called SLT_Demo_Tool.

-	If you would like to evaluate the SLT for yourself without needing to write any code, and with no software dependencies, download the SLT_Demo_Tool folder. This folder contains an installer for windows. Mac & Linux versions are not currently available. Note the installer will automatically download and install MATLAB Runtime if it is not already on your system.

-	If you would like to evaluate the code used to write this tool, or any of the other code used in the accompanying paper, download the entire repository. 

To run the MATLAB source code, the following software dependencies must be installed:
-	MATLAB 2022b
-	MATLAB Signal Processing Toolbox
-	MATLAB Image Processing Toolbox
-	MATLAB Wavelet Toolbox
-	Superlet Transform “nfaslt.m” Available at: https://github.com/TransylvanianInstituteOfNeuroscience/Superlets
