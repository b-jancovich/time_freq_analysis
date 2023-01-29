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

________________________________________________________________________________

This work is licensed under the Creative Commons Attribution 4.0 International License. 
To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

________________________________________________________________________________

The Superlet Transform function "nfaslt.m" was released by the Transylvanian Institute of Neuroscience under the MIT License
Copyright (c) 2020 Harald Bârzan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
