README
SLT_Demo_Tool V1.0

Ben Jancovich
b.jancovich@unsw.edu.au

Center For Marine Science and Innovation
School of Biological, Earth & Environmental Sciences
University of New South Wales

Summary:
This program loads an audio file, computes spectrograms and scalograms using different algorithms & then plots for comparison

Program Description:
This program compares methods for time-frequency analysis of complex, low frequency animal calls. 
The purpose of this program is to demonstrate the Superlets Transform (Moca et al., 2021) as a method for time frequency anaysis in bioacoustics, 
and to illustrate its advantages over conventional methods such as the short-time Fourier transform.

The program is designed to run as a standalone application, with no software dependancies, and no pre-requisite knowledge of coding or signal processing.
It has a rudimentary user interface that allows the user to do the following pre-processing tasks:

	- Load an audio file
	- Trim the file to the time segment of interest
	- Resample the file to the frequency range of interest

The program then asks the user to configure the parameters for each analysis algorithm. 
The recommended value for each parameter is pre-filled in the UI dialog boxes. The algorithms used are:

	- Short-time Fourier transform - implemented using the built-in MATLAB function "spectrogram.m"
	- Continuous wavelet transform - implemented using the built-in MATLAB function "cwt.m"
	- Fractional Adaptive Superresolution Wavelet Transform - implemented using function "nfaslt.m", developed by Moca et al. (2021), 
	released under MIT licence at https://github.com/TransylvanianInstituteOfNeuroscience/Superlets

A UI dialog asks the user to select a location to save the results. 
Two STFT analyses are performed, with different window lengths. 
Results are plotted in a single window, and UI dialog asks whether to play the audio file over the currently selected windows audio device.
A final UI element then asks if the user would like to re-run the analysis on the same (trimmed and resampled) audio file, 
providing an opportunity to test different analysis algorithm parameters.

____________________________________________________________________________
