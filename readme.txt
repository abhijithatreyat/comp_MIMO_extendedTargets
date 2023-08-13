%%README.txt

This code simulates one frame of a MIMO radar system.

The system uses an analog beamformed phased array transmitter and fully digital array receiver. 
Each frame consists of N_beacon subframes, and each subframe contains N_chirp chirps.
On each subframe, the transmitter modulates its spatial domain signal with a different compressive beacon.
The beacon responses in each range-doppler bin are used to estimate angle of departure compressively.
Angle of arrival is directly estimated since the receiver is fully digital (for now, receiver is single element).
%%%%
%%%%
Main function to run is MIMO_Radar_Simulation_script.m
Helper functions:
extract_RD.m - 2D NOMP function for range and doppler
extractSpectrum_MR.m - 1D NOMP function for angle estimation
extract_targets.m - Oversampled FFT extraction of Range and doppler
radar_car_model.m - Converts the point cloud of a car into range, doppler and azimuth values
generate_received_signal.m - Generate the received signal in 1 frame from range, doppler values 
CAD folder contains the point cloud info about the cars in different scenarios


Contributions
This repository is a collaborative effort that brings together the expertise of multiple contributors. The following authors have made significant contributions:

Abhijith Atreya: Owner of the repository. Wrote the core simulation script, which orchestrates the entire simulation process. Additionally, contributed to various helper functions that enable the simulation to model and analyze a MIMO radar system effectively.

Maryam Eslami Rasekh: Original author who introduced the concept of compressive MIMO radar using NOMP algorithm. Her groundbreaking work has led to the development of the extractSpectrum_MR.m function, which leverages 1D NOMP techniques for angle estimation.

Anant Gupta:  The 2D NOMP algorithm integrated into the following code has been inspired by his prior work. The function 'extract_RD.m' extends the concept of NOMP to the range-doppler domain, offering a comprehensive 2D NOMP approach for range and Doppler estimation.
