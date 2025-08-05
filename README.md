This is code to estimate the magnetic fields generated in the FLARE vacuum vessel. There are two steps:
 - calculate the eddy currents in the vacuum vessel for a given PF waveform (using CoilSystem and Industance Solver)
 - calculte the fields resulting from those currents (using Transient Response) 

 All classes and raw data files assume SI units. The sample script shows the main functionality. 