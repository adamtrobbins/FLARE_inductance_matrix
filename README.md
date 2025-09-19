This is code to estimate the magnetic fields generated in the FLARE vacuum vessel. There are two steps:
 - calculate the eddy currents in the vacuum vessel for a given PF waveform (using CoilSystem and Industance Solver)
 - calculte the fields resulting from those currents (using Transient Response) 

 All classes and raw data files assume SI units. The sample script shows the main functionality. 

 This model works reasonably well inside the vessel, less so for the flux leaking through the vessel: while the functional form is reasonable, the magnitude is too low compared to external probe measurements (a problem shared by ANSYS results). 
