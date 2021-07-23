# IMS-HPLC-ethanol-correlation

This toolbox contains functions for the IMS data analysis. 
IMS data can be loaded and pre-processed. 
Furthermore, a correlation to HPLC data is established by different approaches (PCA and peak identification / integration).
The module ethanol_correlation is the main code and IMS_functions contains functions that are executed in the main code.
For reasons of secrecy, two important functions are missing in the toolbox read_mea() and alignment(). 
Also, the IMS data was too large to upload to github.
The missing objects (data and functions) can be obtained from the Mannheim University of Applied Sciences
read_mea() is the function to read the binary .mea data and alignment() is the function to normalize all spectra to the RIP peak. 
The code works if you copy alignment() into IMS_functions and if read_mea() is in the same folder as the other two modules (ethanol_correlation and IMS_functions).
Furthermore, the paths for loading metadata, offline data, and IMS data need to be adjusted in the code.
The folder Offline provides HPLC data corresponding to the IMS data.
