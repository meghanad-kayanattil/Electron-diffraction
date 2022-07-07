
# Gas-phase/liquid-phase electron-diffraction data analysis

This repository contains all the data analysis scripts used for the processing, cleaning and analysis of the electron diffraction images which were taken as part of my PhD project in Max-Planck Institute for the structure and dynamics of matter.

Since majority of the work is part of a manuscript, updated notebooks will be made public as soon as possible.

An example of the simulated electron diffraction pattern of CCl4:

![diffraction_modified_scattering_glycerol_isolated](https://github.com/meghanad-kayanattil/Electron-diffraction/blob/main/Itot_ccl4_simulated_kirk.png)

## Prerequisite

Essential libraries required for the simulation are, numpy, pandas, matplolib, pyFAI, there are other standard libraries like os, cv2 etc used for file handling and image processing in different notebooks as well, but for the simulation of electron diffraction pattern they are not required. 

## Usage/Examples

There are a few analysis/simulation which can be done
using the notebooks in this repository. 

The main notebook which can be used for simulating
gas phase/liquid phase electron diffraction pattern 
can be found [here](https://github.com/meghanad-kayanattil/Electron-diffraction/blob/main/Simulation%20of%20electron%20diffraction%20pattern%20Gas.ipynb)

There are examples and proof of principle implimentations
for small molecules included in the notebook

Edits to the scattering functions will be required if the molecule
contains atoms which are not present in the scattering function defenition.
This can be done following the description in the notebook. 

Depending on the detector and diffractometer specifications sub-functions
need to be modified as per the descretion of the user.


## References

 - For fast azimuthal integration I used [pyFAI](https://pyfai.readthedocs.io/en/master/)
 - Kirkland scattering factors for atoms from [here](https://link.springer.com/chapter/10.1007/978-1-4419-6533-2_11)
 - Lobato et. al. scattering factors for atoms from [here](https://onlinelibrary.wiley.com/iucr/doi/10.1107/S205327331401643X)
 - Very fundamnetal discussions on gas-phase electron diffraction can be found [here](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.8.231) 
