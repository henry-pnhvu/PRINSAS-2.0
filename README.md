# PRINSAS 2.0

## Overview
PRINSAS 2.0 is a GUI-based Python tool for fitting the **Polydisperse Spherical Pore (PDSP) model** to Small Angle Scattering (SAS) data. It extracts pore size distribution, specific surface area (SSA), and porosity from scattering intensity measurements.

## Model Basis
PRINSAS 2.0 is based on the **PDSP model**, which describes the differential scattering cross-section of a system composed of polydisperse spherical pores. The model incorporates Debye’s law for two-phase systems and is mathematically represented as:

$$
\frac{d\Sigma}{d\Omega} (Q) = (\rho_1^* - \rho_2^* )^2  \frac{\phi(1-\phi)}{\bar{V}} \int_{R_{min}}^{R_{max}} V_r^2 f(r) F_{sph} (Qr) dr
$$


where:
- $$\frac{d\Sigma}{d\Omega} (Q) $$ is the differential scattering cross-section,
- $$\rho_1^* $$ and $$\rho_2^* $$ are the scattering length densities of the two phases,
- $$\phi$$ is the total porosity,
- $$\bar{V}$$ is the average pore volume,
- $$V_r$$ is the volume of a pore of radius $$r$$,
- $$f(r)$$ is the pore size distribution,
- $$F_{sph}(Qr)$$ is the form factor for a spherical pore.

## Installation and Execution
### Running the Software
There are two methods to run PRINSAS 2.0:

1. **Using Python (Cross-platform)**:
   - Install Python (version 3.7 or later) and required dependencies.
   - Clone or download this repository.
   - Run the script:
     ```bash
     python run_PRINSAS.py
     ```

2. **Using the Standalone Executable (Windows Only)**:
   - Download the ZIP file from the [GitHub Releases page](https://github.com/henry-pnhvu/PRINSAS-2.0/releases).
   - Extract the ZIP file to a folder on your computer.
   - Locate the `.exe` file inside the extracted folder and run it.  
     **Note:** Windows Defender may require manual approval.

## Further Information
For detailed instructions on using PRINSAS 2.0, including parameter selection and result interpretation, please refer to the **User Manual** included in this repository.

For detailed information regarding the **mathematical model**, **computational implementation**, and **validation** of PRINSAS 2.0, please refer to the accompanying **manuscript**.

## Contact
For fit discussions, bug reports, or feature requests, please contact henry@pnhvu.com.

---
PRINSAS 2.0 is an open-source tool designed to improve accessibility and usability for researchers working with scattering data. Contributions and feedback are welcome!
