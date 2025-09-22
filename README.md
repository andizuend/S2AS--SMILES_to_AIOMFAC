# S2AS tool – convert SMILES to AIOMFAC input files
The **S2AS tool** generates a valid AIOMFAC (-web) model input file for any system characterized by a list of molecular organic components in terms of their Simplified Molecular Input Line Entry System ([SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System "SMILES")) descriptors. To do so, the S2AS code relies on a customized list of SMILES arbitrary target specification ([SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification "SMARTS")) patterns to detect and match individual AIOMFAC subgroups and to treat exceptions consistently. The general motivation and details about the S2AS tool and its applications can be found in an associated publication by Amaladhasan et al., (*in prep.*). The abbreviation AIOMFAC stands for Aerosol Inorganic–Organic Mixtures Functional groups Activity Coefficients (see also [AIOMFAC website](https://aiomfac.lab.mcgill.ca "AIOMFAC") and the [AIOMFAC code repository](https://github.com/andizuend/AIOMFAC)).

## Dependencies
- Open Babel and its Python bindings (see installation instructions below)
- epam Indigo package

## Installation instructions (updated as of September 2025)
Note that the following steps are outlined for a Windows 64-bit installation; similar steps may need to be taken on a Linux machine, but details will differ.
#### (1) Python and pip installation / updates
- Download & install Python 3.10.10 64-bit on Windows, e.g. from [here](https://www.python.org/downloads/windows/). 
I recommend doing a customized installation of 3.10.10 (there pick all options, such as "install for all users"). Note: there have been issues with certain newer versions like python 3.11 and the python bindings of openbabel on Windows (to be revisited in future).
- Make sure that the python 'pip' is installed and up to date. In a Windows command prompt with administrator rights (click Run as Administrator) type:  
	`python -m pip install --upgrade pip`

#### (2) Open Babel installation
- Download & install the Open Babel v3.1.1 (GUI) 64-bit for Windows. Download the [executable installer file](https://github.com/openbabel/openbabel/releases/); a newer version may be available at your time of installation (you could try the installation with that newer version).

- Important follow-up step: delete the `BABEL_DATADIR` variable in the environment variables of Windows. Type in a search field: `edit environment variables` and open the app for editing such variables.

#### (3) Python bindings to openbabel (and pybel)
- To install Python openbabel bindings on Windows, I recommend using a precompiled wheel (`.whl` file).\
  Try: `pip install openbabel-wheel`, see also: [https://pypi.org/project/openbabel-wheel/](https://pypi.org/project/openbabel-wheel/) for dependencies and compatible Python versions.

- Alternative (to tested on Linux; try to install python bindings via `pip install -U openbabel` or:  `pip install openbabel==3.1.1`  for a specific version; see also [https://pypi.org/project/openbabel/](https://pypi.org/project/openbabel/)
