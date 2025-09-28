# S2AS tool – convert SMILES to AIOMFAC input files
The **S2AS tool** generates a valid AIOMFAC model input file for any system characterized by a list of molecular organic components in terms of their Simplified Molecular Input Line Entry System ([SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System "SMILES")) descriptors. To do so, the S2AS code relies on a customized list of SMILES arbitrary target specification ([SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification "SMARTS")) patterns to detect and match individual AIOMFAC subgroups and to treat exceptions consistently. The general motivation and details about the S2AS tool and its applications can be found in an associated publication by Amaladhasan et al., (*in prep.*). The abbreviation AIOMFAC stands for Aerosol Inorganic–Organic Mixtures Functional groups Activity Coefficients (see also [AIOMFAC website](https://aiomfac.lab.mcgill.ca "AIOMFAC") and the [AIOMFAC code repository](https://github.com/andizuend/AIOMFAC)).

----
### Contents of this file
- [Dependencies](#dependencies)
- [Installation instructions](#installation-instructions)
- [Test the installation](#5-test-the-installation)
- [Quick guide to running S2AS](#quick-guide-to-running-s2as)
- [Use as input for AIOMFAC-web](#use-as-input-for-aiomfac-web)
----

## Dependencies
- Open Babel and its Python bindings (see installation instructions below)
- epam Indigo package
- (optional) [MS Visual Studio Community](https://visualstudio.microsoft.com/vs/community/) (for use of included project and solution files `SMILES_to_AIOMFAC.sln`)

## Installation instructions
> [!NOTE] 
> The following steps are first outlined for a Windows 64-bit installation (denoted by steps tagged as [Windows]). Equivalent steps are also shown for installation on a Linux machine (denoted by tag [Linux]). The Linux steps were tested with Ubuntu 24.04 LTS; the details for other Linux distributions may differ. Depending on your operating system, only follow either the [Windows] or the [Linux] steps (not both).

#### (1) Python and pip installation / updates
- [Windows] Download & install Python v3.13.7 or newer for 64-bit on Windows, e.g. from [here](https://www.python.org/downloads/windows/). Version 1.0 of the S2AS tool has been confirmed to work with Python v3.13.7 (as well as the older v3.10.10).
- [Windows] Make sure that the python 'pip' is installed and up to date. In a Windows command prompt with administrator rights (click Run as Administrator) type  
	`python -m pip install --upgrade pip`
- [Linux] Make sure a recent version of python3 is installed (check in a terminal with command `python --version` or `python3 --version`). For convenience make sure that python refers to python3 via command `sudo apt install python-is-python3`
- [Linux] Open a command terminal, navigate to the project folder (e.g. S2AS__SMILES_to_AIOMFAC) and generate a virtual environment (venv) for the project.
  - execute `sudo apt install python3-venv` (if not already installed on your system)
  - create a new venv with command `python3 -m venv .venv`
  - activate your new venv so Python is using it; command: `source .venv/bin/activate`

#### (2) Open Babel
- [Windows] Download & install the Open Babel v3.1.1 (GUI) 64-bit for Windows. Download the [executable installer file](https://github.com/openbabel/openbabel/releases/); a newer version may be available at your time of installation (you could try the installation with that newer version).

#### (3) Python bindings to openbabel (and pybel)
- [Windows] To install Python openbabel bindings on Windows, I recommend using a precompiled wheel (`.whl` file).\
  Try: `pip install openbabel-wheel`, see also: [https://pypi.org/project/openbabel-wheel/](https://pypi.org/project/openbabel-wheel/) for dependencies and compatible Python versions.
- [Linux] To install Python openbabel execute command `pip install openbabel-wheel`

#### (4) Indigo toolkit
[Windows] and [Linux] For use of the Indigo toolkit within the S2AS program, the related package and python bindings need to be installed:
- Run in a command prompt with administrator rights:  `pip install epam.indigo`

#### (5) Test the installation
- [Windows] and [Linux] In a command prompt type: `obabel`\
  The version of the Open Babel software installed should be displayed – if not, check the installation steps above again and/or use an older, compatible version of Open Babel and the Python bindings (successful tests involved Open Babel v3.1.1 GUI 64bit on Windows and openbabel-wheel v3.1.1.22).
- Test with a .bat script file (on Windows) or equivalent Python command:
	- [Windows] Included in the main folder of this repository is a `Run_SMILES_to_AIOMFAC_input.bat` script, which can be edited by any text editor. It executes the S2AS python program for a particular input file case. [Windows] Try running the .bat file (double click to run). If the script executes successfully, the installation above was successful.
	- [Windows] One potential issue with running the .bat files is that your default Python version is different from the one installed above. If so, either run the .bat file in a dedicated Python environment for an installed version of Python compatible with Open Babel.
	- [Linux] For an equivalent test on Linux (or in a [Windows] command prompt): navigate in the command terminal to the `S2AS__SMILES_to_AIOMFAC` folder and execute the following command: \
 ```python SMILES_to_AIOMFAC_input.py ./InputFiles/smiles_1409.txt```.\
 As with the .bat example, this runs the S2AS program for the smiles_1409.txt input file located in folder `InputFiles`. Related output files end up in folder `OutputFiles`.

- [Windows] Optional check when using the included Microsoft Visual Studio Community solution and project files (SMILES_to_AIOMFAC.sln); this requires that you have MS Visual Studio installed with Python support.
	- Open the `SMILES_to_AIOMFAC.sln` program from the main folder (S2AS__SMILES_to_AIOMFAC) in MS Visual Studio. Check that Python 3.xx (64-bit) is set as the Python Environment (the installed version compatible with the Open Babel bindings outlined in step 3). Also check that the epam.indigo package is part of the Python environment set to be used for the current Visual Studio project (if not, add it via the environment-specific package manager in Visual Studio).
	- Run the program with a test input file like "InputFiles/smiles_0001.txt" see related settings near source code line 73 in file `SMILES_to_AIOMFAC_input.py`;
	- If the S2AS program runs without errors in Visual Studio, then all necessary packages should be present; else, check for missing packages and/or potential causes of errors (such as conflict with another installed Python version or environment).

----
# Quick guide to running S2AS
If all you wish to do is to run the S2AS program for your own system of components, this is relatively straightforward task. The following inputs (files) need to be provided.

#### SMILES input file:
1. For your system, i.e. the list of organic components, generate a SMILES input file containing one SMILES string per line, as shown in the example files in folder `InputFiles`. Have a look at the `smiles_0001.txt` example file for its structure and use it as a template.
2. Name and save your file, e.g. as `smiles_0002.txt` or similar, the 4-digit number part can vary, but the prefix `smiles_`is needed.
3. Copy this file into the `InputFiles` folder. Do not list the SMILES of water (or other inorganic components) in that file; S2AS v1.0 only supports the parsing of organic molecules (see details in the article by Amaladhasan et al.).

#### Concentrations input file:
4. Provide a text file containing the concentration information for each component (i.e. each line in the SMILES file). Use the file `input_concentrations_0001.txt` from folder `InputFiles` as a template.
5. Note the very first concentration value stated should be that of water, since the S2AS-generated AIOMFAC input file assumes that water is part of the mixture. Then, the second concentration value will refer to the first organic component (corresponding to the first line in the `smiles_????.txt` file).
6. As stated in the concentration file header, on line 3, you can state the unit of concentration used as either [mol/m^3] or [molec/cm^3]. Only these two options are supported (convert from other units when necessary).
7. Save your `input_concentrations_????.txt` file using the same 4-digit number (in place of ????) as for the `smiles_????.txt` file.
> [!NOTE]
> The concentration values are not needed by the S2AS tool; however, in the case that an invalid or unwanted/unsupported SMILES code is detected (e.g. [O-][O+]=O for ozone or [O]O for hydroperoxyl from a mechanism output), the offending SMILES (row) will be removed and the SMILES list updated. In that case, the corresponding row in the `input_concentrations_????.txt` file will be removed too. Cleaned-up versions of the `smiles_????.txt` file and, if necessary, the `input_concentrations_????.txt` file, will be copied into the `InputFiles` folder as part of the S2AS program execution. For reference, the original SMILES list will also be saved there, renamed as `smiles_????_orig.txt`.

## Running your customized case
Once the two input files have been prepared, you are ready to run the Python script with your case as command line argument or to modify the file name in the .bat script (see also [Test the installation](#5-test-the-installation) about that).
- Option 1: (on Windows)
	- use a text editor like Notepad++ to open and edit the file `Run_SMILES_to_AIOMFAC_input.bat`. In that file, change the file number on the script line `SET _fileN=smiles_1409.txt` to that of your case, e.g. `0002`.
	- Save the file and close the text editor.
 	- Double-click the `Run_SMILES_to_AIOMFAC_input.bat` file to execute this script in a command terminal window. The `SMILES_to_AIOMFAC_input.py` will be executed with your input files and, if everything was set properly, you should see a message in the terminal confirming the processing of your list of SMILES.
  - The generated output will be found as a file labelled `input_????.txt` in folder `OutputFiles`. It is named with prefix `input_` since it is an input file for AIOMFAC-web. You can view that file in a text editor, it states the AIOMFAC subgroups describing each input component. It also includes the AIOMFAC notation for water as component 1.
- Option 2:
	- Open a command terminal and navigate to your local `S2AS__SMILES_to_AIOMFAC` repository directory.
 	- Execute the following command: `python ./SMILES_to_AIOMFAC_input.py ./InputFiles/smiles_????.txt`, where the command line argument 1 is the relative path to and  actual name of your `smiles_????.txt` file. As with Option 1, successful execution of the S2AS program will result in the generation of a file labelled `input_????.txt` in folder `OutputFiles`.
- Option 3: (Visual Studio or other IDE)
	- Alternatively, if code debugging is of interest and/or structures of the molecules and the mapping of AIOMFAC subgroups should be visualized, you can run the Visual Studio solution `SMILES_to_AIOMFAC.sln` (see also [Test the installation](#5-test-the-installation)).
 	- In this case, the file name of your customized SMILES input file should be entered on source code line 73 of `SMILES_to_AIOMFAC_input.py` (within Visual Studio).
  	- Visualizations and additional debbuging output can be enabled on line 51 of the `SMILES_to_AIOMFAC_input.py` file; setting `debugging_verbose = True`.
  	- In the latter case, additional `.png` files will be generated as part of the outputs in folder `OutputFiles`. These images show the structure of a molecule and the step-by-step matching of specific AIOMFAC subgroups that were detected during the parsing of a molecule.
  
### Use as input for AIOMFAC-web
> [!TIP]
> For smaller systems (less than ~50 components), the generated AIOMFAC-web-style input files (`input_????.txt`) can, e.g., be used for mixture activity coefficient computations with the [AIOMFAC-web model](https://aiomfac.lab.mcgill.ca/index.html). For that purpose, the `input_????.txt` should be further modified in a text editor to add, near the end, rows of temperature and component mole fractions (or mass fractions) of interest. The S2AS program writes a single, generic line of composition and temperature inputs for AIOMFAC-web – modify that line and add similar rows (as needed). An example file for such mixture composition inputs is shown [here](https://aiomfac.lab.mcgill.ca/examples/input_examp1.txt). Once your `input_????.txt` file has been edited to include the desired composition points, it can be [uploaded to AIOMFAC-web](https://aiomfac.lab.mcgill.ca/file_upload.html) and run online for generating AIOMFAC predictions of nonideal mixing.
>
> The AIOMFAC thermodynamic mixing model for single-liquid-phase systems and related code are also accessible via [the AIOMFAC repository](https://github.com/andizuend/AIOMFAC). That version operates with S2AS-generated input files.
>
> The generated AIOMFAC-web input file and the related `input_concentrations_????.txt` file can be used as inputs for the [2D polarity–volatility space and lumping framework](https://github.com/andizuend/2D_Polarity_Volatility_lumping). In fact, that repository relies in part on the S2AS tool.
