0.) Note that the following steps are outlined for a Windows 64-bit installation; similar steps may need to be taken on a Linux machine.


1.) Python and pip installation / updates:
------------------------------------------
- download & install python 3.9 or 3.10.10 64-bit on Windows from here: https://www.python.org/downloads/windows/; I recommend doing a customized installation of 3.10.10 (there pick all options, such as "install for all users"). [Note: there are issues with newer versions like 3.11 and the python bindings of openbabel on Windows; to be revisited in future]

- make sure that the python 'pip' is installed and up to date:
	in a Windows command prompt with administrator rights (click Run as Administrator) type:  
	python -m pip install --upgrade pip
	
	
2.) Openbabel installation:
---------------------------
- download & install openbabel 3.1.1 (GUI) 64-bit version for Windows (download the installer file from https://github.com/openbabel/openbabel/releases/ ); [a newer one may exist]

- IMPORTANT step after: delete BABEL_DATADIR variable in the environment variables of Windows; (type in search box: "edit environment variables" and open the control panel for that.


3.) Python bindings to openbabel (and pybel):
---------------------------------------------
- to install Python openbabel bindings on Windows, I recommend using a precompiled .whl wheel downloaded from https://www.lfd.uci.edu/~gohlke/pythonlibs/#openbabel; (for reasons, see this related post: https://stackoverflow.com/questions/42151301/how-do-i-install-openbabel-for-python-3-6-in-windows-10).
  Choose an amd64 version that is compatible with the current Python version (check for compatible versions in command prompt via:  pip debug --verbose ).
  Then change the directory in the cmd prompt (run as admin) to your download directory and install using a pip command like:  pip install openbabel-3.1.1-cp310-cp310-win_amd64.whl

- [Alternative (to test again in future); try to install python bindings (didn't work in past): run in a command prompt with administrator rights: pip install -U openbabel or:  pip install openbabel==3.1.1  for a specific version; see also: https://pypi.org/project/openbabel/]

- [Alternative (added April 2025): check pip install openbabel-wheel , see also: https://pypi.org/project/openbabel-wheel/]


4.) Indigo toolkit installation
-------------------------------
For use of the Indigo toolkit, the indigo package and python bindings need to be installed:
- run in a command prompt with administrator rights:  
	pip install epam.indigo
	
	
5.) Check for successful installations:
----------------------------------------
in a command prompt type: obabel   [the version of the open Babel software installed should be shown]
open Tools_for_SMILES_conversion program in Visual Studio; 
check that Python 3.8 (64-bit) is set as the Python Environment (or in future a newer version);
check that the above Python environment includes epam.indigo as well as openbabel;
(1) try running the program with a test input file like "InputFiles/smiles_0002.txt" see settings near code line 60;
(2) run Run_SMILES_to_AIOMFAC_input_AND_Vaporpressure_Prediction.bat from the Tools_for_SMILES_conversion folder.
if things run without errors, then all necessary packages should be present -- else, check for missing packages and/or potential causes of errors (such as conflict with another installed Python version; One potential issue with running the .bat files is that your default Python version is different from the one installed above. If so, either run the .bat file in a dedicated Python environment for a compatible version with openbabel installed or uninstall the newer versions of Python (unless you really need those).
