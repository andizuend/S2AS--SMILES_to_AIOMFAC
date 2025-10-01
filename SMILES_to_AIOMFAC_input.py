########################################################################### 
# 
#  'SMILES_to_AIOMFAC_input.py'  (S2AS tool)
#
#  :: Purpose ::
#  Python script to process a list of organic molecular structures in the 
#  form of SMILES read from an input file to determine the subgroup 
#  description of the molecules for use as input with the AIOMFAC model.
#  The determined subgroup fingerprints of the molecules are stored in array
#  cpsubs and an AIOMFAC-web style input file is written based on this.
#  By default, the AIOMFAC input file will contain also water (see param.).
#
#  :: Requirements ::
#  Aside from a Python installation, the epam indigo toolkit needs to be 
#  installed, e.g. via command line:  pip install epam.indigo
#
#  :: Specific Settings ::
#  Some specific parameters, e.g. for the treatment of radical species, can 
#  be set near the top of this file in the Switches & Parameters section.
#
#  :: Authors & Copyright ::                                         
#  Dalrin Ampritta Amaladhasan, Andreas Zuend, Dan Hassan-Barthaux 
#  Dept. Atmospheric and Oceanic Sciences, McGill University      
#  -> latest changes: 2025-10-01   
#                                                                                   
#  :: License ::                                                                    
#  This program is free software: you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software        
#  Foundation, either version 3 of the License, or (at your option) any later       
#  version.                                                                         
#  The AIOMFAC model code is distributed in the hope that it will be useful, but    
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or    
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more    
#  details.                                                                         
#  You should have received a copy of the GNU General Public License along with this
#  program. If not, see <http://www.gnu.org/licenses/>.                                                        
#
###########################################################################

import os.path
import sys
import time
from shutil import copy
try:
   from indigo import *
except ModuleNotFoundError:
   print("** ERROR **: The 'indigo package' was not imported; it is required to run this program.")
   print("             Please make sure you loaded the correct (virtual) Python environment.")
   print("             See the README file and follow the instructions there.\n")
   time.sleep(2)  # sleep for 5 seconds
   sys.exit(11)   # exit the script with a non-zero status code
#from indigo import *
from indigo import renderer as renderer
import ModSmilesTools
from SMARTS_query_list import SMARTS_AIOMFAC
from S2AS_mapping import detAIOMFACsubgs

#--- Switches & Parameters ----------------------------------------------------
debugging_verbose = False  # set to True to write detailed information to the terminal window and creation of images; set False to avoid most output;
replaceRadicals = True     # if set True, component SMILES indicating radical atoms (e.g. -C[O] or -C[O.]) will be 
                           # modified into non-radical species (by adding implicit H to reach non-radical valence);
replaceRdbOO = True        # if set True, replace SMILES containing C=[O+][O-] with a similar one with C-O-OH hydroperoxy acid group;
waterAsComp01 = True       # add water to the generated AIOMFAC input file (while H2O is not part of the SMILES list);
image_file_type = "png"    # use png, svg or pdf as desired for structure images in (debugging) output;
#------------------------------------------------------------------------------

#initialize some global variables:
timeStart = time.perf_counter()
indigo = Indigo()
renderer = renderer.IndigoRenderer(indigo)
SMILESlist = []
gpfile2 = []

print(' ')
cmdline = sys.argv[0:]
if len(cmdline) > 1:
    smiles_file = sys.argv[1]
    debugging_verbose = False
else: 
    # for debugging, use the following input file path instead of the command line input
    smiles_file = "InputFiles/smiles_0001.txt"
    #smiles_file = "InputFiles/smiles_1409.txt"

if debugging_verbose:
   debugRenderedStructures = True   #set to True for debugging and showing of highlighted subgroup in a molecule (file will be overwritten);
   showRenderedStructures = True    #set to True for debugging and/or graphical output of png images of the structures using the indigo renderer;
else:
   debugRenderedStructures = False
   showRenderedStructures = False 

if debugging_verbose:
   print('The relative path and name of the SMILES input file is : ', str(smiles_file))
   print('')

# determine the relative path to output folder using the location of this  *.py file as the starting point.
# use of the os.path functions is necessary to ensure proper path strings when calling this .py file 
# from a non-local directory, e.g. the 2D lumping framework program
pyfilepath = os.path.abspath(cmdline[0])
locpath = os.path.dirname(pyfilepath)
outpath = os.path.abspath(locpath +'/OutputFiles')

#clean up output figure folder to avoid presence of misleading older files:
dirfilelist = os.listdir(outpath)
for structfile in dirfilelist:
   if structfile.endswith('.png') and structfile.index('Structure_') != None:
      os.remove(os.path.join(outpath, structfile))
   elif structfile.endswith('.pdf') and structfile.index('Structure_') != None:
      os.remove(os.path.join(outpath, structfile))
   elif structfile.endswith('.svg') and structfile.index('Structure_') != None:
      os.remove(os.path.join(outpath, structfile))

AbsPathFName = os.path.abspath(smiles_file)
with open(AbsPathFName, "r") as inpfile:
   for line in inpfile:
      SMILESlist.append(line)
 
#initialize lists:
ii = len(SMILESlist)
is_invalidSMILES = [True]*ii
is_PureAlcohol = [False]*ii
istatus = [-1]*ii
n_except = [0]*ii
highLAtomsList = []        #list of highlighted atoms
igatoms_index = []         #list of ignored atoms in a structure (used in the indigo substructure matching below) 
neigh_atom_list = [] 
cpsubs = [[0 for col in range(0,173)] for row in range(len(SMILESlist))]   #2D array with as many rows as rows in SMILESlist and 172 columns (for the 172 different AIOMFAC subgroups) 


# +++ first main task: processing list of SMILES (i.e. SMILES cleaning) +++++++++++++++++++++++
Bad = []
for ind,molecule in enumerate(SMILESlist):
   st = str(molecule)
   # (1) run string processing on input SMILES to remove unneccessary attributes like / or \ or revise radical atoms; 
   # they are unresolved for AIOMFAC description of structures.
   st = st.replace('/','')
   st = st.replace('\\','')
   if replaceRadicals:           #replacing only radical atoms, but not charged atoms;
      st = st.replace('[O]','O')
      st = st.replace('[O.]','O')
      st = st.replace('[C]','C')
      st = st.replace('[C.]','C')
      st = st.replace('[N]','N')
      st = st.replace('[N.]','N')
      st = st.replace('[S]','S')
      st = st.replace('[S.]','S')
   if replaceRdbOO:
      st = st.replace('=[O+][O-]','OO')
      st = st.replace('[O-][O+]=','OO') 

   SMILESlist[ind] = st.strip()    #strip() removes leading and trailing whitespace, tabs, etc.
   molec = SMILESlist[ind]

   # (2) use special smarts to check validity of SMILES and whether a molecule contains 
   #     a hydroxyl group and no other/different non-carbon functionalities; 
   validS, errortext = ModSmilesTools.verifySMILES(molec)
   is_invalidSMILES[ind] = not validS
   is_PureAlcohol[ind] = ModSmilesTools.PureAlcoholCheck(molec)

   is_PureAlc = is_PureAlcohol[ind]

   info = -77
   if (is_invalidSMILES[ind]):
      print('')
      print('input SMILES is: ' + indigo.loadMolecule(molec).canonicalSmiles() )

      print('WARNING: issues found with SMILES: ' + errortext)    #for debugging
      print('continue anyways... or manually stop execution and correct input file') 
      input()

   elif len(SMILESlist[ind]) < 1:
      is_invalidSMILES[ind] = True

   else:
      #(3) call function for main task of matching AIOMFAC subgroups via SMARTS:
      cpsubs[ind][:], istatus[ind], n_except[ind] = detAIOMFACsubgs(ind, molec, 
                      is_PureAlc, debugging_verbose, debugRenderedStructures, 
                      image_file_type, showRenderedStructures, igatoms_index, 
                      neigh_atom_list, highLAtomsList)

   if istatus[ind] > 0:
      Bad.append(ind)

   if ind > 1 and ind % 1000 == 0:  #modulo test passed;
      print(f'status update: processed SMILES at index: {ind}')

    
# RealBad = [SMILESlist[i] for i in Bad]

# clean up output array in case of gaps in the SMILES list:
if any(is_invalidSMILES):
   ind = is_invalidSMILES.index(True)
else:
   ind = -1

while ind > -1:
   del is_invalidSMILES[ind]
   del SMILESlist[ind]           #remove item/row in list
   del cpsubs[ind]
   del is_PureAlcohol[ind]
   del istatus[ind]
   del n_except[ind]
   if any(is_invalidSMILES):
      ind = is_invalidSMILES.index(True)
   else:
      ind = -1
#+++++ end of first main task +++++++++++++++++++++++++++++++++++++++++++


# save original SMILES input file with extension '_orig.txt' if it does not exist yet 
# and the modified SMILES list as a file with the same name as the input file which contained the SMILES:
AbsPathFNameOrig = AbsPathFName.replace('.txt','_orig.txt')

if os.path.isfile(AbsPathFNameOrig) != True:
   copy(AbsPathFName, AbsPathFNameOrig)

fileA = open(AbsPathFName, 'w')
for item in SMILESlist:
   fileA.write("{}\n".format(item))
fileA.close()

# generate an 'AIOMFAC input file' name in folder ./OutputFiles:
head, smilesfilename1 = os.path.split(smiles_file)
if smilesfilename1[0:7] == 'smiles_':
   AIOMFACfname = smilesfilename1.replace('smiles','input')
else:
   AIOMFACfname = smilesfilename1.replace('_SMILES','')

# use function writeAIOMFACfile with above data:
_ = ModSmilesTools.writeAIOMFACfile(AIOMFACfname, outpath, SMILESlist, cpsubs, waterAsComp01)

timeStop = time.perf_counter()

# explicitly clean up indigo resources before exit
del renderer
del indigo

print('')
print('SMILES processing for AIOMFAC done. The generated AIOMFAC input file,', AIOMFACfname, ', is located')
print(f'in directory: {outpath} \n')

print(f'Number of SMILES processed: {len(SMILESlist)}\n')

print(f'Average number of SMARTS special mappings and exception cases encountered: {sum(n_except) / len(SMILESlist):0.3f} exceptions per SMILES.\n')

print(f'The SMILES processing and file creation took a time of ~ {timeStop - timeStart:0.1f} sec.\n')