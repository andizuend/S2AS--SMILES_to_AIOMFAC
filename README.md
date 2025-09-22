# S2AS tool – convert SMILES to AIOMFAC input files

The **S2AS tool** generates a valid AIOMFAC (-web) model input file for any system characterized by a list of molecular organic components in terms of their Simplified Molecular Input Line Entry System ([SMILES](https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System "SMILES")) descriptors. To do so, the S2AS code relies on a customized list of SMILES arbitrary target specification ([SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification "SMARTS")) patterns to detect and match individual AIOMFAC subgroups and to treat exceptions consistently. The general motivation and details about the S2AS tool and its applications can be found in an associated publication by Amaladhasan et al., (*in prep.*). The abbreviation AIOMFAC stands for Aerosol Inorganic–Organic Mixtures Functional groups Activity Coefficients (see also [AIOMFAC website](https://aiomfac.lab.mcgill.ca "AIOMFAC")).

## Dependencies
- epam Indigo package
