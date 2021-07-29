# VSP2WOPWOP

This program interfaces OpenVSP with PSU-WOPWOP by generating all of the necessary 
patch and functional data input files needed to run noise predictions for hovering rotors. 
VSP2WOPWOP is intended to (1) minimize the learning curve faced by new researchers 
using PSU-WOPWOP, (2) provide a means for attaining rapid, yet accurate acoustic predictions, 
(3) provide a framework for conducting parametric studies consisting of various blade geometries and loading conditions. 

This program extracts the geometric properties of a rotor blade and airfoil polar data from OpenVSP's 
DegenGeom analysis and XFoil's output files, respectively. These quantities are then referenced to
compute the aerodynamic blade loads using blade element momentum theory (BEMT). Lastly, the input 
files are formatted and written to a directory of the user's choosing. 

## Getting Started

VSP2WOPWOP is comprised of an input (input.py) and a main executable module (VSP2WOPWOP) along with the core functions. The main executable module coordinates when the
the core functions are executed and handles any directory management tasks. The core functions are located in the functions directory, and a description of all these functions is provided below.
The program can be ran directly from a command line by navigating to the program's directory and executing the command "python VSP2WOPWOP". 
Alternatively, a Python IDE can be used to run the program. This will return several dictionaries, which are accessible to the user
and contain a variety of computed geometric and loading parameters. 

## File Descriptions
- input.py:

The input file is used to set up the cases for which to generate the patch and functional data files. VSP2WOPWOP supports any number of DegenGeom files
and loading conditions, which are specified in the comma-delimited lists in this file (dataFileNames and CT). A folder that corresponds 
to each DegenGeom file will be automatically generated in the user specified directory. This folder will contain the geometry and compact 
geometry files as well as all the generated functional data files, which contain the loading information. This is really the only script that the user needs to modify. 

- VSP2WOPWOP.py:

This is the core script that coordinates when each function is executed and sets up the directories to which the generated files would be written to.
Three dictionaries titled XsecPolar, geomParams, and loadParams are returned after running this module. XsecPolar contains the airfoil polars 
for each airfoil cross section along with several other related quantities such as the maximum lift coefficient and the lift curve slope. The geomParams dictionary contains
a variety calculated geometric properties, such as the chord, twist, and solidity distributions. The loadParams dictionary contains information on the trimmed rotor and blade loads. 
Any entrees in these dictionaries can be accessed by calling the name of the dictionary followed by the key that corresponds to that entree in single quotations and enclosed in brackets (e.g. geomParams[ 'chordDist' ]).


- DegenGeom.py:

This function parses through the DegenGeom .csv file attained from OpenVSP and places the degenerate geometry components into a dictionary.
This function can work independently from the rest of the program and supports any number of arbitrary components, 
which do not necessarily need to be rotor blades. Each component along with its degenerate geometry information is placed in separate dictionary entrees. 

- GeomProcess.py:

This function analyzes the blade geometry to determine various geometric parameters, such as the blade solidity, twist, tapper distributions.
All of these quantities are then assembled into a dictionary. It should be noted that the method used to compute the blade twist is sensitive 
to the orientation of the blade geometry in OpenVSP, please refer to the OpenVSP session in the test case to see the proper orientation. 

- NodeCenteredNorms.py:

This function computes the face-area scaled node-centered normal vectors needed for the geometry patch files.
Even though the DegenGeom file contains the face centered normals, the version of PSU-WOPWOP that I have been using, while developing this program
throws an error when I used the face-centered normal vectors, so I wrote this function to bypass this bug.

- polarRead.py:

This function parses through the XFoil polar files and detects the nominal and maximum lift coefficient along with the corresponding angles of attack. The lift curve slope is also computed based on the
angle of attack interval specified in the input file. The lift coefficient polars and the region over which the lift curve slope is computed are plotted when the "check" variable in the input module is set equal to one. These quantities are then assembled into a dictionary. 

- loading.py:

This function trims the rotor to the desired thrust condition, which is specified in the input file, b adjusting the collective pitch setting.  
The aerodynamic blade loads are then computed using BEMT. The user can choose whether to include Prandtl's tip loss formulation.
The computed quantities are also assembled into a dictionary.

- GeomPatchFileWrite.py:

This function formats and writes a binary patch file for a constant structured geometry. The units, comments, and zoneName variables
are most likely the only parameters that the user may wish to modify.

- CompactGeomPatchFileWrite.py:

This function writes a compact geometry binary patch file upon which the loading vectors are prescribed. The units, comments and zoneName variables
are most likely the only parameters that the user may wish to modify.

- LoadingPatchFileWrite.py:

This function writes a compact constant loading binary functional data file. The comments and zoneName variables
are most likely the only parameters that the user may wish to modify.

## Requirements

In order for the program to run the user must provide the DegenGeom .csv file that is generated by OpenVSP. 
The program only supports uniform blade geometries, thus only a single blade needs to be modeled in OpenVSP. 
The geometry and loading of the single blade is duplicated in the PSU-WOPWOP namelist file to correspond to the actual number of blades.
The orientation of the blade is also important, please refer to the OpenVSP session in the test case to see the proper orientation. 
VSP2WOPWOP also requires at least one airfoil polar output file generated by XFoil. 

## Dependencies

VSP2WOPWOP requires the user to have Python 3 installed as well as the numpy and matplotlib modules. 

## Test Case

A test case of the Boeing Model 360 rotor has been developed and validated against experimental data. The folder titled "TestCase"
contains the corresponding OpenVSP session, DegenGeom file, airfoil polars, and the paper that was used to validate the predicted results. 
The written input files for PSU-WOPWOP are placed in the "PatchFiles" folder and the "WOPWOP_Case_Results" contains all the
necessary files to run the case, including the PSU-WOPWOP executable file. The output files from PSU-WOPWOP as well as the figures 
comparing the predicted and measured pressure time histories are included here.

## Future Implementations

In upcoming versions of VSP2WOPWOP, I intend to incorporate predictive capabilities for the blade loads in axial and forward flight. 
I will also add the option for the user to select whether to use the collective pitch or the rotational rate as the trim variable. Additionally,
I will write a function which generates the necessary patch files for PSU-WOPWOP to run a broadband noise prediction. 

## Authors

* **Daniel Weitsman** - [dfw5266@psu.edu]