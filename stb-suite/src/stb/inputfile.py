#!/usr/bin/env python3

#################################################
#     Siesta Tool Box - Suite                   #
# Developed by Dr. Carlos M. O. Bastos          #
#      bastoscmo.github.io                      #
#################################################
    
VERSION = "1.8.0"    

from time import sleep
import numpy as np
import math
import sys
import os
import shutil
import argparse

COLORS = {                          
    'reset': '\033[0m',             
    'cyan': '\033[96m',   
    'blue': '\033[94m',   
    'green': '\033[92m',  
    'yellow': '\033[93m', 
    'red': '\033[91m',    
    'bold': '\033[1m',    
    'underline': '\033[4m'
}                         
                          
def color_text(text: str, color: str) -> str:
    """Retorna texto formatado com cor ANSI"""
    return f"{COLORS[color]}{text}{COLORS['reset']}"
                          
def show_intro() -> None: 
    """Exibe a introdução estilizada da STB-SUITE"""
    os.system('cls' if os.name == 'nt' else 'clear')

    logo = color_text(r"""
.----------------.  .----------------.  .----------------.
| .--------------. || .--------------. || .--------------. |
| |    _______   | || |  _________   | || |   ______     | |
| |   /  ___  |  | || | |  _   _  |  | || |  |_   _ \    | |
| |  |  (__ \_|  | || | |_/ | | \_|  | || |    | |_) |   | |
| |   '.___`-.   | || |     | |      | || |    |  __'.   | |
| |  |`\____) |  | || |    _| |_     | || |   _| |__) |  | |
| |  |_______.'  | || |   |_____|    | || |  |_______/   | |
| |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'
 """, 'cyan')

    description = [
        "Siesta ToolBox Suite",
        "A comprehensive toolkit for SIESTA DFT simulations",
        f"Version {VERSION} | University of Brasilia - 2025",
        "Developed by Dr. Carlos M. O. Bastos"
    ]

    print(logo)
    print("\n" + "="*60)
    for line in description:
        print(line.center(60))
        sleep(0.2)
    print("="*60 + "\n")
    return




# --- Internal Relax Template ---
CALC_RELAX_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include structure.fdf

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.01 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH 
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [6  6  6]
Mesh.CutOff           320  Ry  
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K

## Optional Parameters
  # SCF.MustConverge true 
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   false

## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun            CG
MD.VariableCell         true
MD.MaxForceTol          0.01 eV/Ang
MD.Steps                150

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates 
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize 
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.  #(Use previous Density Matrix if available)

## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.  
  # ON.UseSaveLWF  .false.  

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.

## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Relax Template ---


# --- Internal Total Energy Template ---
CALC_TOTAL_ENERGY_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  # %block Supercell
  # %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [4  13  1]
Mesh.CutOff           320  Ry   
FilterCutoff          150  Ry

## Optional Parameters
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K

## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Total Energy Template ---


# --- Internal AIMD Template ---
CALC_AIMD_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

%block Supercell
   2   0   0
   0   2   0
   0   0   2
%endblock Supercell

## Optional Parameters
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       SZ
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry    
FilterCutoff          150  Ry

## Optional Parameters
  
  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-4  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.
## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## STRUCTURE RELAXATION (MOLECULAR DYNAMICS - MD)
## ===================================================================
##Essential Parameters
MD.TypeOfRun  Nose
MD.InitialTimeStep  1
MD.FinalTimeStep  5000
MD.LengthTimeStep  1.0  fs
MD.InitialTemperature  300.0  K
MD.TargetTemperature   300.0  K
MD.NoseMass  30.0  Ry*fs**2

## Optional Parameters
  # MD.RelaxCellOnly true
  # Constant.Volume false
  # MD.MaxStressTol 1 GPa
  # MD.MaxDispl 0.2 Bohr
  # MD.PreconditionVariableCell 5 Ang
  # MD.Broyden.History.Steps 5
  # MD.Broyden.Cycle.On.Mait true
  # Target.Pressure 0 GPa
  # MD.RemoveIntramáscularPressure false

## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
MD.UseSaveXV            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF  .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteMDHistory          .true.
WriteCoorStep           .true.
WriteForces             .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of AIMD Template ---


# --- Internal Band Structure Template ---
CALC_BANDS_TEMPLATE = """
## ===================================================================
## SYSTEM DEFINITION
## ===================================================================
##Essential Parameters
SystemLabel      siesta
SystemName       siesta

## Include the structure file generated with STB-SUITE
%include struct.fdf

## Optional Parameters
  #%block Supercell  %endblock Supercell
  # %block Zmatrix
  # %endblock Zmatrix
  # ZM.UnitsLength Bohr
  # ZM.UnitsAngle  rad

## ===================================================================
## BASIS SET DEFINITION
## ===================================================================
##Essential Parameters
PAO.BasisType       split
PAO.BasisSize       DZP
PAO.EnergyShift     0.02 Ry

## Optional Parameters
  # PAO.SplitNorm 0.15
  # PAO.SplitNormH
  # PAO.SplitTailNorm true
  # PAO.SplitValence.Legacy false
  # PAO.FixSplitTable false
  # PAO.EnergyCutoff 20 Ry
  # PAO.EnergyPolCutoff 20 Ry
  # PAO.ContractionCutoff 0
  # %block PAO.Basis
  # %endblock PAO.Basis

## ===================================================================
## K-POINT SAMPLING (BRILLOUIN ZONE)
## ===================================================================
##Essential Parameters
kgrid.MonkhorstPack   [1  1  1]
Mesh.CutOff           320  Ry 
FilterCutoff          150  Ry

## Optional Parameters

  # %block kgrid.MonkhostPack
  # %endblock MonkhostPack
  # ChangeKgridInMD true
  # TimeReversalSymmetryForKpoint true

## ===================================================================
## EXCHANGE-CORRELATION (XC) FUNCTIONAL
## ===================================================================
##Essential Parameters
XC.Functional       GGA
XC.Authors          PBE

## Optional Parameters
  # %block XC.mix
  # %endblock XC.mix

## ===================================================================
## SPIN POLARIZATION
## ===================================================================
##Essential Parameters
Spin                non-polarized

## Optional Parameters
  # Spin.Fix false
  # Spin.Total 0
  # SingleExcitation False
  # Spin.OrbitStrength 1.0

## ===================================================================
## SELF-CONSISTENT-FIELD (SCF)
## ===================================================================
##Essential Parameters
MaxSCFIterations        300
SCF.Mixer.Weight        0.1
SCF.DM.Tolerance        1.0d-5  eV
SCF.Mixer.History       6
ElectronicTemperature   300 K
Diag.ParallelOverK      .true.
## Optional Parameters
  # SCF.MustConverge true
  # SCF.Mix Hamiltonian
  # SCF.Mix.Spin all
  # SCF.Mixer.Method Pulay
  # SCF.Mixer.Variant original
  # SCF.Mixer.Kick 0
  # SCF.Mixer.Kick.Weight
  # SCF.Mixer.Restart 0
  # SCF.Mixer.Restart.Save 1
  # SCF.Mixer.Linear.After -1
  # SCF.Mixer.Linear.After.Weight
  # SCF.Mixer.History 10

## ===================================================================
## VAN DER WAALS CORRECTIONS (GRIMME DFT-D3)
## ===================================================================
##Essential Parameters
DFTD3                   .false.
## Optional Parameters
  # DFTD3.UseXCDefaults true
  # DFTD3.BJdamping true
  # DFTD3.s6 1.0
  # DFTD3.rs6 1.0
  # DFTD3.s8 1.0
  # DFTD3.rs8 1.0
  # DFTD3.alpha 14.0
  # DFTD3.a1 0.4
  # DFTD3.a2 5.0
  # DFTD3.2bodyCutOff 60.0 bohr
  # DFTD3.BodyCutoff 40.0 bohr
  # DFTD3.CoordinationCutoff 10.0 bohr

## ===================================================================
## BAND STRUCTURE AND DENSITY OF STATES
## ===================================================================
##Essential Parameters
#___
#IMPORTANT: k-path file created by STB-SUITE.
#For manual k-path use optional
#Parameters "%block BandLines" and comment the bellow line.
%include kpath_bs.fdf
#___

%PDOS.kgrid_Monkhorst_Pack
 20     0    0  0.0
  0    20    0  0.0
  0     0   20  0.0
%end PDOS.kgrid_Monkhorst_Pack

%block  ProjectedDensityOfStates
-20.0  20.0  0.01  8000  eV
%endblock  ProjectedDensityOfStates

## Optional Parameters
#BandLinesScale  ReciprocalLatticeVectors
# %block BandLines
#  Important: Include here the k-path
# %endblock BandLines


## ===================================================================
## NUMERICAL SOLUTION (ELECTRONIC STRUCTURE)
## ===================================================================
## Optional Parameters (All parameters in this section were optional)
  # SolutionMethod diagon
  # NumberOfEigenStates
  # Diag.Use2D true
  # Diag.ProcessorY (can be ~sqrt(N processors))
  # Diag.BlockSize
  # Diag.Algorithm Divide-and-Conquer
  # Diag.Memory 1
  # Diag.UpperLower lower
  # OccupationFunction FD
  # OccupationMPOrder 1

## ===================================================================
## STARTING THE CALCULATION (Using previous results)
## ===================================================================
##Essential Parameters
DM.UseSaveDM            .true.
## Optional Parameters (Set to false or inactive)
  # MD.UseSaveXV  .true.
  # MD.UseSaveZM  .false.
  # MD.UseSaveCG  .false.
  # ON.UseSaveLWF .false.

## ===================================================================
## OUTPUT FILES
## ===================================================================
##Essential Parameters (Active output flags)
WriteCoorInitial        .true.
WriteCoorXmol           .true.
WriteMDXmol             .true.
## Optional Parameters (Interface flags)
  #WriteMDHistory          .true.
  #WriteCoorStep           .true.
  #WriteForces             .true.
  #SaveHS                 .true.
  #SaveRho                .true.
  #CDF.Save               .true.
  #CDF.Compress           9

## Optional Parameters (Other output flags)
  # LongOutPut false
  # WriteKpoints false
  # WriteOrbMom false
  # SCF.Write.Extra false
  # WriteEigenvalues false
  # On.UseSaveLWF false

## Optional Parameters (Wave Function Output)
  # WriteKBands false
  # WriteBands false
  # WFS.Write.For.Bands false
  # WFS.Band.Min 1
  # WFS.Band.Max 30
  # WaveFuncKPointsScale pi/a
"""
# --- End of Band Structure Template ---


# --- Calculation Logic Functions ---

def parse_structure_fdf(filename):
    """
    Parses a .fdf (Siesta) file and returns the lattice vectors and
    the list of chemical species symbols.
    Modified to 'raise Exception' instead of 'sys.exit'.
    """
    in_lattice_block = False
    in_species_block = False
    all_values = []
    species_list = []
    lattice_constant = 1.0 # Default value

    try:
        with open(filename, 'r') as f:
            for line in f:
                cleaned_line = line.split('#', 1)[0].strip()
                if not cleaned_line: continue
                parts = cleaned_line.split()
                lower_line = cleaned_line.lower()

                if lower_line.startswith('latticeconstant'):
                    try:
                        lattice_constant = float(parts[1])
                    except (IndexError, ValueError):
                        raise ValueError(f"'LatticeConstant' line malformed: {line}")
                    continue 

                if lower_line == '%block latticevectors':
                    in_lattice_block = True; continue 
                if lower_line == '%endblock latticevectors':
                    in_lattice_block = False; continue 
                if lower_line == '%block chemicalspecieslabel':
                    in_species_block = True; continue
                if lower_line == '%endblock chemicalspecieslabel':
                    in_species_block = False; continue

                if in_lattice_block:
                    for part in parts:
                        try: all_values.append(float(part))
                        except ValueError: pass
                
                if in_species_block:
                    try:
                        symbol = parts[2] # ex: 'C', 'B', 'N'
                        if symbol not in species_list:
                            species_list.append(symbol)
                    except IndexError: pass
        
        if len(all_values) != 9:
            raise ValueError(f"Expected 9 values in LatticeVectors block, found {len(all_values)}.")
            
        lattice = np.array(all_values).reshape(3, 3) * lattice_constant
        
        if not species_list:
            raise ValueError("No chemical species found in ChemicalSpeciesLabel block.")

        return lattice, species_list

    except FileNotFoundError:
        raise FileNotFoundError(f"Structure file '{filename}' not found.")
    except Exception as e:
        # Re-raise other errors
        raise e

def compute_monkhorts(cella, cellb, cellc, k_density):
    """
    Calculates the reciprocal vectors and the number of Monkhorst-Pack divisions.
    Modified to 'raise' errors.
    """
    volume = np.dot(cella, np.cross(cellb, cellc))
    
    if abs(volume) < 1e-9:
        raise ValueError("Cell volume is zero. Check lattice vectors.")
        
    b1 = 2 * np.pi * np.cross(cellb, cellc) / volume
    b2 = 2 * np.pi * np.cross(cellc, cella) / volume
    b3 = 2 * np.pi * np.cross(cella, cellb) / volume

    lengths = [np.linalg.norm(b) for b in (b1, b2, b3)]
    divisions = [max(1, math.ceil(length / k_density)) for length in lengths]
    return divisions

def copy_pseudopotentials(species_list, pp_path):
    """
    Copies the required .psml files from the 'pp_path' folder to
    the current working directory.
    Returns a list of warnings/errors.
    """
    warnings = []
    if not os.path.isdir(pp_path):
        warnings.append(f"Warning: PP path '{pp_path}' is not a valid directory. PPs not copied.")
        return warnings

    print("\n--- Copying Pseudopotentials ---")
    copied_files = 0
    for symbol in species_list:
        pp_filename = f"{symbol}.psml"
        source_file = os.path.join(pp_path, pp_filename)
        dest_file = os.path.join(os.getcwd(), pp_filename)
        
        if os.path.exists(source_file):
            try:
                shutil.copy2(source_file, dest_file)
                print(f"  ✅ Copied: {pp_filename}")
                copied_files += 1
            except Exception as e:
                warnings.append(f"  [ERROR] Failed to copy {pp_filename}: {e}")
        else:
            warnings.append(f"  [WARNING] File not found: {source_file}")

    if copied_files == len(species_list):
        print("All pseudopotentials copied successfully.")
    else:
        warnings.append("Warning: Some PP files were not found or could not be copied.")
    print("----------------------------------")
    return warnings

# --- Main Generation Logic Function ---
def generate_calculation(struct_file, chosen_mode, pp_path):
    """
    Main function that executes the file generation logic.
    """
    try:
        # --- 1. Input Validation ---
        if not os.path.exists(struct_file):
            print(f"Error: Structure file '{struct_file}' not found.")
            return # Abort the function

        if not pp_path:
            print("ℹ️  PP path is blank. Skipping pseudopotential copy.")
        else:
            print(f"ℹ️  PP path set to: {pp_path}")

        print(f"ℹ️  Mode selected: {chosen_mode}")

        # --- 2. Start Generation Logic ---
        print(f"\nProcessing '{struct_file}'...")
        
        output_file = "calc.fdf" # Fixed output name
        
        # --- Select Template ---
        template_string = ""
        if chosen_mode in ['relax', 'relax+d3']:
            template_string = CALC_RELAX_TEMPLATE
            print("ℹ️  Using 'Relax' template.")
        elif chosen_mode in ['total_energy', 'total_energy+d3']:
            template_string = CALC_TOTAL_ENERGY_TEMPLATE
            print("ℹ️  Using 'Total Energy' template.")
        elif chosen_mode in ['aimd', 'aimd+d3']:
            template_string = CALC_AIMD_TEMPLATE
            print("ℹ️  Using 'AIMD' template.")
        elif chosen_mode in ['bands', 'bands+d3']:
            template_string = CALC_BANDS_TEMPLATE
            print("ℹ️  Using 'Bands' template.")
        
        # --- Set D3 flag ---
        d3_flag = ".false."
        if chosen_mode in ['relax+d3', 'total_energy+d3', 'aimd+d3', 'bands+d3']:
            d3_flag = ".true."
            print("ℹ️  DFT-D3 (van der Waals) correction will be ENABLED.")
        else:
            # 'relax', 'total_energy', 'aimd', or 'bands'
            print("ℹ️  DFT-D3 (van der Waals) correction will be DISABLED.")
        
        d3_line_new = f"DFTD3                   {d3_flag}"
        template_lines = template_string.splitlines(keepends=True)

        # Parse structure and get species
        lattice, species = parse_structure_fdf(struct_file)
        print(f"  Species found: {', '.join(species)}")

        # --- K-Grid Conditional Logic ---
        replace_kgrid = False
        kgrid_line_new = "" # Initialize
        
        if chosen_mode in ['aimd', 'aimd+d3']:
            # For AIMD, do nothing, keep the K-grid from the template
            replace_kgrid = False
            print("ℹ️  AIMD mode selected. K-grid will NOT be modified.")
        else:
            # For Total Energy, Relax, and Bands, calculate the K-grid
            replace_kgrid = True
            print("ℹ️  Calculating K-grid (density = 0.2 1/Å)...")
            k_density = 0.2
            kgrid_divs = compute_monkhorts(lattice[0], lattice[1], lattice[2], k_density)
            kgrid_line_new = f"kgrid.MonkhorstPack   [{kgrid_divs[0]}  {kgrid_divs[1]}  {kgrid_divs[2]}]"
            print(f"  Suggested K-grid: {kgrid_divs[0]} {kgrid_divs[1]} {kgrid_divs[2]}")

        
        # Use only the base file name for the include
        include_line_new = f"%include {os.path.basename(struct_file)}"

        # Process the template and replace lines
        print(f"Writing '{output_file}'...")
        output_lines = []
        for line in template_lines:
            line_stripped_lower = line.strip().lower()
            
            # Replace the structure include line
            if line_stripped_lower.startswith('%include') and 'struct' in line_stripped_lower:
                output_lines.append(include_line_new + '\n')
            
            # Replace the k-grid (if not AIMD)
            elif line_stripped_lower.startswith('kgrid.monkhorstpack'):
                if replace_kgrid:
                    output_lines.append(kgrid_line_new + '\n')
                else:
                    # Keep the original template line (for AIMD)
                    output_lines.append(line)
            
            # Replace the D3 flag
            elif line_stripped_lower.startswith('dftd3'):
                output_lines.append(d3_line_new + '\n')
            
            else:
                # Keep the original template line
                output_lines.append(line)
        
        # Write the output file
        with open(output_file, 'w') as f_output:
            f_output.writelines(output_lines)
        
        print(f"\n✅ File '{output_file}' generated successfully.")

        # --- 5. Copy Pseudopotentials ---
        pp_warnings = []
        if pp_path:
            pp_warnings = copy_pseudopotentials(species, pp_path)
        
        if pp_warnings:
            print("\n--- Pseudopotential Log ---")
            for warning in pp_warnings:
                print(warning)

    except Exception as e:
        # Catch any errors (e.g., File not found, Zero cell volume)
        print(f"\n--- ERROR ---")
        print(f"An error occurred: {e}")
        print("Operation aborted.")

# --- Main Function (Argument Parser) ---
def main():
    """
    Processes the command-line arguments and calls the generation logic.
    """
    # --- Definition of valid modes ---
    mode_list = [
        'total_energy', 'total_energy+d3',
        'relax', 'relax+d3',
        'aimd', 'aimd+d3',
        'bands', 'bands+d3'
    ]

    # --- Configure argparse ---
    parser = argparse.ArgumentParser(
        description="Script to prepare a Siesta input file (.fdf).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Usage examples:\n"
               "  %(prog)s struct.fdf -t relax+d3 -p /path/to/pps\n"
               "  %(prog)s system.fdf -t bands\n"
               "  %(prog)s system.fdf --type aimd\n"
    )
    
    # Argument 1: Structure File (Required, Positional)
    parser.add_argument(
        "structure_file",
        type=str,
        help="Path to the structure file (e.g., struct.fdf)"
    )
    
    # Argument 2: Calculation Mode (Required, with Flag)
    parser.add_argument(
        "-t", "--type",
        type=str,
        required=True, # <-- The -t/--type flag is now mandatory
        choices=mode_list,
        dest="calc_type", # The value will be stored in 'args.calc_type'
        help=f"Calculation mode. Valid choices: {', '.join(mode_list)}"
    )
    
    # Argument 3: PPs Path (Optional, with flag)
    parser.add_argument(
        "-p", "--pp-path",
        type=str,
        default="",
        dest="pp_path",
        help="Path to the pseudopotentials folder - Only psml format (optional)"
    )

    
    parser.add_argument("-v", "--version", action="version",
                        version=f"stb-inputfile {VERSION}")
    parser.add_argument("--no-intro", dest="intro", action="store_false", help="Do not show the introduction")

    args = parser.parse_args()


    if args.intro == True:
        show_intro()

    print("\n" + color_text("Create a input file from structure in fdf format (converted by STB):", 'bold'))
    print("-"*60)


    # Process the arguments provided on the command line
    args = parser.parse_args()

    # Call the main logic function with the processed arguments
    print("--- Siesta Calculation Preparer ---")
    generate_calculation(args.structure_file, args.calc_type, args.pp_path)


# --- Script Entry Point ---
if __name__ == "__main__":
    main()
