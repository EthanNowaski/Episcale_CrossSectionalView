DOI: https://zenodo.org/badge/latestdoi/468540858

CUDA code for developmental biology using Subcellular Element Method

Hardware requirement: 
Nvidia video card that supports SM 2.0+ and CUDA 11.0 

Software environment requirement: 
CMAKE ----- Build system.
CUDA  ----- Provide runtime support for parallel GPU computation.
CGAL  ----- Computational geometry library.
Thrust ---- Build-in library of cuda, similar to STL of C++
Paraview -- (Optional) Visualization software for animation purpose. 

*******************************************
Location of configuration files:
 ./resources
 
Additional parameter-configurating files:
  (1) disc_NXX_X.cfg, X is some number, contains:
        animation (vtk) files storage designation.
        animation files name desgination.
        initial contractile spring coefficient.
  (2) disc_M.cfg
        initial parameter for membrane stiffness, bending stiffness, volume exclusion, simulation time, time step size, etc.
  (3) ECM_input.cfg
        ECM stiffness, adhesion with membrane node, volume exclusion with membrane node.
  (4) discMain_M.cpp
        growth speed, anisotropic contractile spring coefficients, etc.
  (5) SceCells.cu
        cell cycle length
        
*******************************************
How to identify the flow of algorithm calling and simulation proceeding?

1. Look for : int main(int argc, char* argv[])
2. Find : simuDomain.runAllLogic_M
3. Go to SimulationDomainGPU.cu and SimulationDomainGPU.h in /src/SceGPU subfolder 
4. Find : SimulationDomainGPU::runAllLogic_M(double & dt...)
5. You will see that in this function, the followings key algorithms are executed:
   (a) nodes.sceForcesDisc_M
     (i) Go to SceNodes.cu and SceNodes.h in /src/SceGPU subfolder for details
   (b)	cells.runAllCellLogicsDisc_M
     (i) Go to SceCells.cu and SceCells.h in /src/SceGPU subfolder for details


For (a):
  This function deals with: (a) Identifying node-to-node correspondence for contractile springs and cell-cell adhesion, (b) Volume exclusion (via Morse-like interaction) between membrane nodes.
  In SceNodes::sceForcesDisc_M, the following key algorithms are executed:
    (I) prepareSceForceComputation_M
    (II) Setting up the coefficients of the apical and basal contractile springs.
    (III) applySceForcesDisc_M
      (i) Identify membrane nodes connected by cell-cell adhesion forces.
      (ii) Identify the membrane nodes connected by contractile springs.
      (iii) Calculate forces due to membrane-to-membrane volume exclusion (or repulsion) using the function : AddForceDisc_M
    (IV) processMembrAdh_M
      (i) Calculate adhesion between membrane nodes using : applyMembrAdh_M

  For (b):
     (I) Assignment for the number of contractile spring for each cell.
     (II) Speed counter for growth progress and etc are set.
     In cells.runAllCellLogicsDisc_M, the following key algorithms are executed:
         (III) Populating nucleus node for dividing cells.
         (IV) computeApicalLoc
         (V) computeBasalLoc
         (VI) computeIndividualCellHeight
         (VII) eCMCellInteraction
         (VIII) computeCenterPos_M2: computes centers of cells.
         (IX) computeInternalAvgPos_M(): computes average nucleus position of each cell.
         (X) applySceCellDisc_M: volume exclusion between membrane and internal (nucleus) nodes.
         (XI) applyMembContraction2: calculate contractile spring force.
         (XII) applyMemForce_M
         (XIII) applyVolumeConstraint
         (XIV) growAtRandom_M_Ver2: advance the growth progress counter with specific conditions related to neighboring 
               cells.
         (XV) allComponentsMove: update position of each node based on previously calculated forces.
 

*******************************************
To run simulation on UCR HPCC cluster: 
   After uploading all the files under this repository, follow the command below.
   (1) cd [The folder where all your files are placed in]
   (2) module load singluarity
   (3) module load centos
   (4) centos.sh
   (5) module load extra; module load GCC; module load cuda;
   (6) (IF it is the first time compiling) module load CGAL; module load cmake;
   (7) (IF it is the first time compiling) cmake .
   (8) make
   (9) After compilation, enter: exit
   (10) sbatch -p gpu --gres=gpu:1 --time=X:00:00 EpiScale_run.sh; (X here is the number of hours you want to keep the simulation running)

*******************************************
Initial conditions:
  The initial conditions can be found in the resources folder. 
  All of the files with the following key words are associated with either the flat or curved tissue shape.
    (1) N01G00 & N01_0 : for flat wing disc
    (2) N02G00 & N02_0 : for curved wing disc

*******************************************
Examples of different perturbations or mechanism that can be tested using this codecan be found in the 
Examples_of_perturbations folder:
  Inside this folder, there is a list of subfolders that contain different examples which will be listed below.
  To run any of the provided examples, the user must replace the files inside any of the following subfolders with the 
  files inside ./src/srcGPU.
    (1) CombinedPerturbation_CellProlifANDCytoskeletalRegulators: 
            * The tissue is split into 2 domains Anterior (A)|Posterior (P). 
            * Cell proliferation and cytoskeletal regulators are perturned only on the posterior side of the tissue. 
                  (i) Case I: Increase proliferation and actomyosin contractility on the posterior side.
                  (ii) Case II: Increase proliferation and ECM stiffness on the posterior side.
                  (iii) Case III: Increase proliferation and cell-ECM adhesion on the posterior side.
    (2) IncreaseORDecrease_ProlifOnPosteriorTissueCompartment:
            * The tissue is split into 2 domains Anterior (A)|Posterior (P). 
            * In this example one can increase or decrease the proliferation rate in the posterior compartment. 
    (3) PatternedProliferationRate
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents the 
              two lateral sides.
            * In this example one can increase or decrease proliferation in the medial or lateral ends to test the 
              role of patterned proliferation.
    (4) Patterning_ActomyosinContractility
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents the 
              two lateral sides.
            * Apical and basal contractility is higher in the medial region. However, the code can be modified to 
              increase or decrease apical and basal contractility in the M or L domains.
    (5) Patterning_BasalECMStiffness
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents the 
              two lateral sides.
            * ECM stiffness is decreased by 50% in lateral sides as compared to the medial domain. 
            * The code can be easily modified to increase or decrease ECM stiffness in any of the domains.
    (6) Patterning_cellECMAdhesion
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents 
              the two lateral sides.
            * Cell-ECM adhesion is increased by 50% in the lateral sides as compared to the medial domain. 
            * The code can be easily modified to increase or decrease cell-ECM adhesion in any of the domains.
    (7) Patterning_cellMembraneTension
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents the 
              two lateral sides.
            * The code can be easily modified to increase or decrease apical, basal and lateral membrane tension 
              in any of the 3 domains.
    (8) RoleofCellProlif_and_patternedCellECMAdhesion
            * The tissue is split into 3 domains L|M|L where M represents the medial domain and L represents the 
              two lateral sides.
            * In this example we test the case of uniform proliferation across the tissue but higher cell-ECM adhesion 
              in the medial domain.
            * The code can be modified to pattern proliferation as well as increase or decrease cell-ECM adhesion in 
              any of the domains.
    (9) VaryCellPressure_by_ChangingCellVolume
            * In this example, cell pressure is varied by changing the target cell volume in the Lagrange Multiplier.
            * This can be tested using (1) a linear gradient or (2) a step function. Both examples have been provided.

********************************************
Output files generated by the model simulations:
    (1) The print out file that is generated for every job submission. This file will be named according to what the user 
        decides to label the output in line 7 of EpiScale_run.sh. 
            * For example, if the user sets output=Test, then the print out file will be called Test. 
            * This file includes important information such as number of cells, cell order, number of ECM nodes, etc.
    (2) dataOutput folder: Stores all txt files with a "detailedStat_'' prefix. These files include information about every 
        cell such as its rank, number of neighbors, cell area, cell center, cell pressure and the position of the apical 
        node, basal node and nucleus center. 
            * For an example, see the provided sample file detailedStat_example.txt.
    (3) ECMFolder: Includes txt files with information about the ECM nodes along with the vtk files for visualizing the 
        ECM nodes.
    (4) animation folder: Includes vtk files for viewing simulation results. These files store the information for 
        visualizing the nuclei and membrane nodes.
            * machine1 folder: if the flat tissue initial condition is used, the vtk files will be stored here. 
            * machine2 folder: if the curved tissue initial condition is used, the vtk files will be stored here.

Note: The vtk files inside the ECMFolder and animation folders are required for visualizing the results and creating simulation videos.
   
********************************************
[Data processing]
The data processing files can be found inside the DataProcessing folder:
    * ExtractApiBasNucInfo.m: used to extract the apical node, basal node, and nucleus center information.
    * CalcCellHeight.m: used to calculate the cell height using the apical and basal nodes of each individual cell.
    * CalcNuclearPositioning.m: used to calculate the nuclear positioning using the apical node, basal node and 
     nucleus center of each individual cell.

Instructions on how to use these files:
1. The ExtractApiBasNucInfo.m must be ran first. It requires inputs from (1) the "detailedStat_" files inside the dataOutput
   folder and (2) the printout file created. More specifically:
     (i) The input for cell information has to be extracted manually from files with a "detailedStat_" prefix. 
         Every triplet (row-wise, i.e. row1 ~ row3) represents the apical node, basal node, and nucleus center 
         position, respectively. For the detail on what information is presented in the files with "detailedStat_" prefix, 
         see detailedStat_example.txt in the dataOutput folder.
     (ii) For simulations involving growth, the order of cells must be entered to acquire the correct calculations. 
          The cell order can be found in the printout created (see EpiScale_run.sh). 
2. Then the files CalcCellHeight.m and CalcNuclearPositioning.m can be ran to extract the cell height and nuclear positioning 
   for any desired simulation.
