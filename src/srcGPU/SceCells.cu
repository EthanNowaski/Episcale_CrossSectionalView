// to do list:
	//1-center of the tissue is manually intorduced in the obtainTwoNewCenter function to recognize if the mothe cell is in front or behind
	//2- If two division occur at the exact same time, there is a chance of error in detecting mother and daughter cells
	//3- In function processMemVec the lateralBefore  and LateralAfter are not carefully assigned. If the code wanted to be used again for the cases where division is happening, this should be revisited.
	//4- If the code wanted to be used again for the case where node deletion is active then the function for calculating cell pressure (void SceCells::calCellPressure()) need to be revisited.
	//5- two bool variables subcellularPolar and cellularPolar are given values inside the code. Although for now it is always true, it is better to be input parameters.
	//6-the value of L0 in the function calAndAddMM_ContractAdh is directly inside the function. It should be an input of the code 
//7- In the function calAndAddMM_ContractRepl, the values of Morse potential are equal to the values of sceIIDiv_M[i] in the input file. it should be an input of the code.
//8- In the function calBendMulti_Mitotic the equlibrium angle for bending stifness is pi it should be an input for the code
//Notes:
	// 1- Currently the nucleus position is desired location not an enforced position. So, all the functions which used "nucleusLocX" & "nucleusLocY" are not active. Instead two variables "nucleusDesireLocX" & "nucleusDesireLocY" are active and internal avg position represent where the nuclei are located.
	// 2- NucleusLocPercent info is currently updated in "copyFirstCellArr_M" using the mother cellrank and daughter cellrank. For now, the daughter cell
	//    will inherit the same percentage as the mother cell. - Kevin
	// 3- We currently force all cell area to be positive after its execution. However, this is only an adhoc method since modifying the CalTriArea anyway will lead 
	//	  to seg fault. Will need to figure out why eventually. - Kevin
	// 4- To accommodate non-convex cell shape during division, the "isMotherCellBehind" is set to be true all the time. Since we are not tracking lineage in anyway at the moment
	//	  , it is acceptable. - Kevin
	// 5- "calApicalBasalRingIds" function is not updated for the newest division algorithm (for non-convex cell) using basal and apical loc of each cell. Since we are not computing contractile ring forces,
	//	  , it is unnecessary. However, if we do indeed need to compute the force, it has to be updated accordingly. - Kevin
	// 6- The 'contractileSpringGrowthProgressSpeed' and similar constants must be recalculated if the min and max distance from the cell center is adjusted for mitotic rounding. - Kevin

#include "SceCells.h"
#include <cmath>
#include <numeric>
//# define debugModeECM 
double epsilon = 1.0e-12;

__constant__ double membrEquLen;
__constant__ double membrStiff;
__constant__ double membrStiff_Mitotic; //Ali June 30
__constant__ double kContractMemb ; 
__constant__ double pI;
__constant__ double minLength;
__constant__ double minDivisor;
__constant__ uint maxAllNodePerCell;
__constant__ uint maxMembrPerCell;
__constant__ uint maxIntnlPerCell;
__constant__ double bendCoeff;
__constant__ double bendCoeff_Mitotic;//AAMIRI

__constant__ double sceIB_M[5];
__constant__ double sceIBDiv_M[5];
__constant__ double sceII_M[5];
__constant__ double sceN_M[5];  //Ali 
__constant__ double sceIIDiv_M[5];
__constant__ double sceNDiv_M[5]; //Ali 
__constant__ double grthPrgrCriEnd_M;
__constant__ double F_Ext_Incline_M2 ;  //Ali



namespace patch{
	template <typename  T> std::string to_string (const T& n) 
	{
	std:: ostringstream stm ; 
	stm << n ; 
	return stm.str() ; 
	}
}


//Ali &  Abu June 30th
__device__
double calMembrForce_Mitotic(double& length, double& progress, double mitoticCri, double adhereIndex) {
/*	if (adhereIndex==-1) {

		if (progress <= mitoticCri) {
			return (length - membrEquLen) * membrStiff;
		} 
		else {
 			return (length - membrEquLen) *(membrStiff+ (membrStiff_Mitotic-membrStiff)* (progress-mitoticCri)/(1.0-mitoticCri));
		}
	}
*/
//	else { 

		if (progress <= mitoticCri) {
			return (length - membrEquLen) * membrStiff;
		} 
		else {
 			return (length - membrEquLen) *(membrStiff+ (membrStiff_Mitotic-membrStiff)* (progress-mitoticCri)/(1.0-mitoticCri));
		}
 
  //       }

}
//
//Ali
__device__
double calMembrForce_Actin(double& length, double kAvg) {
			return ((length - membrEquLen) * kAvg);
}
__device__
double calMembrForce_ActinMitotic(double& length, double kAvg) {
			return ((length - 0.5*membrEquLen) * kAvg);
}

__device__
double CalMembrLinSpringEnergy(double& length, double kAvg) {
			return  (0.5*kAvg *(length - membrEquLen)*(length - membrEquLen)) ;
}



 __device__
double DefaultMembraneStiff() {

	int kStiff=membrStiff ; 
	return kStiff;
		 
}


__device__
double CalExtForce(double  curTime) {
		return min(curTime * F_Ext_Incline_M2,10.0);
}
//Ali
__device__
double obtainRandAngle(uint& cellRank, uint& seed) {
	thrust::default_random_engine rng(seed);
	// discard n numbers to avoid correlation
	rng.discard(cellRank);
	thrust::uniform_real_distribution<double> u0Pi(0, 2.0 * pI);
	double randomAngle = u0Pi(rng);
	return randomAngle;
}

__device__ uint obtainNewIntnlNodeIndex(uint& cellRank, uint& curActiveCount) {
	return (cellRank * maxAllNodePerCell + maxMembrPerCell + curActiveCount);
}

//AAMIRI
__device__ uint obtainLastIntnlNodeIndex(uint& cellRank, uint& curActiveCount) {
	return (cellRank * maxAllNodePerCell + maxMembrPerCell + curActiveCount );
}

//AAMIRI
__device__ uint obtainMembEndNode(uint& cellRank, uint& activeMembrNodeThis) {
	return (cellRank * maxAllNodePerCell + activeMembrNodeThis - 1 );
}


__device__
bool isAllIntnlFilled(uint& currentIntnlCount) {
	if (currentIntnlCount < maxIntnlPerCell) {
		return false;
	} else {
		return true;
	}
}

//AAMIRI

__device__
int obtainRemovingMembrNodeID(uint &cellRank, uint& activeMembrNodes, uint& seed) {
	thrust::default_random_engine rng(seed);
	// discard n numbers to avoid correlation
	rng.discard(activeMembrNodes);
	thrust::uniform_int_distribution<double> dist(0, activeMembrNodes-1);
	int randomNode = dist(rng);
	return (cellRank * maxAllNodePerCell + randomNode);
}


//AAMIRI
__device__
bool isAllIntnlEmptied(uint& currentIntnlCount) {
	if (currentIntnlCount > 0) {
		return false;
	} else {
		return true;
	}
}


//AAMIRI
__device__
bool isAllMembrEmptied(uint& currentMembrCount) {
	if (currentMembrCount > 0) {
		return false;
	} else {
		return true;
	}
}

__device__
bool longEnough(double& length) {
	if (length > minLength) {
		return true;
	} else {
		return false;
	}
}

__device__
double compDist2D(double &xPos, double &yPos, double &xPos2, double &yPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2));
}


void SceCells::distributeBdryIsActiveInfo() {
	thrust::fill(nodes->getInfoVecs().nodeIsActive.begin(),
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile, true);
}


//void SceCells::UpdateTimeStepByAdaptiveMethod( double adaptiveLevelCoef,double minDt,double maxDt, double & dt) {

	//double energyPrime=( energyCell.totalNodeEnergyCell +eCM.energyECM.totalEnergyECM - energyCell.totalNodeEnergyCellOld - eCM.energyECM.totalEnergyECMOld)/dt ; 

	//eCM.energyECM.totalEnergyPrimeECM=( eCM.energyECM.totalEnergyECM  - eCM.energyECM.totalEnergyECMOld)/dt ; 
	//dt=dt ; // max (minDt, maxDt/sqrt( 1 +adaptiveLevelCoef*pow(eCM.energyECM.totalEnergyPrimeECM,2))) ; 
	//dt=max (minDt, maxDt/sqrt( 1 +pow(adaptiveLevelCoef*eCM.energyECM.totalEnergyPrimeECM,2))) ; 
//}


void SceCells::distributeProfileIsActiveInfo() {
	thrust::fill(
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile,
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile
					+ nodes->getAllocPara().currentActiveProfileNodeCount,
			true);
}

void SceCells::distributeECMIsActiveInfo() {
	uint totalNodeCountForActiveECM = allocPara.currentActiveECM
			* allocPara.maxNodePerECM;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveECM);
	thrust::fill(
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosECM,
			nodes->getInfoVecs().nodeIsActive.begin()
					+ totalNodeCountForActiveECM + allocPara.startPosECM, true);
}

void SceCells::distributeCellIsActiveInfo() {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::transform(
			thrust::make_transform_iterator(countingBegin,
					ModuloFunctor(allocPara.maxNodeOfOneCell)),
			thrust::make_transform_iterator(countingEnd,
					ModuloFunctor(allocPara.maxNodeOfOneCell)),
			thrust::make_permutation_iterator(
					cellInfoVecs.activeNodeCountOfThisCell.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosCells,
			thrust::less<uint>());
}

void SceCells::distributeCellGrowthProgress() {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::copy(
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingEnd,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			nodes->getInfoVecs().nodeGrowPro.begin() + allocPara.startPosCells);
}

void MembrPara::initFromConfig() {
	membrEquLenCPU = globalConfigVars.getConfigValue("MembrEquLen").toDouble();
	membrStiffCPU = globalConfigVars.getConfigValue("MembrStiff").toDouble();
	membrStiff_Mitotic = globalConfigVars.getConfigValue("MembrStiff_Mitotic").toDouble();  //Ali June30
	kContractMemb = globalConfigVars.getConfigValue("KContractMemb").toDouble();  //Ali 
	membrGrowCoeff_Ori =
			globalConfigVars.getConfigValue("MembrGrowCoeff").toDouble();
	membrGrowLimit_Ori =
			globalConfigVars.getConfigValue("MembrGrowLimit").toDouble();
	membrGrowCoeff = membrGrowCoeff_Ori;
	membrGrowLimit = membrGrowLimit_Ori;
        //Ali
        F_Ext_Incline = 
			globalConfigVars.getConfigValue("FExtIncline").toDouble();
        //Ali
	membrBendCoeff =
			globalConfigVars.getConfigValue("MembrBenCoeff").toDouble();

//AAMIRI
	membrBendCoeff_Mitotic =
			globalConfigVars.getConfigValue("MembrBenCoeff_Mitotic").toDouble();

	adjustLimit =
			globalConfigVars.getConfigValue("MembrAdjustLimit").toDouble();
	adjustCoeff =
			globalConfigVars.getConfigValue("MembrAdjustCoeff").toDouble();

	growthConst_N =
			globalConfigVars.getConfigValue("MembrGrowthConst").toDouble();
	initMembrCt_N =
			globalConfigVars.getConfigValue("InitMembrNodeCount").toInt();
	initIntnlCt_N =
			globalConfigVars.getConfigValue("InitCellNodeCount").toInt();
}

SceCells::SceCells() {
	//curTime = 0 + 55800.0;//AAMIRI // Ali I comment that out safely on 04/04/2017
    }

void SceCells::growAtRandom(double d_t) {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	// randomly select growth direction and speed.
	randomizeGrowth();

	//std::cout << "after copy grow info" << std::endl;
	updateGrowthProgress();
	//std::cout << "after update growth progress" << std::endl;
	decideIsScheduleToGrow();
	//std::cout << "after decode os schedule to grow" << std::endl;
	computeCellTargetLength();
	//std::cout << "after compute cell target length" << std::endl;
	computeDistToCellCenter();
	//std::cout << "after compute dist to center" << std::endl;
	findMinAndMaxDistToCenter();
	//std::cout << "after find min and max dist" << std::endl;
	computeLenDiffExpCur();
	//std::cout << "after compute diff " << std::endl;
	stretchCellGivenLenDiff();
	//std::cout << "after apply stretch force" << std::endl;
	cellChemotaxis();
	//std::cout << "after apply cell chemotaxis" << std::endl;
	addPointIfScheduledToGrow();
	//std::cout << "after adding node" << std::endl;
}

/**
 * Use the growth magnitude and dt to update growthProgress.
 */
void SceCells::updateGrowthProgress() {
	thrust::transform(cellInfoVecs.growthSpeed.begin(),
			cellInfoVecs.growthSpeed.begin() + allocPara.currentActiveCellCount,
			cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));
}

/**
 * Decide if the cells are going to add a node or not.
 * Use lastCheckPoint and growthProgress to decide whether add point or not
 */
void SceCells::decideIsScheduleToGrow() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			PtCondiOp(miscPara.growThreshold));
}

/**
 * Calculate target length of cell given the cell growth progress.
 * length is along the growth direction.
 */
void SceCells::computeCellTargetLength() {
	thrust::transform(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.expectedLength.begin(),
			CompuTarLen(bioPara.cellInitLength, bioPara.cellFinalLength));
}

/**
 * Compute distance of each node to its corresponding cell center.
 * The distantce could be either positive or negative, depending on the pre-defined
 * growth direction.
 */
void SceCells::computeDistToCellCenter() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(), CompuDist());
}

/**
 * For nodes of each cell, find the maximum and minimum distance to the center.
 * We will then calculate the current length of a cell along its growth direction
 * using max and min distance to the center.
 */
void SceCells::findMinAndMaxDistToCenter() {
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.smallestDistance.begin(), thrust::equal_to<uint>(),
			thrust::minimum<double>());

	// for nodes of each cell, find the maximum distance from the node to the corresponding

	// cell center along the pre-defined growth direction.
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.biggestDistance.begin(), thrust::equal_to<uint>(),
			thrust::maximum<double>());

}

/**
 * Compute the difference for cells between their expected length and current length.
 */
void SceCells::computeLenDiffExpCur() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.lengthDifference.begin(), CompuDiff());
}

/**
 * Use the difference that just computed and growthXDir&growthYDir
 * to apply stretching force (velocity) on nodes of all cells
 */
void SceCells::stretchCellGivenLenDiff() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			ApplyStretchForce(bioPara.elongationCoefficient));
}
/**
 * This is just an attempt. Cells move according to chemicals.
 */
void SceCells::cellChemotaxis() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			ApplyChemoVel(bioPara.chemoCoefficient));
}
/**
 * Adjust the velocities of nodes.
 * For example, velocity of boundary nodes must be zero.
 */
void SceCells::adjustNodeVel() {
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(
			totalNodeCountForActiveCells + allocPara.startPosCells);

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin(),
							countingIterBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin(),
							countingIterBegin)) + totalNodeCountForActiveCells
					+ allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			VelocityModifier(allocPara.startPosProfile,
					allocPara.currentActiveProfileNodeCount));
}
/**
 * Move nodes according to the velocity we just adjusted.
 */
void SceCells::moveNodes() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells + allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			SaxpyFunctorDim2(dt));
}
/**
 * Add a point to a cell if it is scheduled to grow.
 * This step does not guarantee success ; If adding new point failed, it will not change
 * isScheduleToGrow and activeNodeCount;
 */
void SceCells::addPointIfScheduledToGrow() {

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), countingBegin,
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), countingBegin,
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			AddPtOp(allocPara.maxNodeOfOneCell, miscPara.addNodeDistance,
					miscPara.minDistanceToOtherNode,
					growthAuxData.nodeIsActiveAddress,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress, time(NULL),
					miscPara.growThreshold));
}
//Ali commented this constructor in 04/04/2017 // this constructor is not active
SceCells::SceCells(SceNodes* nodesInput,
		std::vector<uint>& numOfInitActiveNodesOfCells,
		std::vector<SceNodeType>& cellTypes) :
		countingBegin(0), initIntnlNodeCount(
				nodesInput->getAllocPara().maxNodeOfOneCell / 2), initGrowthProgress(
				0.0) {
	curTime = 0.0 + 55800.0;//AAMIRI
        std ::cout << "I am in SceCells constructor with polymorphism shape "<<InitTimeStage<<std::endl ; 
	initialize(nodesInput);

	copyInitActiveNodeCount(numOfInitActiveNodesOfCells);

	thrust::device_vector<SceNodeType> cellTypesToPass = cellTypes;
	setCellTypes(cellTypesToPass);

	distributeIsActiveInfo();
}


SceCells::SceCells(SceNodes* nodesInput,SceECM* eCMInput, Solver * solver,
		std::vector<uint>& initActiveMembrNodeCounts,
		std::vector<uint>& initActiveIntnlNodeCounts,
		std::vector<double> &initGrowProgVec, 
		std::vector<ECellType> &eCellTypeV1, 
		double InitTimeStage) {
        curTime=InitTimeStage ; 
        std ::cout << "I am in SceCells constructor with number of inputs "<<InitTimeStage<<std::endl ; 

	tmpDebug = false;
	aniDebug = false;
	membrPara.initFromConfig();
	shrinkRatio = globalConfigVars.getConfigValue("ShrinkRatio").toDouble();
	centerShiftRatio =
			globalConfigVars.getConfigValue("CenterShiftRatio").toDouble();
	double simulationTotalTime =
			globalConfigVars.getConfigValue("SimulationTotalTime").toDouble();

	double simulationTimeStep =
			globalConfigVars.getConfigValue("SimulationTimeStep").toDouble();
	int TotalNumOfOutputFrames =
			globalConfigVars.getConfigValue("TotalNumOfOutputFrames").toInt();

	std ::cout << "I am in SceCells constructor with zero element "<<InitTimeStage<<std::endl ;
	isInitNucPercentCalculated=false ; 
	isBasalActinPresent=true ;
	isCellGrowSet=false ;
	cout <<" Basal actinomyosin is active on pouch cells" << endl ; 
    addNode=true ;
	cout << " addNode boolean is initialized " <<addNode <<endl ; 


	relaxCount=0 ;
	freqPlotData=int ( (simulationTotalTime-InitTimeStage)/(simulationTimeStep*TotalNumOfOutputFrames) ) ; 

	memNewSpacing = globalConfigVars.getConfigValue("MembrLenDiv").toDouble();
	cout << "relax count is initialized as" << relaxCount << endl ; 
	initialize_M(nodesInput, eCMInput, solver);
	copyToGPUConstMem();
	copyInitActiveNodeCount_M(initActiveMembrNodeCounts,
			initActiveIntnlNodeCounts, initGrowProgVec, eCellTypeV1);
}

void SceCells::initCellInfoVecs() {
	cellInfoVecs.growthProgress.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.expectedLength.resize(allocPara.maxCellCount,
			bioPara.cellInitLength);
	cellInfoVecs.lengthDifference.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.smallestDistance.resize(allocPara.maxCellCount);
	cellInfoVecs.biggestDistance.resize(allocPara.maxCellCount);
	cellInfoVecs.activeNodeCountOfThisCell.resize(allocPara.maxCellCount);
	cellInfoVecs.lastCheckPoint.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.isDividing.resize(allocPara.maxCellCount);
	cellInfoVecs.cellTypes.resize(allocPara.maxCellCount, MX);
	cellInfoVecs.isScheduledToGrow.resize(allocPara.maxCellCount, false);
	cellInfoVecs.centerCoordX.resize(allocPara.maxCellCount);
	cellInfoVecs.centerCoordY.resize(allocPara.maxCellCount);
	cellInfoVecs.centerCoordZ.resize(allocPara.maxCellCount);
	cellInfoVecs.cellRanksTmpStorage.resize(allocPara.maxCellCount);
	cellInfoVecs.growthSpeed.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.growthXDir.resize(allocPara.maxCellCount);
	cellInfoVecs.growthYDir.resize(allocPara.maxCellCount);
	cellInfoVecs.isRandGrowInited.resize(allocPara.maxCellCount, false);
}

void SceCells::initCellInfoVecs_M() {
	cellInfoVecs.daughterCellProduced.resize(allocPara_m.maxCellCount, 0);
	cellInfoVecs.distFromNucleus_normal.resize(allocPara_m.maxCellCount,0.0);
	cellInfoVecs.distFromNucleus_normal_apical.resize(allocPara_m.maxCellCount,0.0);
	cellInfoVecs.individualCellHeight.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.individualCellHeight_springHeight.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.distFromBasalLoc.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.distFromApicalLoc.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.contractActomyo_multip_perCell.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.contractActomyo_multip_apical_perCell.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.numApicalVec.resize(allocPara_m.maxCellCount, 0);
	cellInfoVecs.numBasalVec.resize(allocPara_m.maxCellCount, 0);
	cellInfoVecs.cellRankVec.resize(allocPara_m.maxCellCount, 0);
	//std::cout << "max cell count = " << allocPara_m.maxCellCount << std::endl;
	cellInfoVecs.Cell_Damp.resize(allocPara_m.maxCellCount, 36.0);   //Ali
	cellInfoVecs.growthProgress.resize(allocPara_m.maxCellCount, 0.0); //A&A
	cellInfoVecs.contractileSpringGrowthProgress.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.cellAreaGrowthProgress.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.cellAreaGrowthProgressNonMitotic.resize(allocPara_m.maxCellCount, 0.0);
        cellInfoVecs.growthProgressOld.resize(allocPara_m.maxCellCount, 0.0);//Ali
	cellInfoVecs.Cell_Time.resize(allocPara_m.maxCellCount, 0.0); //Ali
	cellInfoVecs.expectedLength.resize(allocPara_m.maxCellCount,
			bioPara.cellInitLength);
	cellInfoVecs.lengthDifference.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.smallestDistance.resize(allocPara_m.maxCellCount);
	cellInfoVecs.biggestDistance.resize(allocPara_m.maxCellCount);
	cellInfoVecs.activeMembrNodeCounts.resize(allocPara_m.maxCellCount);
	cellInfoVecs.activeIntnlNodeCounts.resize(allocPara_m.maxCellCount);
	cellInfoVecs.lastCheckPoint.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.isDividing.resize(allocPara_m.maxCellCount);
	cellInfoVecs.isEnteringMitotic.resize(allocPara_m.maxCellCount, false);  //A&A
	//cellInfoVecs.isRemoving.resize(allocPara.maxCellCount);//AAMIRI
	cellInfoVecs.isScheduledToGrow.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.isScheduledToShrink.resize(allocPara_m.maxCellCount, false);//AAMIRI
	cellInfoVecs.isCellActive.resize(allocPara_m.maxCellCount, false);//AAMIRI
	cellInfoVecs.centerCoordX.resize(allocPara_m.maxCellCount);
	cellInfoVecs.InternalAvgX.resize(allocPara_m.maxCellCount);
	//cellInfoVecs.InternalAvgIniX.resize(allocPara_m.maxCellCount);
	cellInfoVecs.tmpShiftVecX.resize(allocPara_m.maxCellCount);
	cellInfoVecs.centerCoordY.resize(allocPara_m.maxCellCount);
	cellInfoVecs.InternalAvgY.resize(allocPara_m.maxCellCount);
	//cellInfoVecs.InternalAvgIniY.resize(allocPara_m.maxCellCount);
	cellInfoVecs.tmpShiftVecY.resize(allocPara_m.maxCellCount);
	cellInfoVecs.centerCoordZ.resize(allocPara_m.maxCellCount);
	cellInfoVecs.apicalLocX.resize(allocPara_m.maxCellCount);  //Ali 
	cellInfoVecs.apicalLocY.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.basalLocX.resize(allocPara_m.maxCellCount);  //Ali 
	cellInfoVecs.basalLocY.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.eCMNeighborId.resize(allocPara_m.maxCellCount,-1); //Ali 
	cellInfoVecs.nucleusLocX.resize(allocPara_m.maxCellCount);  //Ali 
	cellInfoVecs.nucleusDesireLocX.resize(allocPara_m.maxCellCount);  //Ali 
	cellInfoVecs.nucleusDesireLocY.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.nucDesireDistApical.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.cellCenterDesireDistApical.resize(allocPara_m.maxCellCount); //Kevin
	cellInfoVecs.cellCenterPerturbedLocLocX.resize(allocPara_m.maxCellCount); //Kevin
	cellInfoVecs.cellCenterPerturbedLocLocY.resize(allocPara_m.maxCellCount); //Kevin
	cellInfoVecs.nucleusLocY.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.nucleusLocPercent.resize(allocPara_m.maxCellCount); //Ali 
	cellInfoVecs.apicalNodeCount.resize(allocPara_m.maxCellCount,0); //Ali 
	cellInfoVecs.basalNodeCount.resize(allocPara_m.maxCellCount,0); //Ali 
	cellInfoVecs.ringApicalId.resize(allocPara_m.maxCellCount,-1); //Ali 
	cellInfoVecs.ringBasalId.resize(allocPara_m.maxCellCount,-1); //Ali 
	cellInfoVecs.sumLagrangeFPerCellX.resize(allocPara_m.maxCellCount,0.0); //Ali 
	cellInfoVecs.sumLagrangeFPerCellY.resize(allocPara_m.maxCellCount,0.0); //Ali 
        cellInfoVecs.HertwigXdir.resize(allocPara_m.maxCellCount,0.0); //A&A 
	cellInfoVecs.HertwigYdir.resize(allocPara_m.maxCellCount,0.0); //A&A 

	cellInfoVecs.cellRanksTmpStorage.resize(allocPara_m.maxCellCount);
	cellInfoVecs.cellRanksTmpStorage1.resize(allocPara_m.maxCellCount);
	cellInfoVecs.growthSpeed.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.growthXDir.resize(allocPara_m.maxCellCount);
	cellInfoVecs.growthYDir.resize(allocPara_m.maxCellCount);
	cellInfoVecs.isRandGrowInited.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.isMembrAddingNode.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.isMembrRemovingNode.resize(allocPara_m.maxCellCount, false); // Ali
	cellInfoVecs.maxTenIndxVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.minTenIndxVec.resize(allocPara_m.maxCellCount); //Ali
	cellInfoVecs.maxTenRiVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.maxDistToRiVec.resize(allocPara_m.maxCellCount); //Ali
	cellInfoVecs.maxTenIndxTypeVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.minDistToRiVec.resize(allocPara_m.maxCellCount); //Ali
	cellInfoVecs.maxTenRiMidXVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.maxTenRiMidYVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.aveTension.resize(allocPara_m.maxCellCount);
	cellInfoVecs.membrGrowProgress.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.membrGrowSpeed.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.cellAreaVec.resize(allocPara_m.maxCellCount, 0.0);
    cellInfoVecs.cellPerimVec.resize(allocPara_m.maxCellCount, 0.0);//AAMIRI
    cellInfoVecs.cellPressure.resize(allocPara_m.maxCellCount, 0.0);//Ali
    cellInfoVecs.sumF_MI_M_N.resize(allocPara_m.maxCellCount, 0.0);//Ali
    cellInfoVecs.sumLagrangeFN.resize(allocPara_m.maxCellCount, 0.0);//Ali
    cellInfoVecs.eCellTypeV2.resize(allocPara_m.maxCellCount, notActive);//Ali 
    //cellInfoVecs.eCellTypeV2Host.resize(allocPara_m.maxCellCount, notActive);//Ali 
    cellInfoVecs.cellRoot.resize(allocPara_m.maxCellCount, -1);//Ali

	thrust:: sequence (cellInfoVecs.cellRoot.begin(),cellInfoVecs.cellRoot.begin()+allocPara_m.currentActiveCellCount) ; //Ali
        std::cout << "initial number of active cells is " <<allocPara_m.currentActiveCellCount <<std::endl;
	    std::cout <<"last cell rank used in the cell root is " <<cellInfoVecs.cellRoot[allocPara_m.currentActiveCellCount-1] << endl ;   
}

void SceCells::initCellNodeInfoVecs() {
	cellNodeInfoVecs.cellRanks.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeXPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeYPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeZPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.distToCenterAlongGrowDir.resize(
			allocPara.maxTotalCellNodeCount);
}

void SceCells::initCellNodeInfoVecs_M() {
	std::cout << "max total node count = " << allocPara_m.maxTotalNodeCount
			<< std::endl;
	cellNodeInfoVecs.cellRanks.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeXPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeYPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeZPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.distToCenterAlongGrowDir.resize(
			allocPara_m.maxTotalNodeCount);

	cellNodeInfoVecs.activeLocXApical.resize(allocPara_m.maxTotalNodeCount); //Ali 
	cellNodeInfoVecs.activeLocYApical.resize(allocPara_m.maxTotalNodeCount); //Ali

	cellNodeInfoVecs.activeLocXBasal.resize(allocPara_m.maxTotalNodeCount); //Ali 
	cellNodeInfoVecs.activeLocYBasal.resize(allocPara_m.maxTotalNodeCount); //Ali 
}

void SceCells::initGrowthAuxData() {
	growthAuxData.nodeIsActiveAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[allocPara.startPosCells]));
	growthAuxData.nodeXPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[allocPara.startPosCells]));
	growthAuxData.nodeYPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[allocPara.startPosCells]));
	growthAuxData.randomGrowthSpeedMin = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMin").toDouble();
	growthAuxData.randomGrowthSpeedMax = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMax").toDouble();
	growthAuxData.randGenAuxPara = globalConfigVars.getConfigValue(
			"RandomGenerationAuxPara").toDouble();

	if (controlPara.simuType == SingleCellTest) {
		growthAuxData.fixedGrowthSpeed = globalConfigVars.getConfigValue(
				"FixedGrowthSpeed").toDouble();
	}
}

void SceCells::initGrowthAuxData_M() {
	growthAuxData.nodeIsActiveAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[allocPara_m.bdryNodeCount]));
	growthAuxData.nodeXPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[allocPara_m.bdryNodeCount]));
	growthAuxData.nodeYPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[allocPara_m.bdryNodeCount]));
	growthAuxData.adhIndxAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeAdhereIndex[allocPara_m.bdryNodeCount]));

	growthAuxData.memNodeType1Address = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().memNodeType1[allocPara_m.bdryNodeCount])); //Ali 
	growthAuxData.randomGrowthSpeedMin_Ori = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMin").toDouble();
	growthAuxData.randomGrowthSpeedMax_Ori = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMax").toDouble();
	growthAuxData.randomGrowthSpeedMin = growthAuxData.randomGrowthSpeedMin_Ori;
	growthAuxData.randomGrowthSpeedMax = growthAuxData.randomGrowthSpeedMax_Ori;
	growthAuxData.grthPrgrCriVal_M_Ori = globalConfigVars.getConfigValue(
			"GrowthPrgrCriVal").toDouble();
	growthAuxData.grthProgrEndCPU = globalConfigVars.getConfigValue(
			"GrowthPrgrValEnd").toDouble();
}

void SceCells::initialize(SceNodes* nodesInput) {
	nodes = nodesInput;
	controlPara = nodes->getControlPara();
	readMiscPara();
	readBioPara();

	allocPara = nodesInput->getAllocPara();

	// max internal node count must be even number.
	assert(allocPara_m.maxIntnlNodePerCell % 2 == 0);

	initCellInfoVecs();
	initCellNodeInfoVecs();
	initGrowthAuxData();

	distributeIsCellRank();
}

void SceCells::initialize_M(SceNodes* nodesInput, SceECM *eCMInput, Solver *solver) {
	std::cout << "Initializing cells ...... " << std::endl;
	//std::cout.flush();
	nodes = nodesInput; //pointer assigned
	eCMPointerCells=eCMInput ; //pointer assigned 
	solverPointer=solver ; 
	allocPara_m = nodesInput->getAllocParaM();
	// max internal node count must be even number.
	assert(allocPara_m.maxIntnlNodePerCell % 2 == 0);

	//std::cout << "break point 1 " << std::endl;
	//std::cout.flush();
	controlPara = nodes->getControlPara();  // It copies the controlPara from nstance of class SceNodes to the instance of class of SceCells
	//std::cout << "break point 2 " << std::endl;
	//std::cout.flush();
	readMiscPara_M();
	//std::cout << "break point 3 " << std::endl;
	//std::cout.flush();
	initCellInfoVecs_M();

	//std::cout << "break point 4 " << std::endl;
	//std::cout.flush();
	readBioPara();
	//std::cout << "break point 5 " << std::endl;
	//std::cout.flush();

	//std::cout << "break point 6 " << std::endl;
	//std::cout.flush();
	initCellNodeInfoVecs_M();
	//std::cout << "break point 7 " << std::endl;
	//std::cout.flush();
	initGrowthAuxData_M();
	//std::cout << "break point 8 " << std::endl;
	//std::cout.flush();

}

void SceCells::copyInitActiveNodeCount(
		std::vector<uint>& numOfInitActiveNodesOfCells) {
	thrust::copy(numOfInitActiveNodesOfCells.begin(),
			numOfInitActiveNodesOfCells.end(),
			cellInfoVecs.activeNodeCountOfThisCell.begin());
}

void SceCells::allComponentsMove() {
	adjustNodeVel();
	moveNodes();
}

/**
 * Mark cell node as either activdistributeIsActiveInfo()e or inactive.
 * left part of the node array will be active and right part will be inactive.
 * the threshold is defined by array activeNodeCountOfThisCell.
 * e.g. activeNodeCountOfThisCell = {2,3} and  maxNodeOfOneCell = 5
 */
void SceCells::distributeIsActiveInfo() {
//std::cout << "before distribute bdry isActive" << std::endl;
	distributeBdryIsActiveInfo();
//std::cout << "before distribute profile isActive" << std::endl;
	distributeProfileIsActiveInfo();
//std::cout << "before distribute ecm isActive" << std::endl;
	distributeECMIsActiveInfo();
//std::cout << "before distribute cells isActive" << std::endl;
	distributeCellIsActiveInfo();
}

void SceCells::distributeIsCellRank() {
	uint totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingCellEnd(
			totalNodeCountForActiveCells);
	std::cerr << "totalNodeCount for active cells "
			<< totalNodeCountForActiveCells << std::endl;

//thrust::counting_iterator<uint> countingECMEnd(countingECMEnd);

// only computes the cell ranks of cells. the rest remain unchanged.
	thrust::transform(countingBegin, countingCellEnd,
			nodes->getInfoVecs().nodeCellRank.begin() + allocPara.startPosCells,
			DivideFunctor(allocPara.maxNodeOfOneCell));
	std::cerr << "finished cellRank transformation" << std::endl;
}

/**
 * This method computes center of all cells.
 * more efficient then simply iterating the cell because of parallel reducing.
 */
void SceCells::computeCenterPos() {
	uint totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	uint totalNumberOfActiveNodes = thrust::reduce(
			cellInfoVecs.activeNodeCountOfThisCell.begin(),
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin(),
							cellNodeInfoVecs.activeZPoss.begin())), isTrue());

	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalNumberOfActiveNodes,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin(),
							cellNodeInfoVecs.activeZPoss.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())),
			thrust::equal_to<uint>(), CVec3Add());
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.activeNodeCountOfThisCell.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())), CVec3Divide());
}

/**
 * 2D version of cell division.
 * Division process is done by creating two temporary vectors to hold the node information
 * that are going to divide.
 *
 * step 1: based on lengthDifference, expectedLength and growthProgress,
 *     this process determines whether a certain cell is ready to divide and then assign
 *     a boolean value to isDivided.
 *
 * step 2. copy those cells that will divide in to the temp vectors created
 *
 * step 3. For each cell in the temp vectors, we sort its nodes by its distance to the
 * corresponding cell center.
 * This step is not very effcient when the number of cells going to divide is big.
 * but this is unlikely to happen because cells will divide according to external chemical signaling
 * and each will have different divide progress.
 *
 * step 4. copy the right part of each cell of the sorted array (temp1) to left part of each cell of
 * another array
 *
 * step 5. transform isActive vector of both temp1 and temp2, making only left part of each cell active.
 *
 * step 6. insert temp2 to the end of the cell array
 *
 * step 7. copy temp1 to the previous position of the cell array.
 *
 * step 8. add activeCellCount of the system.
 *
 * step 9. mark isDivide of all cells to false.
 */

void SceCells::divide2DSimplified() {
	bool isDivisionPresent = decideIfGoingToDivide();
	if (!isDivisionPresent) {
		return;
	}
	copyCellsPreDivision();
	sortNodesAccordingToDist();
	copyLeftAndRightToSeperateArrays();
	transformIsActiveArrayOfBothArrays();
	addSecondArrayToCellArray();
	copyFirstArrayToPreviousPos();
	updateActiveCellCount();
	markIsDivideFalse();
}

bool SceCells::decideIfGoingToDivide() {
// step 1
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lengthDifference.begin(),
							cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lengthDifference.begin(),
							cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin()))
					+ allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isDividing.begin(),
							cellInfoVecs.growthProgress.begin())),
			CompuIsDivide(miscPara.isDivideCriticalRatio,
					allocPara.maxNodeOfOneCell));
// sum all bool values which indicate whether the cell is going to divide.
// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toBeDivideCount = thrust::reduce(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin() + allocPara.currentActiveCellCount,
			(uint) (0));
	if (divAuxData.toBeDivideCount > 0) {
		return true;
	} else {
		return false;
	}
}

void SceCells::copyCellsPreDivision() {
// step 2 : copy all cell rank and distance to its corresponding center with divide flag = 1
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
			* allocPara.maxNodeOfOneCell;
	divAuxData.tmpIsActiveHold1 = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, true);
	divAuxData.tmpDistToCenter1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpCellRankHold1 = thrust::device_vector<uint>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpXValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpZValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

	divAuxData.tmpCellTypes = thrust::device_vector<SceNodeType>(
			divAuxData.nodeStorageCount);

	divAuxData.tmpIsActiveHold2 = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpDistToCenter2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpXValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpZValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

// step 2 , continued
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeCellType.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeCellType.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(cellInfoVecs.isDividing.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpCellRankHold1.begin(),
							divAuxData.tmpDistToCenter1.begin(),
							divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin(),
							divAuxData.tmpCellTypes.begin())), isTrue());
}

/**
 * performance wise, this implementation is not the best because I can use only one sort_by_key
 * with speciialized comparision operator. However, This implementation is more robust and won't
 * compromise performance too much.
 */
void SceCells::sortNodesAccordingToDist() {
//step 3
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		thrust::sort_by_key(
				divAuxData.tmpDistToCenter1.begin()
						+ i * allocPara.maxNodeOfOneCell,
				divAuxData.tmpDistToCenter1.begin()
						+ (i + 1) * allocPara.maxNodeOfOneCell,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								divAuxData.tmpXValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell,
								divAuxData.tmpYValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell,
								divAuxData.tmpZValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell)));
	}
}

/**
 * scatter_if() is a thrust function.
 * inputIter1 first,
 * inputIter1 last,
 * inputIter2 map,
 * inputIter3 stencil
 * randomAccessIter output
 */
void SceCells::copyLeftAndRightToSeperateArrays() {
//step 4.
	thrust::scatter_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold1.end(),
							divAuxData.tmpYValueHold1.end(),
							divAuxData.tmpZValueHold1.end())),
			make_transform_iterator(countingBegin,
					LeftShiftFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					IsRightSide(allocPara.maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold2.begin(),
							divAuxData.tmpYValueHold2.begin(),
							divAuxData.tmpZValueHold2.begin())));
}

void SceCells::transformIsActiveArrayOfBothArrays() {
	thrust::transform(countingBegin,
			countingBegin + divAuxData.nodeStorageCount,
			divAuxData.tmpIsActiveHold1.begin(),
			IsLeftSide(allocPara.maxNodeOfOneCell));
	thrust::transform(countingBegin,
			countingBegin + divAuxData.nodeStorageCount,
			divAuxData.tmpIsActiveHold2.begin(),
			IsLeftSide(allocPara.maxNodeOfOneCell));
	if (divAuxData.toBeDivideCount != 0) {
		std::cout << "before insert, active cell count in nodes:"
				<< nodes->getAllocPara().currentActiveCellCount << std::endl;
	}
}

void SceCells::addSecondArrayToCellArray() {
/// step 6. call SceNodes function to add newly divided cells
	nodes->addNewlyDividedCells(divAuxData.tmpXValueHold2,
			divAuxData.tmpYValueHold2, divAuxData.tmpZValueHold2,
			divAuxData.tmpIsActiveHold2, divAuxData.tmpCellTypes);
}

void SceCells::copyFirstArrayToPreviousPos() {
	thrust::scatter(
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActiveHold1.begin(),
							divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActiveHold1.end(),
							divAuxData.tmpXValueHold1.end(),
							divAuxData.tmpYValueHold1.end(),
							divAuxData.tmpZValueHold1.end())),
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(countingBegin,
									divAuxData.tmpCellRankHold1.begin())),
					CompuPos(allocPara.maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells)));

	/**
	 * after dividing, the cell should resume the initial
	 * (1) node count, which defaults to be half size of max node count
	 * (2) growth progress, which defaults to 0
	 * (3) last check point, which defaults to 0
	 */
	thrust::scatter_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(initIntnlNodeCount, initGrowthProgress,
							initGrowthProgress)),
			thrust::make_zip_iterator(
					thrust::make_tuple(initIntnlNodeCount, initGrowthProgress,
							initGrowthProgress))
					+ allocPara.currentActiveCellCount, countingBegin,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())), isTrue());

// TODO: combine this one with the previous scatter_if to improve efficiency.
	thrust::fill(
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount
					+ divAuxData.toBeDivideCount,
			allocPara.maxNodeOfOneCell / 2);

}

void SceCells::updateActiveCellCount() {
	allocPara.currentActiveCellCount = allocPara.currentActiveCellCount
			+ divAuxData.toBeDivideCount;
	NodeAllocPara para = nodes->getAllocPara();
	para.currentActiveCellCount = allocPara.currentActiveCellCount;
	nodes->setAllocPara(para);
}

void SceCells::markIsDivideFalse() {
	thrust::fill(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin() + allocPara.currentActiveCellCount,
			false);
}

void SceCells::readMiscPara() {
	miscPara.addNodeDistance = globalConfigVars.getConfigValue(
			"DistanceForAddingNode").toDouble();
	miscPara.minDistanceToOtherNode = globalConfigVars.getConfigValue(
			"MinDistanceToOtherNode").toDouble();
	miscPara.isDivideCriticalRatio = globalConfigVars.getConfigValue(
			"IsDivideCrticalRatio").toDouble();
// reason for adding a small term here is to avoid scenario when checkpoint might add many times
// up to 0.99999999 which is theoretically 1.0 but not in computer memory. If we don't include
// this small term we might risk adding one more node.
	int maxNodeOfOneCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toInt();
	miscPara.growThreshold = 1.0 / (maxNodeOfOneCell - maxNodeOfOneCell / 2)
			+ epsilon;
}

void SceCells::readMiscPara_M() {
	miscPara.addNodeDistance = globalConfigVars.getConfigValue(
			"DistanceForAddingNode").toDouble();
	miscPara.minDistanceToOtherNode = globalConfigVars.getConfigValue(
			"MinDistanceToOtherNode").toDouble();
	miscPara.isDivideCriticalRatio = globalConfigVars.getConfigValue(
			"IsDivideCrticalRatio").toDouble();
// reason for adding a small term here is to avoid scenario when checkpoint might add many times
// up to 0.99999999 which is theoretically 1.0 but not in computer memory. If we don't include
// this small term we might risk adding one more node.
	int maxIntnlNodePerCell = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();
	miscPara.growThreshold = 1.0
			/ (maxIntnlNodePerCell - maxIntnlNodePerCell / 2) + epsilon;
	miscPara.prolifDecayCoeff = globalConfigVars.getConfigValue(
			"ProlifDecayCoeff").toDouble();
}

void SceCells::readBioPara() {
	if (controlPara.simuType != Disc_M) {
		bioPara.cellInitLength = globalConfigVars.getConfigValue(
				"CellInitLength").toDouble();
		std::cout << "break point 1 " << bioPara.cellInitLength << std::endl;
		std::cout.flush();
		bioPara.cellFinalLength = globalConfigVars.getConfigValue(
				"CellFinalLength").toDouble();
		std::cout << "break point 2 " << bioPara.cellFinalLength << std::endl;
		std::cout.flush();
		bioPara.elongationCoefficient = globalConfigVars.getConfigValue(
				"ElongateCoefficient").toDouble();

		std::cout << "break point 3 " << bioPara.elongationCoefficient
				<< std::endl;
		std::cout.flush();

	}

	if (controlPara.simuType == Beak) {
		std::cout << "break point 4 " << std::endl;
		std::cout.flush();
		bioPara.chemoCoefficient = globalConfigVars.getConfigValue(
				"ChemoCoefficient").toDouble();
	}
	//std::cin >> jj;
}

void SceCells::randomizeGrowth() {
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)) + allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin())),
			AssignRandIfNotInit(growthAuxData.randomGrowthSpeedMin,
					growthAuxData.randomGrowthSpeedMax,
					allocPara.currentActiveCellCount,
					growthAuxData.randGenAuxPara));
}

/**
 * To run all the cell level logics.
 * First step we got center positions of cells.
 * Grow.
 */
void SceCells::runAllCellLevelLogicsDisc(double dt) {
	this->dt = dt;
//std::cerr << "enter run all cell level logics" << std::endl;
	computeCenterPos();
//std::cerr << "after compute center position." << std::endl;

	if (nodes->getControlPara().controlSwitchs.stab == OFF) {
		growAtRandom(dt);
		//grow2DTwoRegions(dt, region1, region2);
		//std::cerr << "after grow cells" << std::endl;
		//distributeIsActiveInfo();
		//std::cerr << "after distribute is active info." << std::endl;
		divide2DSimplified();
		//std::cerr << "after divide 2D simplified." << std::endl;
		distributeIsActiveInfo();
		//std::cerr << "after distribute is active info." << std::endl;
		distributeCellGrowthProgress();
	}

	allComponentsMove();
	//std::cerr << "after all components move." << std::endl;
}

//Ali void SceCells::runAllCellLogicsDisc_M(double dt) {
void SceCells::runAllCellLogicsDisc_M(double & dt, double Damp_Coef, double InitTimeStage, 
										double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_ECM, double timeRatio_Crit_Division,
											double volume_Increase_Target_Ratio, double volume_Increase_Scale, double postDivision_restorationRateScale, int cycle,
											double distFromNucleus_max, double distFromNucleus_min, double distFromNucleus_normalMax1,double distFromNucleus_normalMax2,double distFromNucleus_normalMax3,
											double distFromNucleus_normalMax_apical1, double distFromNucleus_normalMax_apical2, double distFromNucleus_normalMax_apical3,
											double percentage_before_timeRatio_Crit_Division_scaling, double growthProgressSpeed,
											int maxApicalBasalNodeNum, double maxLengthToAddMemNodes, double mitoRndActomyoStrengthScaling, double thresholdToIntroduceNewCell, double mitoticThreshold) {   //Ali

// double time_Ratio = current_Time/total_Time;
// if (timeRatio == timeRatio_Crit_Division){
// 	std::cout<<"timeRatio = "<<timeRatio<<std::endl;
// }

	// if (cycle < 0){
	// 	divAuxData.cellRank_division = -1;
	// 	divAuxData.cellRank_division2 = -1;//allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle >= 10){
	// 	divAuxData.cellRank_division = -1;
	// 	divAuxData.cellRank_division2 = -1;//allocPara_m.currentActiveCellCount;
	// }

	// if (cycle == 0 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 31;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 1 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 25;//33;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 2 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 6;//27;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 3 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 10;//29;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 4 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 61;//35;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 5 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 30;//86;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 6 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 87;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 7 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 5;//88;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 8 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 42;//89;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else if (cycle == 9 && timeRatio < timeRatio_Crit_Division){
	// 	divAuxData.cellRank_division = 37;//90;
	// 	divAuxData.cellRank_division2 = allocPara_m.currentActiveCellCount;
	// }
	// else{
	// 	divAuxData.cellRank_division = -1;
	// 	divAuxData.cellRank_division2 = -1;//allocPara_m.currentActiveCellCount;
	// }

#ifdef debugModeECM 
	cudaEvent_t start1, start2, start3, start4, start5, start6, start7, start8, start9, start10, start11, start12, start13, stop;
	float elapsedTime1, elapsedTime2, elapsedTime3, elapsedTime4, elapsedTime5, elapsedTime6,  elapsedTime7 , elapsedTime8 ; 
	float elapsedTime9, elapsedTime10, elapsedTime11, elapsedTime12, elapsedTime13  ; 
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&start3);
	cudaEventCreate(&start4);
	cudaEventCreate(&start5);
	cudaEventCreate(&start6);
	cudaEventCreate(&start7);
	cudaEventCreate(&start8);
	cudaEventCreate(&start9);
	cudaEventCreate(&start10);
	cudaEventCreate(&start11);
	cudaEventCreate(&start12);
	cudaEventCreate(&start13);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start1, 0);
#endif

	// std::cout << "     *** 1 ***" << endl;
	this->dt = dt;
    this->Damp_Coef=Damp_Coef ; //Ali 
    this->InitTimeStage=InitTimeStage   ;  //A & A 
	growthAuxData.prolifDecay = exp(-curTime * miscPara.prolifDecayCoeff);
        //cout<< "Current Time in simulation is: "<<curTime <<endl; 
	growthAuxData.randomGrowthSpeedMin = growthAuxData.prolifDecay
			* growthAuxData.randomGrowthSpeedMin_Ori;
	growthAuxData.randomGrowthSpeedMax = growthAuxData.prolifDecay
			* growthAuxData.randomGrowthSpeedMax_Ori;

	bool cellPolar=true ; 
	bool subCellPolar= true  ; 
	// std::cout << "     *** 2 ***" << endl;
	if (curTime>500000) {
	//	eCMPointerCells->SetIfECMIsRemoved(false) ; 
	//	isBasalActinPresent=false ; 
	//	nodes->SetApicalAdhPresence(true) ; 
	}

 	if (curTime==InitTimeStage) {
		lastPrintNucleus=10000000  ; //just a big number 
		outputFrameNucleus=0 ;
	//	computeInternalAvgPos_M();
		nodes->isInitPhase=false ; // This bool variable is not active in the code anymore

		string uniqueSymbolOutput =
		globalConfigVars.getConfigValue("UniqueSymbol").toString();

		std::string cSVFileName = "EnergyExportCell_" + uniqueSymbolOutput + ".CSV";
		ofstream EnergyExportCell ;
		EnergyExportCell.open(cSVFileName.c_str() );
		EnergyExportCell <<"curTime"<<","<<"totalMembrLinSpringEnergyCell" << "," <<"totalMembrBendSpringEnergyCell" <<"," <<
		"totalNodeIIEnergyCell"<<"," <<"totalNodeIMEnergyCell"<<","<<"totalNodeEnergyCell"<< std::endl;

// Adding cell cycle timer information here

		std::vector<double> growthProgressReinitialize;
		growthProgressReinitialize.push_back(-9999.9);//0.0163);//Cell0
		growthProgressReinitialize.push_back(-9999.9);//0.0181);//Cell1
// Ver7
// growthProgressReinitialize.push_back(	 -99999.9);
// growthProgressReinitialize.push_back(    -99999.9);
// growthProgressReinitialize.push_back(   -0.58);
// growthProgressReinitialize.push_back(   -2.69);
// growthProgressReinitialize.push_back(   -1.2);
// growthProgressReinitialize.push_back(    0.1100);
// growthProgressReinitialize.push_back(   -2.978);
// growthProgressReinitialize.push_back(   -1.0070);
// growthProgressReinitialize.push_back(   -2.8052);
// growthProgressReinitialize.push_back(   -2.0354);
// growthProgressReinitialize.push_back(   -1.1824);
// growthProgressReinitialize.push_back(   -0.6577);
// growthProgressReinitialize.push_back(    0.6080);
// growthProgressReinitialize.push_back(   -1.0282);
// growthProgressReinitialize.push_back(    -0.2033);
// growthProgressReinitialize.push_back(   -2.34);
// growthProgressReinitialize.push_back(    0.0127);
// growthProgressReinitialize.push_back(   -0.2035);
// growthProgressReinitialize.push_back(   -1.299);
// growthProgressReinitialize.push_back(   -2.703);
// growthProgressReinitialize.push_back(   -1.562);
// growthProgressReinitialize.push_back(   -2.750);
// growthProgressReinitialize.push_back(   -0.3286);
// growthProgressReinitialize.push_back(   -2.083);
// growthProgressReinitialize.push_back(   -2.79);
// growthProgressReinitialize.push_back(   -1.1567);
// growthProgressReinitialize.push_back(   -0.5034);
// growthProgressReinitialize.push_back(   -1.9003);
// growthProgressReinitialize.push_back(    0.4964);
// growthProgressReinitialize.push_back(   -0.4520);
// growthProgressReinitialize.push_back(   -2.002);
// growthProgressReinitialize.push_back(   -1.000);
// growthProgressReinitialize.push_back(   -1.880);
// growthProgressReinitialize.push_back(    0.3719);
// growthProgressReinitialize.push_back(   -0.7133);
// growthProgressReinitialize.push_back(   -1.172);
// growthProgressReinitialize.push_back(    0.0251);
// growthProgressReinitialize.push_back(   -2.323);
// growthProgressReinitialize.push_back(   -1.960);
// growthProgressReinitialize.push_back(   -0.1294);
// growthProgressReinitialize.push_back(    0.2848);
// growthProgressReinitialize.push_back(   -2.912);
// growthProgressReinitialize.push_back(    0.2526);
// growthProgressReinitialize.push_back(   -2.165);
// growthProgressReinitialize.push_back(   -1.031);
// growthProgressReinitialize.push_back(   -0.7257);
// growthProgressReinitialize.push_back(   -2.087);
// growthProgressReinitialize.push_back(   -1.018);
// growthProgressReinitialize.push_back(    0.0391);
// growthProgressReinitialize.push_back(   -2.1332);
// growthProgressReinitialize.push_back(   -3.2330);
// growthProgressReinitialize.push_back(   -0.3449);
// growthProgressReinitialize.push_back(   -2.0334);
// growthProgressReinitialize.push_back(   -0.0101);
// growthProgressReinitialize.push_back(    0.4452);
// growthProgressReinitialize.push_back(   -2.013);
// growthProgressReinitialize.push_back(    0.0002);
// growthProgressReinitialize.push_back(   -1.048);
// growthProgressReinitialize.push_back(    0.2862);
// growthProgressReinitialize.push_back(   -9999.9);
// growthProgressReinitialize.push_back(   -9999.9);	

growthProgressReinitialize.push_back( 	-9999.9); //Cell 2
growthProgressReinitialize.push_back(   -9999.9); //Cell 3

// //Ver8
// std::cout<<"growthProgress initial profile Ver. 8"<<std::endl;
// growthProgressReinitialize.push_back(   -1.0825);
// growthProgressReinitialize.push_back(   -3.4494);
// growthProgressReinitialize.push_back(   -2.0672);
// growthProgressReinitialize.push_back(   -2.8107);
// growthProgressReinitialize.push_back(   -0.1243);
// growthProgressReinitialize.push_back(   -2.1773);
// growthProgressReinitialize.push_back(   -1.2537);
// growthProgressReinitialize.push_back(   -2.7960);
// growthProgressReinitialize.push_back(   -0.9416);
// growthProgressReinitialize.push_back(   -2.3824);
// growthProgressReinitialize.push_back(   -0.7202);
// growthProgressReinitialize.push_back(   -1.5852);
// growthProgressReinitialize.push_back(   -3.1438);
// growthProgressReinitialize.push_back(   -2.5268);
// growthProgressReinitialize.push_back(    0.3817);
// growthProgressReinitialize.push_back(   -2.8524);
// growthProgressReinitialize.push_back(    0.0097);
// growthProgressReinitialize.push_back(   -1.2120);
// growthProgressReinitialize.push_back(    0.7336);
// growthProgressReinitialize.push_back(   -3.1678);
// growthProgressReinitialize.push_back(   -1.6186);
// growthProgressReinitialize.push_back(   -3.0467);
// growthProgressReinitialize.push_back(    0.5881);
// growthProgressReinitialize.push_back(   -3.4803);
// growthProgressReinitialize.push_back(   -0.2066);
// growthProgressReinitialize.push_back(   -3.1411);
// growthProgressReinitialize.push_back(   -1.8009);
// growthProgressReinitialize.push_back(   -2.3956);
// growthProgressReinitialize.push_back(   -0.0997);
// growthProgressReinitialize.push_back(   -1.6665);
// growthProgressReinitialize.push_back(    0.3703);
// growthProgressReinitialize.push_back(   -2.7272);
// growthProgressReinitialize.push_back(    0.1945);
// growthProgressReinitialize.push_back(   -1.0363);
// growthProgressReinitialize.push_back(   -2.8839);
// growthProgressReinitialize.push_back(    0.1254);
// growthProgressReinitialize.push_back(   -0.8563);
// growthProgressReinitialize.push_back(   -2.0085);
// growthProgressReinitialize.push_back(   -1.3187);
// growthProgressReinitialize.push_back(   -3.1771);
// growthProgressReinitialize.push_back(   -2.4804);
// growthProgressReinitialize.push_back(   -1.7266);
// growthProgressReinitialize.push_back(   -3.2890);
// growthProgressReinitialize.push_back(    0.3365);
// growthProgressReinitialize.push_back(   -1.4138);
// growthProgressReinitialize.push_back(   -2.0647);
// growthProgressReinitialize.push_back(    0.3252);
// growthProgressReinitialize.push_back(   -1.9307);
// growthProgressReinitialize.push_back(   -3.0274);
// growthProgressReinitialize.push_back(   -0.1839);
// growthProgressReinitialize.push_back(   -1.8436);
// growthProgressReinitialize.push_back(   -2.4728);
// growthProgressReinitialize.push_back(   -1.7834);
// growthProgressReinitialize.push_back(   -3.0901);
// growthProgressReinitialize.push_back(    0.5037);
// growthProgressReinitialize.push_back(   -1.0554);
// growthProgressReinitialize.push_back(   -3.2459);

// //Ver9
// std::cout<<"growthProgress initial profile Ver. 9"<<std::endl;
// growthProgressReinitialize.push_back(-1.3519);
// growthProgressReinitialize.push_back(-0.2639);
// growthProgressReinitialize.push_back(-3.1180);
// growthProgressReinitialize.push_back(0.7373);
// growthProgressReinitialize.push_back(-1.8542);
// growthProgressReinitialize.push_back(-0.4155);
// growthProgressReinitialize.push_back(-3.1281);
// growthProgressReinitialize.push_back(0.6909);
// growthProgressReinitialize.push_back(-1.9936);
// growthProgressReinitialize.push_back(-0.6243);
// growthProgressReinitialize.push_back(-3.1490);
// growthProgressReinitialize.push_back(0.7337);
// growthProgressReinitialize.push_back(-1.9932);
// growthProgressReinitialize.push_back(-0.4141);
// growthProgressReinitialize.push_back(-3.4583);
// growthProgressReinitialize.push_back(0.7396);
// growthProgressReinitialize.push_back(-1.4360);
// growthProgressReinitialize.push_back(-0.3950);
// growthProgressReinitialize.push_back(-3.3643);
// growthProgressReinitialize.push_back(0.6351);
// growthProgressReinitialize.push_back(-1.9657);
// growthProgressReinitialize.push_back(-0.8638);
// growthProgressReinitialize.push_back(-3.4182);
// growthProgressReinitialize.push_back(0.6429);
// growthProgressReinitialize.push_back(-2.3591);
// growthProgressReinitialize.push_back(-0.3721);
// growthProgressReinitialize.push_back(-3.4084);
// growthProgressReinitialize.push_back(0.6688);
// growthProgressReinitialize.push_back(-2.2819);
// growthProgressReinitialize.push_back(-0.9605);
// growthProgressReinitialize.push_back(-3.3901);
// growthProgressReinitialize.push_back(0.4828);
// growthProgressReinitialize.push_back(-2.2890);
// growthProgressReinitialize.push_back(-0.6388);
// growthProgressReinitialize.push_back(-3.4209);
// growthProgressReinitialize.push_back(0.5414);
// growthProgressReinitialize.push_back(-1.7997);
// growthProgressReinitialize.push_back(-0.5446);
// growthProgressReinitialize.push_back(-2.8349);
// growthProgressReinitialize.push_back(0.5172);
// growthProgressReinitialize.push_back(-1.5476);
// growthProgressReinitialize.push_back(-0.4978);
// growthProgressReinitialize.push_back(-3.1851);
// growthProgressReinitialize.push_back(0.5574);
// growthProgressReinitialize.push_back(-1.6700);
// growthProgressReinitialize.push_back(-0.6098);
// growthProgressReinitialize.push_back(-3.3590);
// growthProgressReinitialize.push_back(0.5648);
// growthProgressReinitialize.push_back(-1.9684);
// growthProgressReinitialize.push_back(-0.5266);
// growthProgressReinitialize.push_back(-2.9763);
// growthProgressReinitialize.push_back(0.5215);
// growthProgressReinitialize.push_back(-1.8311);
// growthProgressReinitialize.push_back(-0.6464);
// growthProgressReinitialize.push_back(-3.2580);
// growthProgressReinitialize.push_back(0.6640);
// growthProgressReinitialize.push_back(-2.0480);

// //Ver10
// std::cout<<"Growth-progress profile Ver 10 is used"<<std::endl;
// growthProgressReinitialize.push_back( -1.9657);
// growthProgressReinitialize.push_back( -0.9376);
// growthProgressReinitialize.push_back( 0.3190);
// growthProgressReinitialize.push_back( -2.3226);
// growthProgressReinitialize.push_back( -0.9976);
// growthProgressReinitialize.push_back( 0.6312);
// growthProgressReinitialize.push_back( -2.4786);
// growthProgressReinitialize.push_back( -0.4733);
// growthProgressReinitialize.push_back( 0.6356);
// growthProgressReinitialize.push_back( -2.4616);
// growthProgressReinitialize.push_back( -1.3734);
// growthProgressReinitialize.push_back( -0.0702);
// growthProgressReinitialize.push_back( -2.4089);
// growthProgressReinitialize.push_back( -1.4015);
// growthProgressReinitialize.push_back( -0.0385);
// growthProgressReinitialize.push_back( -2.4974);
// growthProgressReinitialize.push_back( -1.2876);
// growthProgressReinitialize.push_back( 0.1130);
// growthProgressReinitialize.push_back( -2.4161);
// growthProgressReinitialize.push_back( -1.3892);
// growthProgressReinitialize.push_back( -0.0828);
// growthProgressReinitialize.push_back( -2.4240);
// growthProgressReinitialize.push_back( -1.1015);
// growthProgressReinitialize.push_back( 0.0612);
// growthProgressReinitialize.push_back( -2.2057);
// growthProgressReinitialize.push_back( -1.1241);
// growthProgressReinitialize.push_back( 0.7000);
// growthProgressReinitialize.push_back( -2.2458);
// growthProgressReinitialize.push_back( -1.2208);
// growthProgressReinitialize.push_back( -0.1874);
// growthProgressReinitialize.push_back( -2.2313);
// growthProgressReinitialize.push_back( -1.2195);
// growthProgressReinitialize.push_back( 0.1352);
// growthProgressReinitialize.push_back( -2.3796);
// growthProgressReinitialize.push_back( -0.9151);
// growthProgressReinitialize.push_back( 0.7106);
// growthProgressReinitialize.push_back( -1.9618);
// growthProgressReinitialize.push_back( -0.7738);
// growthProgressReinitialize.push_back( 0.4064);
// growthProgressReinitialize.push_back( -1.7924);
// growthProgressReinitialize.push_back( -0.6859);
// growthProgressReinitialize.push_back( 0.5891);
// growthProgressReinitialize.push_back( -1.9121);
// growthProgressReinitialize.push_back( -0.5693);
// growthProgressReinitialize.push_back( 0.7430);
// growthProgressReinitialize.push_back( -2.2218);
// growthProgressReinitialize.push_back( -1.0864);
// growthProgressReinitialize.push_back( 0.1956);
// growthProgressReinitialize.push_back( -2.1195);
// growthProgressReinitialize.push_back( -0.8659);
// growthProgressReinitialize.push_back( 0.2405);
// growthProgressReinitialize.push_back( -2.4792);
// growthProgressReinitialize.push_back( -1.1303);
// growthProgressReinitialize.push_back( 0.3634);
// growthProgressReinitialize.push_back( -2.2099);
// growthProgressReinitialize.push_back( -1.0004);
// growthProgressReinitialize.push_back( 0.2268);

// growthProgressReinitialize.push_back(   -9999.9);
// growthProgressReinitialize.push_back(   -9999.9);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0140);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0178);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0192);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0109);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0028);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0030);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0052);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0168);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0051);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0163);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0049);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0186);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0070);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0039);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0050);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0123);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0095);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0070);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0166);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0117);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0110);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0183);
// 		growthProgressReinitialize.push_back(-9999.9);//0.0057);
// 		for (int w = 0; w < allocPara_m.maxCellCount; w++){
// 			if (w < allocPara_m.currentActiveCellCount){
// 				cellInfoVecs.growthProgress[w] = growthProgressReinitialize[w];
// 			}
// 			else{
// 				cellInfoVecs.growthProgress[w] = 0.0;
// 			}
// 		}

// //Ver11
// std::cout<<"Growth-progress profile Ver 11 is used"<<std::endl;
// growthProgressReinitialize.push_back(	 0.8821);
// growthProgressReinitialize.push_back(   -0.5595);
// growthProgressReinitialize.push_back(   -0.0401);
// growthProgressReinitialize.push_back(    0.4962);
// growthProgressReinitialize.push_back(   -0.6362);
// growthProgressReinitialize.push_back(   -0.0459);
// growthProgressReinitialize.push_back(    0.6200);
// growthProgressReinitialize.push_back(   -0.6327);
// growthProgressReinitialize.push_back(   -0.0860);
// growthProgressReinitialize.push_back(    0.6928);
// growthProgressReinitialize.push_back(   -0.7430);
// growthProgressReinitialize.push_back(   -0.1503);
// growthProgressReinitialize.push_back(    0.4501);
// growthProgressReinitialize.push_back(   -0.6755);
// growthProgressReinitialize.push_back(   -0.1654);
// growthProgressReinitialize.push_back(    0.6631);
// growthProgressReinitialize.push_back(   -0.7269);
// growthProgressReinitialize.push_back(   -0.0736);
// growthProgressReinitialize.push_back(    0.6441);
// growthProgressReinitialize.push_back(   -0.5833);
// growthProgressReinitialize.push_back(    0.1423);
// growthProgressReinitialize.push_back(    0.6433);
// growthProgressReinitialize.push_back(   -0.3941);
// growthProgressReinitialize.push_back(    0.1151);
// growthProgressReinitialize.push_back(    0.7180);
// growthProgressReinitialize.push_back(   -0.6606);
// growthProgressReinitialize.push_back(    0.0325);
// growthProgressReinitialize.push_back(    0.8795);
// growthProgressReinitialize.push_back(   -0.6248);
// growthProgressReinitialize.push_back(    0.0403);
// growthProgressReinitialize.push_back(    0.7261);
// growthProgressReinitialize.push_back(   -0.4992);
// growthProgressReinitialize.push_back(    0.0718);
// growthProgressReinitialize.push_back(    0.6641);
// growthProgressReinitialize.push_back(   -0.4828);
// growthProgressReinitialize.push_back(    0.0918);
// growthProgressReinitialize.push_back(    0.7225);
// growthProgressReinitialize.push_back(   -0.5154);
// growthProgressReinitialize.push_back(    0.2171);
// growthProgressReinitialize.push_back(    0.7240);
// growthProgressReinitialize.push_back(   -0.5883);
// growthProgressReinitialize.push_back(    0.1415);
// growthProgressReinitialize.push_back(    0.7030);
// growthProgressReinitialize.push_back(   -0.5916);
// growthProgressReinitialize.push_back(    0.1664);
// growthProgressReinitialize.push_back(    0.6865);
// growthProgressReinitialize.push_back(   -0.3889);
// growthProgressReinitialize.push_back(    0.1533);
// growthProgressReinitialize.push_back(    0.6678);
// growthProgressReinitialize.push_back(   -0.4366);
// growthProgressReinitialize.push_back(    0.1610);
// growthProgressReinitialize.push_back(    0.7406);
// growthProgressReinitialize.push_back(   -0.4220);
// growthProgressReinitialize.push_back(    0.1918);
// growthProgressReinitialize.push_back(    0.7091);
// growthProgressReinitialize.push_back(   -0.4620);
// growthProgressReinitialize.push_back(    0.0381);
// growthProgressReinitialize.push_back(   -9999.9);
// growthProgressReinitialize.push_back(   -9999.9);

//Ver11
std::cout<<"Growth-progress profile Ver 11 is used"<<std::endl;
growthProgressReinitialize.push_back(	-9999.9);//Cell4
growthProgressReinitialize.push_back(   -0.5595);//Cell5
growthProgressReinitialize.push_back(   -0.0401);//Cell6
growthProgressReinitialize.push_back(    0.4962);//Cell7
growthProgressReinitialize.push_back(   -0.6362);//Cell8
growthProgressReinitialize.push_back(   -0.0459);//Cell9
growthProgressReinitialize.push_back(    0.6200);//Cell10
growthProgressReinitialize.push_back(   -0.6327);//Cell11
growthProgressReinitialize.push_back(   -0.0860);//Cell12
growthProgressReinitialize.push_back(    0.6928);//Cell13
growthProgressReinitialize.push_back(   -0.7430);//Cell14
growthProgressReinitialize.push_back(   -0.1503);//Cell15
growthProgressReinitialize.push_back(    0.4501);//Cell16
growthProgressReinitialize.push_back(   -0.6755);//Cell17
growthProgressReinitialize.push_back(   -0.1654);//Cell18
growthProgressReinitialize.push_back(    0.6631);//Cell19
growthProgressReinitialize.push_back(   -0.7269);//Cell20
growthProgressReinitialize.push_back(   -0.0736);//Cell21
growthProgressReinitialize.push_back(    0.6441);//Cell22
growthProgressReinitialize.push_back(   -0.5833);//Cell23
growthProgressReinitialize.push_back(    0.1423);//Cell24
growthProgressReinitialize.push_back(    0.6433);//Cell25
growthProgressReinitialize.push_back(   -0.3941);//Cell26
growthProgressReinitialize.push_back(    0.1151);//Cell27
growthProgressReinitialize.push_back(    0.7180);//Cell28
growthProgressReinitialize.push_back(   -0.6606);//Cell29
growthProgressReinitialize.push_back(    0.0325);//Cell30
growthProgressReinitialize.push_back(    0.8900);//Cell31
growthProgressReinitialize.push_back(   -0.6248);//Cell32
growthProgressReinitialize.push_back(    0.0403);//Cell33
growthProgressReinitialize.push_back(    0.7261);//Cell34
growthProgressReinitialize.push_back(   -0.4992);//Cell35
growthProgressReinitialize.push_back(    0.0718);//Cell36
growthProgressReinitialize.push_back(    0.6641);//Cell37
growthProgressReinitialize.push_back(   -0.4828);//Cell38
growthProgressReinitialize.push_back(    0.0918);//Cell39
growthProgressReinitialize.push_back(    0.7225);//Cell40
growthProgressReinitialize.push_back(   -0.5154);//Cell41
growthProgressReinitialize.push_back(    0.2171);//Cell42
growthProgressReinitialize.push_back(    0.7240);//Cell43
growthProgressReinitialize.push_back(   -0.5883);//Cell44
growthProgressReinitialize.push_back(    0.1415);//Cell45
growthProgressReinitialize.push_back(    0.7030);//Cell46
growthProgressReinitialize.push_back(   -0.5916);//Cell47
growthProgressReinitialize.push_back(    0.1664);//Cell48
growthProgressReinitialize.push_back(    0.6865);//Cell49
growthProgressReinitialize.push_back(   -0.3889);//Cell50
growthProgressReinitialize.push_back(    0.1533);//Cell51
growthProgressReinitialize.push_back(    0.6678);//Cell52
growthProgressReinitialize.push_back(   -0.4366);//Cell53
growthProgressReinitialize.push_back(    0.1610);//Cell54
growthProgressReinitialize.push_back(    0.7406);//Cell55
growthProgressReinitialize.push_back(   -0.4220);//Cell56
growthProgressReinitialize.push_back(    0.1918);//Cell57
growthProgressReinitialize.push_back(    0.7091);//Cell58
growthProgressReinitialize.push_back(   -0.4620);//Cell59
growthProgressReinitialize.push_back(    0.0381);//Cell60
growthProgressReinitialize.push_back(   -9999.9);//Cell61
growthProgressReinitialize.push_back(   -9999.9);//Cell62
		growthProgressReinitialize.push_back(-9999.9);//0.0140);
		growthProgressReinitialize.push_back(-9999.9);//0.0178);
		growthProgressReinitialize.push_back(-9999.9);//0.0192);
		growthProgressReinitialize.push_back(-9999.9);//0.0109);
		growthProgressReinitialize.push_back(-9999.9);//0.0028);
		growthProgressReinitialize.push_back(-9999.9);//0.0030);
		growthProgressReinitialize.push_back(-9999.9);//0.0052);
		growthProgressReinitialize.push_back(-9999.9);//0.0168);
		growthProgressReinitialize.push_back(-9999.9);//0.0051);
		growthProgressReinitialize.push_back(-9999.9);//0.0163);
		growthProgressReinitialize.push_back(-9999.9);//0.0049);
		growthProgressReinitialize.push_back(-9999.9);//0.0186);
		growthProgressReinitialize.push_back(-9999.9);//0.0070);
		growthProgressReinitialize.push_back(-9999.9);//0.0039);
		growthProgressReinitialize.push_back(-9999.9);//0.0050);
		growthProgressReinitialize.push_back(-9999.9);//0.0123);
		growthProgressReinitialize.push_back(-9999.9);//0.0095);
		growthProgressReinitialize.push_back(-9999.9);//0.0070);
		growthProgressReinitialize.push_back(-9999.9);//0.0166);
		growthProgressReinitialize.push_back(-9999.9);//0.0117);
		growthProgressReinitialize.push_back(-9999.9);//0.0110);
		growthProgressReinitialize.push_back(-9999.9);//0.0183);
		growthProgressReinitialize.push_back(-9999.9);//0.0057);
		for (int w = 0; w < allocPara_m.maxCellCount; w++){
			if (w < allocPara_m.currentActiveCellCount){
				cellInfoVecs.growthProgress[w] = growthProgressReinitialize[w];
			}
			else{
				cellInfoVecs.growthProgress[w] = 0.0;
			}
		}

	// Adding actin level (number of contractile springs) information here
	std::vector<double> distFromNucleus_normalMaxVec;
	distFromNucleus_normalMaxVec.push_back(-9999.9999);//CellID0
	distFromNucleus_normalMaxVec.push_back(-9999.9999);//CellID1
	distFromNucleus_normalMaxVec.push_back(0.197835294);
	distFromNucleus_normalMaxVec.push_back(0.197882353);
	distFromNucleus_normalMaxVec.push_back(0.216423529);
	distFromNucleus_normalMaxVec.push_back(0.2248     );
	distFromNucleus_normalMaxVec.push_back(0.223623529);
	distFromNucleus_normalMaxVec.push_back(0.248752941);
	distFromNucleus_normalMaxVec.push_back(0.251717647);
	distFromNucleus_normalMaxVec.push_back(0.281788235);
	distFromNucleus_normalMaxVec.push_back(0.284329412);
	distFromNucleus_normalMaxVec.push_back(0.284094118);
	distFromNucleus_normalMaxVec.push_back(0.282776471);
	distFromNucleus_normalMaxVec.push_back(0.2968     );
	distFromNucleus_normalMaxVec.push_back(0.299294118);
	distFromNucleus_normalMaxVec.push_back(0.3        );
	distFromNucleus_normalMaxVec.push_back(0.292658824);
	distFromNucleus_normalMaxVec.push_back(0.292329412);
	distFromNucleus_normalMaxVec.push_back(0.284847059);
	distFromNucleus_normalMaxVec.push_back(0.292564706);
	distFromNucleus_normalMaxVec.push_back(0.274823529);
	distFromNucleus_normalMaxVec.push_back(0.242776471);
	distFromNucleus_normalMaxVec.push_back(0.247529412);
	distFromNucleus_normalMaxVec.push_back(0.254682353);
	distFromNucleus_normalMaxVec.push_back(0.224658824);
	distFromNucleus_normalMaxVec.push_back(0.234776471);
	distFromNucleus_normalMaxVec.push_back(0.263105882);
	distFromNucleus_normalMaxVec.push_back(0.236988235);
	distFromNucleus_normalMaxVec.push_back(0.216423529);
	distFromNucleus_normalMaxVec.push_back(0.240658824);
	distFromNucleus_normalMaxVec.push_back(0.196941176);
	distFromNucleus_normalMaxVec.push_back(0.2184     );
	distFromNucleus_normalMaxVec.push_back(0.232094118);
	distFromNucleus_normalMaxVec.push_back(0.238070588);
	distFromNucleus_normalMaxVec.push_back(0.248329412);
	distFromNucleus_normalMaxVec.push_back(0.244611765);
	distFromNucleus_normalMaxVec.push_back(0.230070588);
	distFromNucleus_normalMaxVec.push_back(0.221411765);
	distFromNucleus_normalMaxVec.push_back(0.2352     );
	distFromNucleus_normalMaxVec.push_back(0.198305882);
	distFromNucleus_normalMaxVec.push_back(0.201552941);
	distFromNucleus_normalMaxVec.push_back(0.2176     );
	distFromNucleus_normalMaxVec.push_back(0.210964706);
	distFromNucleus_normalMaxVec.push_back(0.233082353);
	distFromNucleus_normalMaxVec.push_back(0.228141176);
	distFromNucleus_normalMaxVec.push_back(0.236564706);
	distFromNucleus_normalMaxVec.push_back(0.244141176);
	distFromNucleus_normalMaxVec.push_back(0.2232     );
	distFromNucleus_normalMaxVec.push_back(0.229741176);
	distFromNucleus_normalMaxVec.push_back(0.220847059);
	distFromNucleus_normalMaxVec.push_back(0.214164706);
	distFromNucleus_normalMaxVec.push_back(0.212376471);
	distFromNucleus_normalMaxVec.push_back(0.218870588);
	distFromNucleus_normalMaxVec.push_back(0.245129412);
	distFromNucleus_normalMaxVec.push_back(0.230964706);
	distFromNucleus_normalMaxVec.push_back(0.252047059);
	distFromNucleus_normalMaxVec.push_back(0.243576471);
	distFromNucleus_normalMaxVec.push_back(0.254870588);
	distFromNucleus_normalMaxVec.push_back(0.245976471);
	distFromNucleus_normalMaxVec.push_back(0.209270588);
	distFromNucleus_normalMaxVec.push_back(0.195529412);
	distFromNucleus_normalMaxVec.push_back(0.194588235);
	distFromNucleus_normalMaxVec.push_back(0.215482353);//CellID62

	std::vector<double> distFromNucleus_normalMaxVec_apical;
	distFromNucleus_normalMaxVec_apical.push_back(-9999.9999);//CellID0
	distFromNucleus_normalMaxVec_apical.push_back(-9999.9999);//CellID1
	distFromNucleus_normalMaxVec_apical.push_back(0.260141176);
	distFromNucleus_normalMaxVec_apical.push_back(0.255152941);
	distFromNucleus_normalMaxVec_apical.push_back(0.236658824);
	distFromNucleus_normalMaxVec_apical.push_back(0.236470588);
	distFromNucleus_normalMaxVec_apical.push_back(0.242776471);
	distFromNucleus_normalMaxVec_apical.push_back(0.2408     );
	distFromNucleus_normalMaxVec_apical.push_back(0.245835294);
	distFromNucleus_normalMaxVec_apical.push_back(0.250211765);
	distFromNucleus_normalMaxVec_apical.push_back(0.257741176);
	distFromNucleus_normalMaxVec_apical.push_back(0.248423529);
	distFromNucleus_normalMaxVec_apical.push_back(0.251058824);
	distFromNucleus_normalMaxVec_apical.push_back(0.247623529);
	distFromNucleus_normalMaxVec_apical.push_back(0.279058824);
	distFromNucleus_normalMaxVec_apical.push_back(0.276235294);
	distFromNucleus_normalMaxVec_apical.push_back(0.261505882);
	distFromNucleus_normalMaxVec_apical.push_back(0.266117647);
	distFromNucleus_normalMaxVec_apical.push_back(0.245788235);
	distFromNucleus_normalMaxVec_apical.push_back(0.249270588);
	distFromNucleus_normalMaxVec_apical.push_back(0.253082353);
	distFromNucleus_normalMaxVec_apical.push_back(0.253694118);
	distFromNucleus_normalMaxVec_apical.push_back(0.277223529);
	distFromNucleus_normalMaxVec_apical.push_back(0.261270588);
	distFromNucleus_normalMaxVec_apical.push_back(0.252376471);
	distFromNucleus_normalMaxVec_apical.push_back(0.259341176);
	distFromNucleus_normalMaxVec_apical.push_back(0.245223529);
	distFromNucleus_normalMaxVec_apical.push_back(0.251435294);
	distFromNucleus_normalMaxVec_apical.push_back(0.262635294);
	distFromNucleus_normalMaxVec_apical.push_back(0.297552941);
	distFromNucleus_normalMaxVec_apical.push_back(0.276564706);
	distFromNucleus_normalMaxVec_apical.push_back(0.283529412);
	distFromNucleus_normalMaxVec_apical.push_back(0.289552941);
	distFromNucleus_normalMaxVec_apical.push_back(0.272141176);
	distFromNucleus_normalMaxVec_apical.push_back(0.2424     );
	distFromNucleus_normalMaxVec_apical.push_back(0.243623529);
	distFromNucleus_normalMaxVec_apical.push_back(0.188235294);
	distFromNucleus_normalMaxVec_apical.push_back(0.257835294);
	distFromNucleus_normalMaxVec_apical.push_back(0.245223529);
	distFromNucleus_normalMaxVec_apical.push_back(0.226682353);
	distFromNucleus_normalMaxVec_apical.push_back(0.208282353);
	distFromNucleus_normalMaxVec_apical.push_back(0.211388235);
	distFromNucleus_normalMaxVec_apical.push_back(0.192470588);
	distFromNucleus_normalMaxVec_apical.push_back(0.1944     );
	distFromNucleus_normalMaxVec_apical.push_back(0.194117647);
	distFromNucleus_normalMaxVec_apical.push_back(0.227623529);
	distFromNucleus_normalMaxVec_apical.push_back(0.224894118);
	distFromNucleus_normalMaxVec_apical.push_back(0.242917647);
	distFromNucleus_normalMaxVec_apical.push_back(0.215152941);
	distFromNucleus_normalMaxVec_apical.push_back(0.208423529);
	distFromNucleus_normalMaxVec_apical.push_back(0.208941176);
	distFromNucleus_normalMaxVec_apical.push_back(0.195717647);
	distFromNucleus_normalMaxVec_apical.push_back(0.204094118);
	distFromNucleus_normalMaxVec_apical.push_back(0.183764706);
	distFromNucleus_normalMaxVec_apical.push_back(0.217223529);
	distFromNucleus_normalMaxVec_apical.push_back(0.222870588);
	distFromNucleus_normalMaxVec_apical.push_back(0.226917647);
	distFromNucleus_normalMaxVec_apical.push_back(0.213223529);
	distFromNucleus_normalMaxVec_apical.push_back(0.221458824);
	distFromNucleus_normalMaxVec_apical.push_back(0.207011765);
	distFromNucleus_normalMaxVec_apical.push_back(0.201035294);
	distFromNucleus_normalMaxVec_apical.push_back(0.191858824);
	distFromNucleus_normalMaxVec_apical.push_back(0.191576471);//CellID62



		for (int cellRank = 0; cellRank < allocPara_m.maxCellCount; cellRank++){
			if (cellRank < allocPara_m.currentActiveCellCount){
				// cellInfoVecs.distFromNucleus_normal[cellRank] = distFromNucleus_normalMax;//(distFromNucleus_normalMax) - (distFromNucleus_normalMax - (-14.0))*(1.0 - pow(distFromNucleus_normalVec[cellRank], 3.0));
				if (cellRank >= 2 && cellRank <= 21){
					cellInfoVecs.distFromNucleus_normal[cellRank] = 1.0*distFromNucleus_normalMaxVec[cellRank];
					cellInfoVecs.distFromNucleus_normal_apical[cellRank] = 1.0*distFromNucleus_normalMaxVec_apical[cellRank];
					if (cellRank == 2){
						std::cout<<"cellInfoVecs.distFromNucleus_normal["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal[cellRank]<<std::endl;
						std::cout<<"cellInfoVecs.distFromNucleus_normal_apical["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal_apical[cellRank]<<std::endl;
					}
				}
				else if (cellRank >= 22 && cellRank <= 42){
					cellInfoVecs.distFromNucleus_normal[cellRank] = (1.0)*distFromNucleus_normalMaxVec[cellRank];
					cellInfoVecs.distFromNucleus_normal_apical[cellRank] = 1.0*distFromNucleus_normalMaxVec_apical[cellRank];
					if (cellRank == 22){
						std::cout<<"cellInfoVecs.distFromNucleus_normal["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal[cellRank]<<std::endl;
						std::cout<<"cellInfoVecs.distFromNucleus_normal_apical["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal_apical[cellRank]<<std::endl;
					}
				}
				else if (cellRank >= 43 && cellRank <= 62){
					cellInfoVecs.distFromNucleus_normal[cellRank] = 1.0*distFromNucleus_normalMaxVec[cellRank];
					cellInfoVecs.distFromNucleus_normal_apical[cellRank] = 1.0*distFromNucleus_normalMaxVec_apical[cellRank];
					if (cellRank == 43){
						std::cout<<"cellInfoVecs.distFromNucleus_normal["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal[cellRank]<<std::endl;
						std::cout<<"cellInfoVecs.distFromNucleus_normal_apical["<<cellRank<<"] = "<<cellInfoVecs.distFromNucleus_normal_apical[cellRank]<<std::endl;
					}
				}
				else{
					cellInfoVecs.distFromNucleus_normal[cellRank] = -99999.9;
					cellInfoVecs.distFromNucleus_normal_apical[cellRank] = 99999.9;
					
					// std::cout<<"infoVecs.contractActomyo_multip["<<i<<"] = "<<infoVecs.contractActomyo_multip[i]<<std::endl;
				}
			}
			else{
				cellInfoVecs.distFromNucleus_normal[cellRank] = -99999.9;
				cellInfoVecs.distFromNucleus_normal_apical[cellRank] = 99999.9;
			}
		}
	}
	curTime = curTime + dt;
	bool tmpIsInitPhase= nodes->isInitPhase ;
	// std::cout << "     *** 3 ***" << endl;

 	if (nodes->isMemNodeTypeAssigned==false) {
    	assignMemNodeType();  // Ali
		cout << " I assigned boolean values for membrane node types " << endl; 
		nodes->isMemNodeTypeAssigned=true ; 
	}
	if (nodes->isMemNodeTypeAssigned_postCellDivision==false){
		cout<<" I begin to reassign boolean values for membrane node types post-growth"<<endl;
		assignMemNodeType();
		cout<<" I reassigned boolean values for membrane node types post-growth"<<endl;
		nodes->isMemNodeTypeAssigned_postCellDivision=true;
	}
	if (nodes->isMemNodeTypeAssigned_postAddNode==false){
		cout<<" I begin to reassign boolean values for membrane node types post-growth"<<endl;
		assignMemNodeType();
		cout<<" I reassigned boolean values for membrane node types post-growth"<<endl;
		nodes->isMemNodeTypeAssigned_postAddNode=true;
	}
	// std::cout << "     *** 4 ***" << endl;
#ifdef debugModeECM
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif

	double contractileSpringGrowthProgressSpeed = 0.5*9.3405e-5;//9.3405e-5;//For minM14max4 //5.7081e-5;//For minM7max4    //7.2649e-05;//For minM7max7
	double cellAreaGrowthProgressSpeed = 0.00022;//0.00011;//0.5*9.984e-5;//9.984e-5;//0.00002088;//0.02*0.000075/(4.0);
	double cellAreaGrowthProgressNonMitoticSpeed = 0.5*9.984e-5;//9.984e-5;//0.00002088;//0.02*0.000075/(4.0);
	double distFromBasalLocSpeed = 0.0001;//0.00005;
	double distFromApicalLocSpeed = 0.0001;//0.00005;
	if (timeRatio == 0){
		std::cout<<"contractileSpringGrowthProgressSpeed = "<<contractileSpringGrowthProgressSpeed<<std::endl;
		std::cout<<"cellAreaGrowthProgressSpeed = "<<cellAreaGrowthProgressSpeed<<std::endl;
		std::cout<<"cellAreaGrowthProgressNonMitoticSpeed = "<<cellAreaGrowthProgressNonMitoticSpeed<<std::endl;
		std::cout<<"Calculating the number of apical node and basal node in each cell"<<std::endl;
		for (int i = 0; i < allocPara_m.maxCellCount*allocPara_m.maxAllNodePerCell; i++){
			int cellRank = i/allocPara_m.maxAllNodePerCell;
			if (nodes->getInfoVecs().memNodeType1[i] == apical1){
				cellInfoVecs.numApicalVec[cellRank] += 1;
			}
			if (nodes->getInfoVecs().memNodeType1[i] == basal1){
				cellInfoVecs.numBasalVec[cellRank] += 1;
			}
		}
		for (int i = 0; i < allocPara_m.maxCellCount; i++){
			cellInfoVecs.cellRankVec[i] = i;
			std::cout<<"Cell["<<i<<"] has "<<cellInfoVecs.numApicalVec[i]<<" apical nodes and "<<cellInfoVecs.numBasalVec[i]<<" basal nodes"<<std::endl;
		}
	}


	for (int w = 0; w < allocPara_m.currentActiveCellCount; w++){
		if (cellInfoVecs.isEnteringMitotic[w] == true){
			cellInfoVecs.contractileSpringGrowthProgress[w] += contractileSpringGrowthProgressSpeed;//0.00015;
			cellInfoVecs.cellAreaGrowthProgress[w] += cellAreaGrowthProgressSpeed; //0.00005;
			cellInfoVecs.cellAreaGrowthProgressNonMitotic[w] += cellAreaGrowthProgressNonMitoticSpeed;//0.00005;
		}
		else{
			cellInfoVecs.contractileSpringGrowthProgress[w] = 0.0;
			cellInfoVecs.cellAreaGrowthProgressNonMitotic[w] += cellAreaGrowthProgressNonMitoticSpeed;//0.00005;
			cellInfoVecs.cellAreaGrowthProgress[w] += cellAreaGrowthProgressSpeed;//0.00005;
			cellInfoVecs.distFromBasalLoc[w] += distFromBasalLocSpeed;
			cellInfoVecs.distFromApicalLoc[w] += distFromApicalLocSpeed;
			// if (w == 86){
			// 	std::cout<<"cellInfoVecs.distFromNucleus_normal[w] : " <<cellInfoVecs.distFromNucleus_normal[w]<<std::endl;
			// 	std::cout<<"cellInfoVecs.distFromNucleus_normal_apical[w] : "<<cellInfoVecs.distFromNucleus_normal_apical[w]<<std::endl;
			// 	std::cout<<"distFromNucleus_normal[w]*_individualCellHeight[w] : "<<cellInfoVecs.distFromNucleus_normal[w]*cellInfoVecs.individualCellHeight[w]<<std::endl;
			// 	std::cout<<"distFromNucleus_normal_apical[w]*_individualCellHeight[w] : "<<cellInfoVecs.distFromNucleus_normal_apical[w]*cellInfoVecs.individualCellHeight[w]<<std::endl;
			// }
		}
	}
	
	 
	
	
	// std::cout<< "	*** 4.5 ***"	<< std::endl;
	int nucRepopuRate = 5000;//2000;
	// if (timeRatio > timeRatio_Crit_Division && relaxCount % nucRepopuRate == 0){
	if (relaxCount % nucRepopuRate == 0){
		// if (relaxCount == 100){
		// 	std::cout<<"Nucleus are being repopuluated every "<<nucRepopuRate<<" time stpes"<<std::endl;
		// }
		
		uint targetCellRank;
		for (int i = 0; i < allocPara_m.currentActiveCellCount; i++){
			if (cellInfoVecs.isEnteringMitotic[i] == false){
				continue;
			}
			targetCellRank = i;
			if (targetCellRank < 2){
				continue;
			}
			else if (targetCellRank > 62 && targetCellRank < 86){
				continue;
			}
			else if (cellInfoVecs.activeIntnlNodeCounts[targetCellRank] >= allocPara_m.maxIntnlNodePerCell){
				continue;
			}

			int startIndex = targetCellRank*allocPara_m.maxAllNodePerCell + allocPara_m.maxMembrNodePerCell;
			int currentNumOfNucleus = 0;
			for (int j = startIndex; j < (targetCellRank+1)*allocPara_m.maxAllNodePerCell; j++){
				if (nodes->getInfoVecs().nodeIsActive[j] == true){
					currentNumOfNucleus+=1;
				}
			}
			int numOfMissingNulceus = allocPara_m.maxIntnlNodePerCell - currentNumOfNucleus;
			// std::cout<<"Currently missing "<<numOfMissingNulceus<<" nucleus particle in cell "<<targetCellRank<<std::endl;

			double total_nucNodeLocX = 0.0;
			double total_nucNodeLocY = 0.0;
			for (int k = 0; k < currentNumOfNucleus; k++){
				total_nucNodeLocX += nodes->getInfoVecs().nodeLocX[startIndex + k];
				total_nucNodeLocY += nodes->getInfoVecs().nodeLocY[startIndex + k];
			}
			total_nucNodeLocX = total_nucNodeLocX/currentNumOfNucleus;
			total_nucNodeLocY = total_nucNodeLocY/currentNumOfNucleus;

			nodes->getInfoVecs().nodeLocX[startIndex + currentNumOfNucleus] = 
								total_nucNodeLocX;
			nodes->getInfoVecs().nodeLocY[startIndex + currentNumOfNucleus] = 
								total_nucNodeLocY;
			nodes->getInfoVecs().nodeIsActive[startIndex + currentNumOfNucleus] = true;

			// if (numOfMissingNulceus != 0){
			// 	for (int k = 0; k < numOfMissingNulceus; k++){
			// 		nodes->getInfoVecs().nodeLocX[startIndex + currentNumOfNucleus + k] = 
			// 							nodes->getInfoVecs().nodeLocX[startIndex + k] + 0.01;
			// 		nodes->getInfoVecs().nodeLocY[startIndex + currentNumOfNucleus + k] = 
			// 							nodes->getInfoVecs().nodeLocY[startIndex + k] + 0.01;
			// 		nodes->getInfoVecs().nodeIsActive[startIndex + currentNumOfNucleus + k] = true;
			// 	}
			// }
			currentNumOfNucleus = 0;
			for (int j = startIndex; j < (targetCellRank+1)*allocPara_m.maxAllNodePerCell; j++){
				if (nodes->getInfoVecs().nodeIsActive[j] == true){
					currentNumOfNucleus+=1;
				}
				cellInfoVecs.activeIntnlNodeCounts[targetCellRank] = currentNumOfNucleus;
			}
			// std::cout<<"number of nucleus nodes in cell["<<targetCellRank<<"] = "<<currentNumOfNucleus<<" after repopulation"<<std::endl;
			// std::cout<<"CellRank = "<<targetCellRank<<"new activeIntnlNodeCounts = "<<cellInfoVecs.activeIntnlNodeCounts[targetCellRank]<<std::endl;
		}
	}
	// std::cout << "     *** 4.75 ***" << endl;
	
	
    // computeApicalLoc(timeRatio, timeRatio_Crit_Division);  //Ali
	// //////
	
	// // std::cout << "     *** 5 ***" << endl;
    // computeBasalLoc();  //Ali
	// std::cout << "     *** 6 ***" << endl;

	// for (int i = 0; i < allocPara_m.currentActiveCellCount; i++){
	// 	cellInfoVecs.individualCellHeight[i] = sqrt((cellInfoVecs.apicalLocX[i] - cellInfoVecs.basalLocX[i])*(cellInfoVecs.apicalLocX[i] - cellInfoVecs.basalLocX[i]) +
	// 													(cellInfoVecs.apicalLocY[i] - cellInfoVecs.basalLocY[i])*(cellInfoVecs.apicalLocY[i] - cellInfoVecs.basalLocY[i]));
	// 	cellInfoVecs.individualCellHeight_springHeight[i] = cellInfoVecs.individualCellHeight[i]*0.3; 
	// 	//Note: what is calculated here really? Well, we want to extend a distance away from the cell center (-ish) position to determine how far
	// 	//      contractile spring will be placed. We now assume that the upper and lower 20% of the cell height is covered by contractile spring,
	// 	//      therefore, leading to the non-spring portion taking 60% (hence 30% above center and 30% below center).
	// 	cellInfoVecs.distFromNucleus_normal[i] = -1.0*cellInfoVecs.individualCellHeight_springHeight[i];
	// 	cellInfoVecs.distFromNucleus_normal_apical[i] = 1.0*cellInfoVecs.individualCellHeight_springHeight[i];
	// }

	
		computeApicalLoc(timeRatio, timeRatio_Crit_Division);
    	computeBasalLoc(); 
	uint recalculateCellHeight = 4;
	if (timeRatio == 0){// || relaxCount%recalculateCellHeight == 0){	
		computeIndividualCellHeight(distFromNucleus_normalMax1,distFromNucleus_normalMax2,distFromNucleus_normalMax3,
							 distFromNucleus_normalMax_apical1,distFromNucleus_normalMax_apical2,distFromNucleus_normalMax_apical3);//This function does the above commented out computation.
		// computeIndividualCellHeight_Ver2();//Ver2 does cell height calculation only, not manipulating the number of contractile springs.
	}
	else if (relaxCount%recalculateCellHeight == 0){
		computeIndividualCellHeight_Ver2();//Ver2 does cell height calculation only, not manipulating the number of contractile springs.
	}
#ifdef debugModeECM
	cudaEventRecord(start3, 0);
	cudaEventSynchronize(start3);
	cudaEventElapsedTime(&elapsedTime2, start2, start3);
#endif

	
	// eCMCellInteraction(cellPolar,subCellPolar,tmpIsInitPhase, timeRatio, timeRatio_Crit_ECM, timeRatio_Crit_Division, relaxCount);
	eCMCellInteraction(cellPolar,subCellPolar,tmpIsInitPhase, timeRatio, timeRatio_Crit_ECM, timeRatio_Crit_Division, relaxCount, mitoticThreshold);

	if (cellInfoVecs.isOneTimeStepPostDivision != false){
		std::cout<<"Condition ' cellInfoVecs.isOneTimeStepPostDivision != false ' is satisfied"<<std::endl;
		std::cout<<"cellInfoVecs.isOneTimeStepPostDivision = "<<cellInfoVecs.isOneTimeStepPostDivision<<std::endl;
		std::cout<<"cellInfoVecs.isPostDivision = "<<cellInfoVecs.isPostDivision<<std::endl;
		/*for (int kk = 0; kk < cellInfoVecs.eCMNeighborId.size(); kk++){
			std::cout<<"eCMNeighborId["<<kk<<"] after division= "<<cellInfoVecs.eCMNeighborId[kk]<<std::endl;
		}*/
	}

	if (cellInfoVecs.isPostDivision == true){
		cellInfoVecs.isOneTimeStepPostDivision = true;
		cellInfoVecs.isPostDivision = false;
		std::cout<<"Condition ' cellInfoVecs.isPostDivision == true ' is satisfied"<<std::endl;
		std::cout<<"cellInfoVecs.isOneTimeStepPostDivision = "<<cellInfoVecs.isOneTimeStepPostDivision<<std::endl;
		std::cout<<"cellInfoVecs.isPostDivision = "<<cellInfoVecs.isPostDivision<<std::endl;
	}
	else if (cellInfoVecs.isPostDivision == false && cellInfoVecs.isOneTimeStepPostDivision == true){
		cellInfoVecs.isTwoTimeStepPostDivision = false;//true;
		cellInfoVecs.isOneTimeStepPostDivision = false;
	}
	else{
		cellInfoVecs.isOneTimeStepPostDivision = false;
		cellInfoVecs.isTwoTimeStepPostDivision = false;
	}

	if (cellInfoVecs.isPostAddMembrNodes == true){
		cellInfoVecs.isPostAddMembrNodes = false;
	}
	
	// if (timeRatio > timeRatio_Crit_Division && nodes->isECMPairPrinted==false){
	// if (timeRatio == timeRatio_Crit_Division){
	// 	// std::cout<<"bdrynodecount = "<<allocPara_m.bdryNodeCount<<std::endl;

	// 	// for (int kk = 0; kk < cellInfoVecs.basalLocX.size(); kk++){
	// 	// 	// std::cout<<"BasalLocX["<<kk<<"] = "<< cellInfoVecs.basalLocX[kk]<<std::endl;
	// 	// 	// std::cout<<"BasalLocY["<<kk<<"] = "<< cellInfoVecs.basalLocY[kk]<<std::endl;
	// 	// }

	// 	for (int kk = 0; kk < (allocPara_m.currentActiveCellCount*allocPara_m.maxAllNodePerCell); kk++){
	// 		if (kk/allocPara_m.maxAllNodePerCell ==86 || kk/allocPara_m.maxAllNodePerCell == 31 || kk/allocPara_m.maxAllNodePerCell == 32 ||kk/allocPara_m.maxAllNodePerCell == 30){
	// 		}
	// 		else{
	// 			continue;
	// 		}
			
	// 	}
	// 	for (int kk = 0; kk < cellInfoVecs.eCMNeighborId.size(); kk++){
	// 		std::cout<<"eCMNeighborId["<<kk<<"] = "<<cellInfoVecs.eCMNeighborId[kk]<<std::endl;
	// 	}
	// 	nodes->isECMPairPrinted=true;
	// }
	// std::cout << "     *** 7 ***" << endl;

#ifdef debugModeECM
	cudaEventRecord(start4, 0);
	cudaEventSynchronize(start4);
	cudaEventElapsedTime(&elapsedTime3, start3, start4);
#endif


	computeCenterPos_M2(); //Ali 
	// std::cout << "     *** 8 ***" << endl;
	computeInternalAvgPos_M(); //Ali // right now internal points represent nucleus
	// std::cout << "     *** 9 ***" << endl;
	//computeNucleusLoc() ;

#ifdef debugModeECM
	cudaEventRecord(start5, 0);
	cudaEventSynchronize(start5);
	cudaEventElapsedTime(&elapsedTime4, start4, start5);
#endif

 	if (isInitNucPercentCalculated==false && controlPara.resumeSimulation==0) {
		computeNucleusIniLocPercent(); //Ali 
		writeNucleusIniLocPercent(); //Ali 
		isInitNucPercentCalculated=true ; 
		cout << " I computed initial location of nucleus positions in percent" << endl; 
	}
	else if (isInitNucPercentCalculated==false && controlPara.resumeSimulation==1){
		readNucleusIniLocPercent(); //Ali 
		isInitNucPercentCalculated=true ; 
		for (int j = 0; j < cellInfoVecs.nucleusLocPercent.size(); j++){
			std::cout<<"nucleusLocPercent["<<j<<"] = "<<cellInfoVecs.nucleusLocPercent[j]<<std::endl;
		}
		cout << " I read initial location of nucleus positions in percent, since I am in resume mode" << endl; 
	}
	// std::cout << "     *** 10 ***" << endl;

	// if (timeRatio == 0){
		computeNucleusDesireLoc() ; // Ali
	// }
	// if (timeRatio == timeRatio_Crit_Division || timeRatio == timeRatio_Crit_Division+0.2){
	// 	std::cout<<"timeRatio = "<<timeRatio<<std::endl;
	// 	std::cout<<"cellInfoVecs.nucDesireDistApical[31] = "<<cellInfoVecs.nucDesireDistApical[31]<<std::endl;
	// 	std::cout<<"cellInfoVecs.centerCoordX[31] = "<<cellInfoVecs.centerCoordX[31]<<std::endl;
	// 	std::cout<<"cellInfoVecs.centerCoordY[31] = "<<cellInfoVecs.centerCoordY[31]<<std::endl;
	// 	std::cout<<"cellInfoVecs.apicalLocX[31] = "<<cellInfoVecs.apicalLocX[31]<<std::endl;
	// 	std::cout<<"cellInfoVecs.apicalLocY[31] = "<<cellInfoVecs.apicalLocY[31]<<std::endl;
	// 	std::cout<<"cellInfoVecs.nucDesireDistApical[87] = "<<cellInfoVecs.nucDesireDistApical[87]<<std::endl;
	// 	std::cout<<"cellInfoVecs.centerCoordX[87] = "<<cellInfoVecs.centerCoordX[87]<<std::endl;
	// 	std::cout<<"cellInfoVecs.centerCoordY[87] = "<<cellInfoVecs.centerCoordY[87]<<std::endl;
	// 	std::cout<<"cellInfoVecs.apicalLocX[87] = "<<cellInfoVecs.apicalLocX[87]<<std::endl;
	// 	std::cout<<"cellInfoVecs.apicalLocY[87] = "<<cellInfoVecs.apicalLocY[87]<<std::endl;
	// }
	// std::cout << "     *** 11 ***" << endl;

#ifdef debugModeECM
	cudaEventRecord(start6, 0);
	cudaEventSynchronize(start6);
	cudaEventElapsedTime(&elapsedTime5, start5, start6);
#endif


//	if (tmpIsInitPhase==false) {
//		updateInternalAvgPosByNucleusLoc_M ();
//	}
	//PlotNucleus (lastPrintNucleus, outputFrameNucleus) ;  
    //BC_Imp_M() ;  //Ali

	applySceCellDisc_M();
	// std::cout << "     *** 12 ***" << endl;

#ifdef debugModeECM
	cudaEventRecord(start7, 0);
	cudaEventSynchronize(start7);
	cudaEventElapsedTime(&elapsedTime6, start6, start7);
#endif
	if (isBasalActinPresent) {
		// cout << " I am applying basal contraction" << endl ; 
		// applyMembContraction(timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_Division, distFromNucleus_max, distFromNucleus_min, percentage_before_timeRatio_Crit_Division_scaling) ;  // Ali
		applyMembContraction2(timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_Division, distFromNucleus_max, distFromNucleus_min,mitoRndActomyoStrengthScaling);// percentage_before_timeRatio_Crit_Division_scaling) ;  // Kevin
	}
	// std::cout << "     *** 13 ***" << endl;
#ifdef debugModeECM
	cudaEventRecord(start8, 0);
	cudaEventSynchronize(start8);
	cudaEventElapsedTime(&elapsedTime7, start7, start8);
#endif


	//	applyNucleusEffect() ;
	//	applyForceInteractionNucleusAsPoint() ; 

	applyMemForce_M(cellPolar,subCellPolar);
	// if (timeRatio > timeRatio_Crit_Division && nodes->isActinLevelDisplayed == false){
	// 	for (int w = 0; w < nodes->getInfoVecs().nodeActinLevel.size(); w++){
	// 		if (w/allocPara_m.maxAllNodePerCell == 31 || w/allocPara_m.maxAllNodePerCell == 86){
	// 			std::cout<<"actin level["<<w/allocPara_m.maxAllNodePerCell<<"] = "<<nodes->getInfoVecs().nodeActinLevel[w]<<std::endl;
	// 		}
	// 	}
	// 	nodes->isActinLevelDisplayed = true;
	// }
	// std::cout << "     *** 14 ***" << endl;

#ifdef debugModeECM
	cudaEventRecord(start9, 0);
	cudaEventSynchronize(start9);
	cudaEventElapsedTime(&elapsedTime8, start8, start9);
#endif



	applyVolumeConstraint(timeRatio, timeRatio_Crit_Division, volume_Increase_Target_Ratio, volume_Increase_Scale, postDivision_restorationRateScale, cycle);  //Ali 
	// std::cout << "     *** 15 ***" << endl;
	if (timeRatio > timeRatio_Crit_Division && nodes->isActinLevelDisplayed == false){
		for (int kk = 0; kk < cellInfoVecs.cellAreaVec.size(); kk++){
			// if (kk == 31 || kk == 86){
			// 	std::cout<<"Cell["<<kk<<"] area = "<<cellInfoVecs.cellAreaVec[kk]<<std::endl;
			// }

			// for (int kk = 0; kk < allocPara_m.maxCellCount; kk++){
			// 	double cellRank = kk;
			// 	uint intnlIndxMemBegin = cellRank * allocPara_m.maxAllNodePerCell;
			// 	uint intnlIndxBegin = cellRank * allocPara_m.maxAllNodePerCell + allocPara_m.maxMembrNodePerCell;
			// 	uint intnlIndxEnd = intnlIndxBegin + cellInfoVecs.activeIntnlNodeCounts[kk];
			// 	std::cout<<"cellRank = "<<cellRank<<", intnlIndxMemBegin = "<<intnlIndxMemBegin<<", intnlIndxBegin = "<<intnlIndxBegin<<", intnlIndxEnd = "<<intnlIndxEnd<<std::endl;
			// }
		}
		// for (int kk = 0; kk < cellInfoVecs.Cell_Damp.size(); kk++){
		// 	std::cout<<"Node Damping for cell ["<<kk<<"] = "<<cellInfoVecs.Cell_Damp[kk]<<std::endl;
		// }
		nodes->isActinLevelDisplayed=true;
	}

#ifdef debugModeECM
	cudaEventRecord(start10, 0);
	cudaEventSynchronize(start10);
	cudaEventElapsedTime(&elapsedTime9, start9, start10);
#endif



	//ApplyExtForces() ; // now for single cell stretching
	//computeContractileRingForces() ; 


//	computeCenterPos_M();    //Ali cmment //
	

 	// if (isCellGrowSet==false) {
		// growAtRandom_M(dt);
		// growAtRandom_M(growthProgressSpeed);
		growAtRandom_M_Ver2(growthProgressSpeed, mitoticThreshold);
		// cout << "I set the growth level. Since the cells are not growing a divising for this simulation I won't go inside this function any more" << endl ;
		// isCellGrowSet=true ;
	// }
	// std::cout << "     *** 16 ***" << endl;


    relaxCount=relaxCount+1 ; 
	// std::cout << "     *** 17 ***" << endl;
	distributeCellGrowthProgress_M();
	// std::cout << "     *** 18 ***" << endl;

    findTangentAndNormal_M();//AAMIRI ADDED May29
	// std::cout << "     *** 19 ***" << endl;
	StoreNodeOldPositions() ; 
	// std::cout << "     *** 20 ***" << endl;
#ifdef debugModeECM
	cudaEventRecord(start11, 0);
	cudaEventSynchronize(start11);
	cudaEventElapsedTime(&elapsedTime10, start10, start11);
#endif



	allComponentsMove_M();
	// std::cout << "     *** 21 ***" << endl;
	int checkForMitosisAndDivision = 250;
	int checkForOverextension = 500;
	bool membrAddingNode = false;

	// Uncomment the line below to turn cell proliferation off 
	// if (1 < 0){
	if (relaxCount % checkForMitosisAndDivision == 0){
		if (allocPara_m.currentActiveCellCount < allocPara_m.maxCellCount){
			// std::cout<<"contractileSpringGrowthProgress[10] = "<<cellInfoVecs.contractileSpringGrowthProgress[10]<<std::endl;
			enterMitoticCheckForDivAxisCal(mitoticThreshold) ;
			
			divide2D_M(volume_Increase_Target_Ratio, timeRatio, thresholdToIntroduceNewCell);
			// std::cout<<"	*** 4.25 *** "<<	std::endl;
			if (relaxCount % checkForOverextension == 0){
				updateMembrGrowthProgress_M(); 
				handleMembrGrowth_M(maxApicalBasalNodeNum, maxLengthToAddMemNodes);
			}
			nodes->adhUpdate=true;
		}
	}

#ifdef debugModeECM
	cudaEventRecord(start12, 0);
	cudaEventSynchronize(start12);
	cudaEventElapsedTime(&elapsedTime11, start11, start12);
#endif



    
	//allComponentsMoveImplicitPart() ;  
	// if (relaxCount%checkForMitosisAndDivision==0) {
	// 	updateMembrGrowthProgress_M();  
	// 	handleMembrGrowth_M();
	// 	// std::cout << "     *** 10 ***" << endl;
	// 	std::cout.flush();
	// 	// relaxCount=0 ; // Ali
	// 	nodes->adhUpdate=true; // Ali 
	// }

# ifdef debugModeECM 
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime12, start12, stop);
	std::cout << "time 1 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime1 << endl ; 
	std::cout << "time 2 spent in cell for moving the membrane node of cells and ECM nodes are: " << elapsedTime2 << endl ; 
	std::cout << "time 3 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime3 << endl ; 
	std::cout << "time 4 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime4 << endl ; 
	std::cout << "time 5 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime5 << endl ; 
	std::cout << "time 6 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime6 << endl ; 
	std::cout << "time 7 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime7 << endl ; 
	std::cout << "time 8 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime8 << endl ; 
	std::cout << "time 9 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime9 << endl ; 
	std::cout << "time 10 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime10 << endl ; 
	std::cout << "time 11 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime11 << endl ; 
	std::cout << "time 12 spent in cell module for moving the membrane node of cells and ECM nodes are: " << elapsedTime12 << endl ; 
#endif

// std::cout << "     *** 22 ***" << endl;
}

void SceCells::runStretchTest(double dt) {
	this->dt = dt;
	computeCenterPos();
	growAlongX(false, dt);
	moveNodes();
}

void SceCells::growAlongX(bool isAddPt, double d_t) {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	setGrowthDirXAxis();

//std::cout << "after copy grow info" << std::endl;
	updateGrowthProgress();
//std::cout << "after update growth progress" << std::endl;
	decideIsScheduleToGrow();
//std::cout << "after decode os schedule to grow" << std::endl;
	computeCellTargetLength();
//std::cout << "after compute cell target length" << std::endl;
	computeDistToCellCenter();
//std::cout << "after compute dist to center" << std::endl;
	findMinAndMaxDistToCenter();
//std::cout << "after find min and max dist" << std::endl;
	computeLenDiffExpCur();
//std::cout << "after compute diff " << std::endl;
	stretchCellGivenLenDiff();

	if (isAddPt) {
		addPointIfScheduledToGrow();
	}
}

void SceCells::growWithStress(double d_t) {
}

std::vector<CVector> SceCells::getAllCellCenters() {
	thrust::host_vector<double> centerX = cellInfoVecs.centerCoordX;
	thrust::host_vector<double> centerY = cellInfoVecs.centerCoordY;
	thrust::host_vector<double> centerZ = cellInfoVecs.centerCoordZ;
	std::vector<CVector> result;
	for (uint i = 0; i < allocPara.currentActiveCellCount; i++) {
		CVector pos = CVector(centerX[i], centerY[i], centerZ[i]);
		result.push_back(pos);
	}
	return result;
}

void SceCells::setGrowthDirXAxis() {
	thrust::fill(cellInfoVecs.growthXDir.begin(),
			cellInfoVecs.growthXDir.begin() + allocPara.currentActiveCellCount,
			1.0);
	thrust::fill(cellInfoVecs.growthYDir.begin(),
			cellInfoVecs.growthYDir.begin() + allocPara.currentActiveCellCount,
			0.0);
	thrust::fill(cellInfoVecs.growthSpeed.begin(),
			cellInfoVecs.growthSpeed.begin() + allocPara.currentActiveCellCount,
			growthAuxData.fixedGrowthSpeed);
}

std::vector<double> SceCells::getGrowthProgressVec() {
	thrust::host_vector<double> growthProVec = cellInfoVecs.growthProgress;
	std::vector<double> result;
	for (uint i = 0; i < allocPara.currentActiveCellCount; i++) {
		result.push_back(growthProVec[i]);
	}
	return result;
}

void SceCells::copyCellsPreDivision_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
			* allocPara_m.maxAllNodePerCell;

	divAuxData.tmpApicalLoc = thrust::device_vector<double>(2,0.0);
	divAuxData.tmpBasalLoc = thrust::device_vector<double>(2,0.0);

	divAuxData.tmpIsActive_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, true);
	divAuxData.tmpNodePosX_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodePosY_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodeType = thrust::device_vector<MembraneType1>(
			divAuxData.nodeStorageCount, notAssigned1); //Ali 
	divAuxData.tmpNodeMemMirrorIndex_M = thrust::device_vector<int>(
			divAuxData.nodeStorageCount, -1);

	divAuxData.tmpCellRank_M = thrust::device_vector<uint>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpNucleusCenterPosX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpNucleusCenterPosY_M = thrust::device_vector<double>(
		divAuxData.toBeDivideCount, 0);

	divAuxData.tmpIntAvgX_M = thrust::device_vector<double>(   //Ali 
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpIntAvgY_M = thrust::device_vector<double>(   //Ali 
			divAuxData.toBeDivideCount, 0);

	divAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodeType1 = thrust::device_vector<MembraneType1>(
			divAuxData.nodeStorageCount, notAssigned1); //Ali

	divAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodeType2 = thrust::device_vector<MembraneType1>(
			divAuxData.nodeStorageCount, notAssigned1); //Ali

        //A&A
        divAuxData.tmpHertwigXdir = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
        divAuxData.tmpHertwigYdir = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
        //A&A

// step 2 , continued // copy node info values ready for division /comment A&A
	thrust::counting_iterator<uint> iStart(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
									+ allocPara_m.bdryNodeCount))		
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(cellInfoVecs.isDividing.begin(),
					make_transform_iterator(iStart,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
							divAuxData.tmpNodePosX_M.begin(),
							divAuxData.tmpNodePosY_M.begin(),
							divAuxData.tmpNodeType.begin(),
							divAuxData.tmpNodeMemMirrorIndex_M.begin())), isTrue());
// step 3 , continued  //copy cell info values ready for division /comment A&A
	thrust::counting_iterator<uint> iBegin(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.HertwigXdir.begin(),
							cellInfoVecs.HertwigYdir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.HertwigXdir.begin(),
							cellInfoVecs.HertwigYdir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpCellRank_M.begin(),
							divAuxData.tmpDivDirX_M.begin(),
							divAuxData.tmpDivDirY_M.begin(),
                            divAuxData.tmpHertwigXdir.begin(),
                            divAuxData.tmpHertwigYdir.begin(),
							divAuxData.tmpCenterPosX_M.begin(),
							divAuxData.tmpCenterPosY_M.begin(),
							divAuxData.tmpIntAvgX_M.begin(),
							divAuxData.tmpIntAvgY_M.begin()
							)), isTrue());

			for (int w = 0; w < cellInfoVecs.isDividing.size(); w++){
				if (cellInfoVecs.isDividing[w] == true){
					divAuxData.tmpApicalLoc[0] = cellInfoVecs.apicalLocX[w];
					divAuxData.tmpApicalLoc[1] = cellInfoVecs.apicalLocY[w];
					divAuxData.tmpBasalLoc[0] = cellInfoVecs.basalLocX[w];
					divAuxData.tmpBasalLoc[1] = cellInfoVecs.basalLocY[w];
				}
			}
}

void SceCells::copyCellsForPrintingOnly_M() {
	uint totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	
	printAuxData.nodeStorageCount = printAuxData.toBeDivideCount
			* allocPara_m.maxAllNodePerCell;
// std::cout<<"copy 1"<<std::endl;
	printAuxData.tmpApicalLoc = thrust::device_vector<double>(2,0.0);
	printAuxData.tmpBasalLoc = thrust::device_vector<double>(2,0.0);
// std::cout<<"copy 2"<<std::endl;
	printAuxData.tmpIsActive_M = thrust::device_vector<bool>(
			printAuxData.nodeStorageCount, true);
	printAuxData.tmpNodePosX_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpNodePosY_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpNodeType = thrust::device_vector<MembraneType1>(
			printAuxData.nodeStorageCount, notAssigned1); //Ali 
	printAuxData.tmpNodeMemMirrorIndex_M = thrust::device_vector<int>(
			printAuxData.nodeStorageCount, -1);
// std::cout<<"copy 3"<<std::endl;	
	printAuxData.tmpCellRank_M = thrust::device_vector<uint>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpDivDirX_M = thrust::device_vector<double>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpDivDirY_M = thrust::device_vector<double>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpNucleusCenterPosX_M = thrust::device_vector<double>(
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpNucleusCenterPosY_M = thrust::device_vector<double>(
		printAuxData.toBeDivideCount, 0);
// std::cout<<"copy 4"<<std::endl;	
	printAuxData.tmpIntAvgX_M = thrust::device_vector<double>(   //Ali 
			printAuxData.toBeDivideCount, 0);
	printAuxData.tmpIntAvgY_M = thrust::device_vector<double>(   //Ali 
			printAuxData.toBeDivideCount, 0);
// std::cout<<"copy 5"<<std::endl;	
	printAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
			printAuxData.nodeStorageCount, false);
	printAuxData.tmpXPos1_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpYPos1_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpNodeType1 = thrust::device_vector<MembraneType1>(
			printAuxData.nodeStorageCount, notAssigned1); //Ali
// std::cout<<"copy 6"<<std::endl;	
	printAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
			printAuxData.nodeStorageCount, false);
	printAuxData.tmpXPos2_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpYPos2_M = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
	printAuxData.tmpNodeType2 = thrust::device_vector<MembraneType1>(
			printAuxData.nodeStorageCount, notAssigned1); //Ali
// std::cout<<"copy 7"<<std::endl;	
        //A&A
        printAuxData.tmpHertwigXdir = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
        printAuxData.tmpHertwigYdir = thrust::device_vector<double>(
			printAuxData.nodeStorageCount, 0.0);
        //A&A
// std::cout<<"copy 8"<<std::endl;	
// step 2 , continued // copy node info values ready for division /comment A&A
	thrust::counting_iterator<uint> iStart(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
							+ allocPara_m.bdryNodeCount))		
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(printAuxData.isDividing.begin(),
					make_transform_iterator(iStart,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(printAuxData.tmpIsActive_M.begin(),
							printAuxData.tmpNodePosX_M.begin(),
							printAuxData.tmpNodePosY_M.begin(),
							printAuxData.tmpNodeType.begin(),
							printAuxData.tmpNodeMemMirrorIndex_M.begin())), isTrue());
// std::cout<<"copy 9"<<std::endl;
// step 3 , continued  //copy cell info values ready for division /comment A&A
	thrust::counting_iterator<uint> iBegin(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.HertwigXdir.begin(),
							cellInfoVecs.HertwigYdir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.HertwigXdir.begin(),
							cellInfoVecs.HertwigYdir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(printAuxData.tmpCellRank_M.begin(),
							printAuxData.tmpDivDirX_M.begin(),
							printAuxData.tmpDivDirY_M.begin(),
                            printAuxData.tmpHertwigXdir.begin(),
                            printAuxData.tmpHertwigYdir.begin(),
							printAuxData.tmpCenterPosX_M.begin(),
							printAuxData.tmpCenterPosY_M.begin(),
							printAuxData.tmpIntAvgX_M.begin(),
							printAuxData.tmpIntAvgY_M.begin()
							)), isTrue());
// std::cout<<"copy 10"<<std::endl;
			// for (int w = 0; w < printAuxData.isDividing.size(); w++){
			// 	if (printAuxData.isDividing[w] == true){
			// 		printAuxData.tmpApicalLoc[0] = cellInfoVecs.apicalLocX[w];
			// 		printAuxData.tmpApicalLoc[1] = cellInfoVecs.apicalLocY[w];
			// 		printAuxData.tmpBasalLoc[0] = cellInfoVecs.basalLocX[w];
			// 		printAuxData.tmpBasalLoc[1] = cellInfoVecs.basalLocY[w];
			// 	}
			// }

}

// void SceCells::copyCellsEnterMitotic() {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;

// 	divAuxData.nodeStorageCount = divAuxData.toEnterMitoticCount
// 			* allocPara_m.maxAllNodePerCell;
// 			std::cout<<"nodeStorageCount = "<<divAuxData.nodeStorageCount<<std::endl;

// 	divAuxData.tmpIsActive_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, true);
// 	divAuxData.tmpNodePosX_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpNodePosY_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpNodeType = thrust::device_vector<MembraneType1>(
// 			divAuxData.nodeStorageCount, notAssigned1); //Ali 
// 	divAuxData.tmpNodeMemMirrorIndex_M = thrust::device_vector<int>(
// 			divAuxData.nodeStorageCount, -1);

// 	divAuxData.tmpCellRank_M = thrust::device_vector<uint>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpDivDirX_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpDivDirY_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpNucleusCenterPosX_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);
// 	divAuxData.tmpNucleusCenterPosY_M = thrust::device_vector<double>(
// 			divAuxData.toEnterMitoticCount, 0);

// 	divAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, false);
// 	divAuxData.tmpXPos1_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpYPos1_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);

// 	divAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, false);
// 	divAuxData.tmpXPos2_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpYPos2_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);

// // step 2 , continued // copy node info values ready for division /comment A&A
// 	thrust::counting_iterator<uint> iStart(0);
// 	// thrust::copy_if(
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(
// 	// 						nodes->getInfoVecs().nodeIsActive.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocX.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocY.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().memNodeType1.begin()
// 	// 								+ allocPara_m.bdryNodeCount)),
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(
// 	// 						nodes->getInfoVecs().nodeIsActive.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocX.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocY.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().memNodeType1.begin()
// 	// 								+ allocPara_m.bdryNodeCount))
// 	// 				+ totalNodeCountForActiveCells,
// 	// 		thrust::make_permutation_iterator(cellInfoVecs.isEnteringMitotic.begin(),
// 	// 				make_transform_iterator(iStart,
// 	// 						DivideFunctor(allocPara_m.maxAllNodePerCell))),
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
// 	// 						divAuxData.tmpNodePosX_M.begin(),
// 	// 						divAuxData.tmpNodePosY_M.begin(),
// 	// 						divAuxData.tmpNodeType.begin())), isTrue());

// 	thrust::copy_if(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().nodeIsActive.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocX.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocY.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().memNodeType1.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeMemMirrorIndex.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().nodeIsActive.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocX.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocY.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().memNodeType1.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
// 									+ allocPara_m.bdryNodeCount))
// 					+ totalNodeCountForActiveCells,
// 			thrust::make_permutation_iterator(cellInfoVecs.isEnteringMitotic.begin(),
// 					make_transform_iterator(iStart,
// 							DivideFunctor(allocPara_m.maxAllNodePerCell))),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
// 							divAuxData.tmpNodePosX_M.begin(),
// 							divAuxData.tmpNodePosY_M.begin(),
// 							divAuxData.tmpNodeType.begin(),
// 							divAuxData.tmpNodeMemMirrorIndex_M.begin())), isTrue());

// // step 3 , continued for cell properties //copy cell info values ready for division /comment A&A
// 	thrust::counting_iterator<uint> iBegin(0);
// 	thrust::copy_if(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
// 							cellInfoVecs.growthYDir.begin(),
// 							cellInfoVecs.centerCoordX.begin(),
// 							cellInfoVecs.centerCoordY.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
// 							cellInfoVecs.growthYDir.begin(),
// 							cellInfoVecs.centerCoordX.begin(),
// 							cellInfoVecs.centerCoordY.begin()))
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.isEnteringMitotic.begin(),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(divAuxData.tmpCellRank_M.begin(),
// 							divAuxData.tmpDivDirX_M.begin(),
// 							divAuxData.tmpDivDirY_M.begin(),
// 							divAuxData.tmpCenterPosX_M.begin(),
// 							divAuxData.tmpCenterPosY_M.begin())), isTrue());
// }
// void SceCells::copyCellsEnterDivision() {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;

// 	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
// 			* allocPara_m.maxAllNodePerCell;
// 			std::cout<<"nodeStorageCount = "<<divAuxData.nodeStorageCount<<std::endl;

// 	divAuxData.tmpIsActive_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, true);
// 	divAuxData.tmpNodePosX_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpNodePosY_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpNodeType = thrust::device_vector<MembraneType1>(
// 			divAuxData.nodeStorageCount, notAssigned1); //Ali 
// 	divAuxData.tmpNodeMemMirrorIndex_M = thrust::device_vector<int>(
// 			divAuxData.nodeStorageCount, -1);

// 	divAuxData.tmpCellRank_M = thrust::device_vector<uint>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpDivDirX_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpDivDirY_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpNucleusCenterPosX_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);
// 	divAuxData.tmpNucleusCenterPosY_M = thrust::device_vector<double>(
// 			divAuxData.toBeDivideCount, 0);

// 	divAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, false);
// 	divAuxData.tmpXPos1_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpYPos1_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);

// 	divAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
// 			divAuxData.nodeStorageCount, false);
// 	divAuxData.tmpXPos2_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);
// 	divAuxData.tmpYPos2_M = thrust::device_vector<double>(
// 			divAuxData.nodeStorageCount, 0.0);

// // step 2 , continued // copy node info values ready for division /comment A&A
// 	thrust::counting_iterator<uint> iStart(0);
// 	// thrust::copy_if(
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(
// 	// 						nodes->getInfoVecs().nodeIsActive.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocX.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocY.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().memNodeType1.begin()
// 	// 								+ allocPara_m.bdryNodeCount)),
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(
// 	// 						nodes->getInfoVecs().nodeIsActive.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocX.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().nodeLocY.begin()
// 	// 								+ allocPara_m.bdryNodeCount,
// 	// 						nodes->getInfoVecs().memNodeType1.begin()
// 	// 								+ allocPara_m.bdryNodeCount))
// 	// 				+ totalNodeCountForActiveCells,
// 	// 		thrust::make_permutation_iterator(cellInfoVecs.isEnteringMitotic.begin(),
// 	// 				make_transform_iterator(iStart,
// 	// 						DivideFunctor(allocPara_m.maxAllNodePerCell))),
// 	// 		thrust::make_zip_iterator(
// 	// 				thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
// 	// 						divAuxData.tmpNodePosX_M.begin(),
// 	// 						divAuxData.tmpNodePosY_M.begin(),
// 	// 						divAuxData.tmpNodeType.begin())), isTrue());

// 	thrust::copy_if(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().nodeIsActive.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocX.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocY.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().memNodeType1.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeMemMirrorIndex.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().nodeIsActive.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocX.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeLocY.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().memNodeType1.begin()
// 									+ allocPara_m.bdryNodeCount,
// 							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
// 									+ allocPara_m.bdryNodeCount))
// 					+ totalNodeCountForActiveCells,
// 			thrust::make_permutation_iterator(cellInfoVecs.isEnteringMitotic.begin(),
// 					make_transform_iterator(iStart,
// 							DivideFunctor(allocPara_m.maxAllNodePerCell))),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
// 							divAuxData.tmpNodePosX_M.begin(),
// 							divAuxData.tmpNodePosY_M.begin(),
// 							divAuxData.tmpNodeType.begin(),
// 							divAuxData.tmpNodeMemMirrorIndex_M.begin())), isTrue());

// // step 3 , continued for cell properties //copy cell info values ready for division /comment A&A
// 	thrust::counting_iterator<uint> iBegin(0);
// 	thrust::copy_if(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
// 							cellInfoVecs.growthYDir.begin(),
// 							cellInfoVecs.centerCoordX.begin(),
// 							cellInfoVecs.centerCoordY.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
// 							cellInfoVecs.growthYDir.begin(),
// 							cellInfoVecs.centerCoordX.begin(),
// 							cellInfoVecs.centerCoordY.begin()))
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.isEnteringMitotic.begin(),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(divAuxData.tmpCellRank_M.begin(),
// 							divAuxData.tmpDivDirX_M.begin(),
// 							divAuxData.tmpDivDirY_M.begin(),
// 							divAuxData.tmpCenterPosX_M.begin(),
// 							divAuxData.tmpCenterPosY_M.begin())), isTrue());
// }

void SceCells::copyCellsEnterDivision() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
			* allocPara_m.maxAllNodePerCell;
			std::cout<<"nodeStorageCount = "<<divAuxData.nodeStorageCount<<std::endl;

	divAuxData.tmpIsActive_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, true);
	divAuxData.tmpNodePosX_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodePosY_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodeType = thrust::device_vector<MembraneType1>(
			divAuxData.nodeStorageCount, notAssigned1); //Ali 
	divAuxData.tmpNodeMemMirrorIndex_M = thrust::device_vector<int>(
			divAuxData.nodeStorageCount, -1);

	divAuxData.tmpCellRank_M = thrust::device_vector<uint>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpNucleusCenterPosX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpNucleusCenterPosY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);

	divAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

	divAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

// step 2 , continued // copy node info values ready for division /comment A&A
	thrust::counting_iterator<uint> iStart(0);
	// thrust::copy_if(
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(
	// 						nodes->getInfoVecs().nodeIsActive.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocX.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocY.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().memNodeType1.begin()
	// 								+ allocPara_m.bdryNodeCount)),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(
	// 						nodes->getInfoVecs().nodeIsActive.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocX.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocY.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().memNodeType1.begin()
	// 								+ allocPara_m.bdryNodeCount))
	// 				+ totalNodeCountForActiveCells,
	// 		thrust::make_permutation_iterator(cellInfoVecs.isEnteringMitotic.begin(),
	// 				make_transform_iterator(iStart,
	// 						DivideFunctor(allocPara_m.maxAllNodePerCell))),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
	// 						divAuxData.tmpNodePosX_M.begin(),
	// 						divAuxData.tmpNodePosY_M.begin(),
	// 						divAuxData.tmpNodeType.begin())), isTrue());

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().memNodeType1.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeMemMirrorIndex.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(cellInfoVecs.isDividing.begin(),
					make_transform_iterator(iStart,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
							divAuxData.tmpNodePosX_M.begin(),
							divAuxData.tmpNodePosY_M.begin(),
							divAuxData.tmpNodeType.begin(),
							divAuxData.tmpNodeMemMirrorIndex_M.begin())), isTrue());

// step 3 , continued for cell properties //copy cell info values ready for division /comment A&A
	thrust::counting_iterator<uint> iBegin(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpCellRank_M.begin(),
							divAuxData.tmpDivDirX_M.begin(),
							divAuxData.tmpDivDirY_M.begin(),
							divAuxData.tmpCenterPosX_M.begin(),
							divAuxData.tmpCenterPosY_M.begin())), isTrue());
}




void SceCells::createTwoNewCellArr_M() {
	divAuxData.tmp1MemActiveCounts.clear();
	divAuxData.tmp1InternalActiveCounts.clear();
	divAuxData.tmp2MemActiveCounts.clear();
	divAuxData.tmp2InternalActiveCounts.clear();
	divAuxData.isMotherCellBehind.clear(); //Ali

	//divDebug();

	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		divAuxData.tmp1IntnlVec.clear();
		divAuxData.tmp2IntnlVec.clear();

		vector<CVector> membrNodes;
		vector<CVector> intnlNodes;
		vector<MembraneType1> nodeTypeIndxDiv ; 
		vector<CVector> cellCenterLine_Basal2Apical;
		vector<CVector> cellCenterLine_Apical2Basal;
		vector<CVector> cellCenterLine_Basal2Apical_leftShift;
		vector<CVector> cellCenterLine_Basal2Apical_rightShift;
		vector<CVector> cellCenterLine_Apical2Basal_leftShift;
		vector<CVector> cellCenterLine_Apical2Basal_rightShift;
		vector<double> cellCenterLine_MirrorLength_Basal2Apical;
		vector<double> cellCenterLine_MirrorLength_Apical2Basal;
		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
			// std::cout<<"HAHA ERROR 1"<<std::endl;
		// obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes,nodeTypeIndxDiv); // Ali
		obtainMembrAndIntnlNodesPlusNodeType2(i, membrNodes, intnlNodes, nodeTypeIndxDiv, 
						cellCenterLine_Basal2Apical, cellCenterLine_Apical2Basal, 
						cellCenterLine_Basal2Apical_leftShift, cellCenterLine_Basal2Apical_rightShift,
						cellCenterLine_Apical2Basal_leftShift, cellCenterLine_Apical2Basal_rightShift,
						cellCenterLine_MirrorLength_Basal2Apical, cellCenterLine_MirrorLength_Apical2Basal);  
			// std::cout<<"HAHA ERROR 2"<<std::endl;
		CVector oldCellCenter = obtainCellCenter(i);
			// std::cout<<"oldCellCenter = "<<oldCellCenter.x<<" "<<oldCellCenter.y<<std::endl;
			// std::cout<<"HAHA ERROR 3"<<std::endl;
		CVector oldNucleusCenter = obtainNucleusCenter(i, intnlNodes);// Kevin
			// std::cout<<"HAHA ERROR 4"<<std::endl;
		// CVector oldIntCenter = obtainIntCenter(i);

                //A&A commented
		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
		//		lenAlongMajorAxis);
                                              
		/*CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
				lenAlongMajorAxis);*/


		CVector divDir, true_divDir;
		// divDir.x = divAuxData.tmpHertwigXdir[i] ; //A&A
		// divDir.y = divAuxData.tmpHertwigYdir[i] ; //A&A 
		divDir.x = divAuxData.tmpBasalLoc[0] - oldNucleusCenter.x;
		divDir.y = divAuxData.tmpBasalLoc[1] - oldNucleusCenter.y;
		true_divDir = divDir.rotateNegativeNintyDeg_XY_CC();
		// std::cout<<"true_divDir = "<<true_divDir.x<<" "<<true_divDir.y<<std::endl;
			// std::cout<<"HAHA ERROR 5"<<std::endl;
		// double lenAlongHertwigAxis = calLengthAlongHertwigAxis(divDir, oldCellCenter, membrNodes);//A&A added
		double lenAlongHertwigAxis = calLengthAlongHertwigAxis(true_divDir, oldNucleusCenter, membrNodes);
			// std::cout<<"HAHA ERROR 6"<<std::endl;
 
		std::vector<VecValT> tmp1Membr, tmp2Membr;
		CVector intCell1Center, intCell2Center;
        // obtain the center of two cell along the shortest distance between the membrane nodes of mother cell. There is also a tuning factor to shift the centers inside the cell "shiftRatio"
		// obtainTwoNewIntCenters(oldIntCenter, divDir, lenAlongHertwigAxis, intCell1Center,
	  	// 		intCell2Center);
		obtainTwoNewIntCenters(oldNucleusCenter, true_divDir, lenAlongHertwigAxis, intCell1Center,
			intCell2Center);
		

		// decide each membrane nodes and internal nodes of mother cell is going to belongs to daugther cell 1 or 2. Also shrink the internal nod position along the aixs connecting mother cell to the internal nodes by a factor given as an input in the name of "Shrink ratio"
		// prepareTmpVec(i, divDir, oldCellCenter, oldIntCenter,tmp1Membr, tmp2Membr);
		prepareTmpVec(i, divDir, oldNucleusCenter, oldNucleusCenter,tmp1Membr, tmp2Membr);
		// for (int n = 0; n < tmp1Membr.size(); n++){
		// 	// std::cout<<"tmp1Membr["<<n<<"].val = "<<tmp1Membr[n].val<<", .vec = "<<tmp1Membr[n].vec.x<<" "<<tmp1Membr[n].vec.y<<" "<<tmp1Membr[n].vec.z<<", type = "<<tmp1Membr[n].type<<std::endl;
		// }
		// for (int m = 0; m < tmp2Membr.size(); m++){
		// 	// std::cout<<"tmp2Membr["<<m<<"].val = "<<tmp2Membr[m].val<<", .vec = "<<tmp2Membr[m].vec.x<<" "<<tmp2Membr[m].vec.y<<" "<<tmp2Membr[m].vec.z<<", type = "<<tmp2Membr[m].type<<std::endl;
		// }
		//create the two new membrane line based on the specified distance. 
		// processMemVec(i, tmp1Membr, tmp2Membr, oldNucleusCenter);
		processMemVec_Ver2(i, tmp1Membr, tmp2Membr, oldNucleusCenter,
						cellCenterLine_Basal2Apical, cellCenterLine_Apical2Basal, 
						cellCenterLine_Basal2Apical_leftShift, cellCenterLine_Basal2Apical_rightShift,
						cellCenterLine_Apical2Basal_leftShift, cellCenterLine_Apical2Basal_rightShift,
						cellCenterLine_MirrorLength_Basal2Apical,cellCenterLine_MirrorLength_Apical2Basal );
        // shift the internal to make sure the center of new daugther cell is exactly similar to what have chosen in the function "obtainTwoNewCenters"
		shiftIntnlNodesByCellCenter(intCell1Center, intCell2Center);
// assemble two new daughter cells information.
		assembleVecForTwoCells(i);
	}
	//divDebug();
}

void SceCells::printCellCenterLine_M(int p) {
	printAuxData.tmp1MemActiveCounts.clear();
	printAuxData.tmp1InternalActiveCounts.clear();
	printAuxData.tmp2MemActiveCounts.clear();
	printAuxData.tmp2InternalActiveCounts.clear();
	printAuxData.isMotherCellBehind.clear(); 

	for (uint i = 0; i < printAuxData.toBeDivideCount; i++) {
		printAuxData.tmp1IntnlVec.clear();
		printAuxData.tmp2IntnlVec.clear();

		vector<CVector> membrNodes;
		vector<CVector> intnlNodes;
		vector<MembraneType1> nodeTypeIndxDiv ; 
		vector<CVector> cellCenterLine_Basal2Apical;
		vector<CVector> cellCenterLine_Apical2Basal;
		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
			// std::cout<<"HAHA ERROR 1"<<std::endl;
		// obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes,nodeTypeIndxDiv); // Ali
		obtainMembrAndIntnlNodesPlusNodeType2_printingOnly(i, membrNodes, intnlNodes, nodeTypeIndxDiv, 
						cellCenterLine_Basal2Apical, cellCenterLine_Apical2Basal);  
		for (int j = 0; j < cellCenterLine_Apical2Basal.size(); j++){
			std::cout<<"cellRank = "<<p<<" "<<cellCenterLine_Apical2Basal[j].x<<" "<<cellCenterLine_Apical2Basal[j].y<<std::endl;
		}
	}
}

//A&A
// void SceCells::findHertwigAxis() {
// 	divAuxData.tmp1MemActiveCounts.clear();
// 	divAuxData.tmp1InternalActiveCounts.clear();
// 	divAuxData.tmp2MemActiveCounts.clear();
// 	divAuxData.tmp2InternalActiveCounts.clear();

// 	//divDebug();

// 	for (uint i = 0; i < divAuxData.toEnterMitoticCount; i++) {
//                 uint cellRank = divAuxData.tmpCellRank_M[i];
// 		vector<CVector> membrNodes;
// 		vector<CVector> intnlNodes;
// 		vector<MembraneType1> nodeTypeIndxDiv ; 
// 		std::pair <int ,int > ringIds ; 

// 		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
// 		obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes,nodeTypeIndxDiv); // Ali 

// 		CVector oldCellCenter = obtainCellCenter(i);// cell center
// 		double lenAlongMajorAxis;
// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis);

// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis); //Ali
// 		CVector divDir = calDivDir_ApicalBasal(oldCellCenter, membrNodes,
// 				lenAlongMajorAxis,nodeTypeIndxDiv); //Ali

//                cellInfoVecs.HertwigXdir[cellRank]=divDir.x ; 
//                cellInfoVecs.HertwigYdir[cellRank]=divDir.y ; 
//        ringIds =calApicalBasalRingIds(divDir, oldCellCenter, membrNodes,nodeTypeIndxDiv); //Ali
//        		   // it is local membrane id ;  
//                cellInfoVecs.ringApicalId[cellRank]=ringIds.first ; 
//                cellInfoVecs.ringBasalId [cellRank]=ringIds.second ;

//                std::cout<<cellInfoVecs.HertwigXdir[cellRank]<<"HertwigXdir Thrust" <<std::endl;  
//                std::cout<<cellInfoVecs.HertwigYdir[cellRank]<<"HertwigYdir Thrust" <<std::endl;  

//                std::cout<<divDir.x<<"HertwigXdir " <<std::endl;  
//                std::cout<<divDir.y<<"HertwigYdir " <<std::endl;  


// 	}
// 	//divDebug();
// }// This is the original code that find HertwigAxis based on cell center

//This is the new code that find HertwigAxis based on nucleus center
// void SceCells::findHertwigAxis() {
// 	divAuxData.tmp1MemActiveCounts.clear();
// 	divAuxData.tmp1InternalActiveCounts.clear();
// 	divAuxData.tmp2MemActiveCounts.clear();
// 	divAuxData.tmp2InternalActiveCounts.clear();

// 	//divDebug();

// 	for (uint i = 0; i < divAuxData.toEnterMitoticCount; i++) {
//                 uint cellRank = divAuxData.tmpCellRank_M[i];
// 		vector<CVector> membrNodes;
// 		vector<CVector> intnlNodes;
// 		vector<MembraneType1> nodeTypeIndxDiv ; 
// 		std::pair <int ,int > ringIds ; 

// 		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
// 		obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes, nodeTypeIndxDiv);
		

// 		// CVector oldCellCenter = obtainCellCenter(i);// cell center
// 		CVector oldCellCenter = obtainNucleusCenter(i, intnlNodes);// Kevin
// 		double lenAlongMajorAxis;
// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis);

// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis); //Ali
// 		CVector divDir = calDivDir_ApicalBasal(oldCellCenter, membrNodes,
// 				lenAlongMajorAxis,nodeTypeIndxDiv); //Ali

//                cellInfoVecs.HertwigXdir[cellRank]=divDir.x ; 
//                cellInfoVecs.HertwigYdir[cellRank]=divDir.y ; 
//        ringIds =calApicalBasalRingIds(divDir, oldCellCenter, membrNodes,nodeTypeIndxDiv); //Ali
//        		   // it is local membrane id ;  
//                cellInfoVecs.ringApicalId[cellRank]=ringIds.first ; 
//                cellInfoVecs.ringBasalId [cellRank]=ringIds.second ;

//                std::cout<<cellInfoVecs.HertwigXdir[cellRank]<<"HertwigXdir Thrust" <<std::endl;  
//                std::cout<<cellInfoVecs.HertwigYdir[cellRank]<<"HertwigYdir Thrust" <<std::endl;  

//                std::cout<<divDir.x<<"HertwigXdir " <<std::endl;  
//                std::cout<<divDir.y<<"HertwigYdir " <<std::endl;  


// 	}
// 	//divDebug();
// }

// void SceCells::findHertwigAxis_useBasalApicalLoc() {
// 	divAuxData.tmp1MemActiveCounts.clear();
// 	divAuxData.tmp1InternalActiveCounts.clear();
// 	divAuxData.tmp2MemActiveCounts.clear();
// 	divAuxData.tmp2InternalActiveCounts.clear();

// 	//divDebug();

// 	for (uint i = 0; i < divAuxData.toEnterMitoticCount; i++) {
//                 uint cellRank = divAuxData.tmpCellRank_M[i];
// 		vector<CVector> membrNodes;
// 		vector<CVector> intnlNodes;
// 		vector<MembraneType1> nodeTypeIndxDiv ; 
// 		std::pair <int ,int > ringIds ; 

// 		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
// 		obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes, nodeTypeIndxDiv);
		

// 		// CVector oldCellCenter = obtainCellCenter(i);// cell center
// 		CVector oldCellCenter = obtainNucleusCenter(i, intnlNodes);// Kevin
// 		double lenAlongMajorAxis;
// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis);

// 		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
// 		//		lenAlongMajorAxis); //Ali
// 		// CVector divDir = calDivDir_ApicalBasal(oldCellCenter, membrNodes,
// 		// 		lenAlongMajorAxis,nodeTypeIndxDiv); //Ali

// 		CVector divDir;
// 		CVector tmpDivDir;
// 		tmpDivDir.x = cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank];
// 		tmpDivDir.y = cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank];
// 		double length = sqrt((cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank])*(cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank]) +
// 						(cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank])*(cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank]));
// 		tmpDivDir.x = tmpDivDir.x/length;
// 		tmpDivDir.y = tmpDivDir.y/length;
// 		divDir = tmpDivDir.rotateNintyDeg_XY_CC();

//                cellInfoVecs.HertwigXdir[cellRank]=divDir.x ; 
//                cellInfoVecs.HertwigYdir[cellRank]=divDir.y ; 
//        ringIds =calApicalBasalRingIds(divDir, oldCellCenter, membrNodes,nodeTypeIndxDiv); //Ali
//        		   // it is local membrane id ;  
//                cellInfoVecs.ringApicalId[cellRank]=ringIds.first ; 
//                cellInfoVecs.ringBasalId [cellRank]=ringIds.second ;

//                std::cout<<cellInfoVecs.HertwigXdir[cellRank]<<"HertwigXdir Thrust" <<std::endl;  
//                std::cout<<cellInfoVecs.HertwigYdir[cellRank]<<"HertwigYdir Thrust" <<std::endl;  

//                std::cout<<divDir.x<<"HertwigXdir " <<std::endl;  
//                std::cout<<divDir.y<<"HertwigYdir " <<std::endl;  


// 	}
// 	//divDebug();
// }

void SceCells::findHertwigAxis_useBasalApicalLoc() {
	divAuxData.tmp1MemActiveCounts.clear();
	divAuxData.tmp1InternalActiveCounts.clear();
	divAuxData.tmp2MemActiveCounts.clear();
	divAuxData.tmp2InternalActiveCounts.clear();
	uint isDividingCount = 0;
	//divDebug();
	for (int i = 0; i < cellInfoVecs.isDividing.size(); i++){
		if (cellInfoVecs.isDividing[i] == true){
			isDividingCount += 1;
		}
	}
	if (isDividingCount > 1){
		std::cout<<"More than one cell undergoing division! Instability very likely to occur!"<<std::endl;
	}
	for (uint i = 0; i < isDividingCount; i++) {
                uint cellRank = divAuxData.tmpCellRank_M[i];
		vector<CVector> membrNodes;
		vector<CVector> intnlNodes;
		vector<MembraneType1> nodeTypeIndxDiv ; 
		std::pair <int ,int > ringIds ; 

		//obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);
		obtainMembrAndIntnlNodesPlusNodeType(i, membrNodes, intnlNodes, nodeTypeIndxDiv);
		

		// CVector oldCellCenter = obtainCellCenter(i);// cell center
		CVector oldCellCenter = obtainNucleusCenter(i, intnlNodes);// Kevin
		double lenAlongMajorAxis;
		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
		//		lenAlongMajorAxis);

		//CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
		//		lenAlongMajorAxis); //Ali
		// CVector divDir = calDivDir_ApicalBasal(oldCellCenter, membrNodes,
		// 		lenAlongMajorAxis,nodeTypeIndxDiv); //Ali

		CVector divDir;
		CVector tmpDivDir;
		tmpDivDir.x = cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank];
		tmpDivDir.y = cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank];
		double length = sqrt((cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank])*(cellInfoVecs.apicalLocX[cellRank] - cellInfoVecs.basalLocX[cellRank]) +
						(cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank])*(cellInfoVecs.apicalLocY[cellRank] - cellInfoVecs.basalLocY[cellRank]));
		tmpDivDir.x = tmpDivDir.x/length;
		tmpDivDir.y = tmpDivDir.y/length;
		divDir = tmpDivDir.rotateNintyDeg_XY_CC();

               cellInfoVecs.HertwigXdir[cellRank]=divDir.x ; 
               cellInfoVecs.HertwigYdir[cellRank]=divDir.y ; 
    //    ringIds =calApicalBasalRingIds(divDir, oldCellCenter, membrNodes,nodeTypeIndxDiv); //Ali
    //    		   // it is local membrane id ;  
    //            cellInfoVecs.ringApicalId[cellRank]=ringIds.first ; 
    //            cellInfoVecs.ringBasalId [cellRank]=ringIds.second ;

               std::cout<<cellInfoVecs.HertwigXdir[cellRank]<<"HertwigXdir Thrust" <<std::endl;  
               std::cout<<cellInfoVecs.HertwigYdir[cellRank]<<"HertwigYdir Thrust" <<std::endl;  

               std::cout<<divDir.x<<"HertwigXdir " <<std::endl;  
               std::cout<<divDir.y<<"HertwigYdir " <<std::endl;  


	}
	//divDebug();
}

void SceCells::copyFirstCellArr_M(double quiescence1, double quiescence1_half) {
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;

	//Ali to preserve the neighbors information of each cell for the copySecondCellArr_M  function if two neighbor cell divide at eaxctly one time step and the order
	// of mother and daughter cells are oppposite the methodology won't work. I think it almost never this situation will happen.
	// cout<<"nodeCellRankFront size = "<<nodes->getInfoVecs().nodeCellRankFront.size()<<std::endl;
	// cout<<"nodeCellRankFrontOld size = "<<nodes->getInfoVecs().nodeCellRankFrontOld.size()<<std::endl;
	// cout<<"nodeCellRankBehind size = "<<nodes->getInfoVecs().nodeCellRankBehind.size()<<std::endl;
	// cout<<"nodeCellRankBehindOld size = "<<nodes->getInfoVecs().nodeCellRankBehindOld.size()<<std::endl;
    thrust::copy (nodes->getInfoVecs().nodeCellRankFront.begin(),nodes->getInfoVecs().nodeCellRankFront.begin()+allocPara_m.currentActiveCellCount,
	              nodes->getInfoVecs().nodeCellRankFrontOld.begin()) ; 
    thrust::copy (nodes->getInfoVecs().nodeCellRankBehind.begin(),nodes->getInfoVecs().nodeCellRankBehind.begin()+allocPara_m.currentActiveCellCount,
	        	  nodes->getInfoVecs().nodeCellRankBehindOld.begin()) ; 
	cout << "Number of cells ready to divide in this time step is " <<divAuxData.toBeDivideCount << endl ; 
	if (divAuxData.toBeDivideCount>1) {
    cout << "Warnining: at Least two cells divided at the same time step chance of error in finding next neighbor of each cell"<<  endl ; 
	}
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		uint cellRank = divAuxData.tmpCellRank_M[i];
		std::cout<<"cellRank undergoing division = "<<cellRank<<std::endl;
		uint cellRankDaughter = allocPara_m.currentActiveCellCount + i; //Ali 

		//WE WILL UPDATE THE NUCLEUSLOCPERCENTAGE HERE!
		cellInfoVecs.nucleusLocPercent[cellRankDaughter] = cellInfoVecs.nucleusLocPercent[cellRank];
		std::cout<<"nucleusLocPercent Mother["<<cellRank<<"] = "<<cellInfoVecs.nucleusLocPercent[cellRank]<<std::endl;
		std::cout<<"nucleusLocPercent Daughter["<<cellRankDaughter<<"] = "<<cellInfoVecs.nucleusLocPercent[cellRankDaughter]<<std::endl;
	
		// std::cout<<"cellRankDaughter = "<<cellRankDaughter<<std::endl;
		uint nodeStartIndx = cellRank * maxAllNodePerCell
				+ allocPara_m.bdryNodeCount;
		// std::cout<<"nodeStartIndx = "<<nodeStartIndx<<std::endl;
		uint tmpStartIndx = i * maxAllNodePerCell;
		uint tmpEndIndx = (i + 1) * maxAllNodePerCell;
		// std::cout<<"tmpStartIndx = "<<tmpStartIndx<<std::endl;
		// std::cout<<"tmpEndIndx = "<<tmpEndIndx<<std::endl;
		thrust::constant_iterator<int> noAdhesion(-1), noAdhesion2(-1);

		// for (int index = 0; index < divAuxData.tmpNodeType1.size(); index++){
		// 	std::cout<<"divAuxData.tmpNodeType1["<<index<<"] = "<<divAuxData.tmpNodeType1[index]<<std::endl;
		// }

		thrust::copy(
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2,divAuxData.tmpNodeType1.begin())) + tmpStartIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2, divAuxData.tmpNodeType1.begin())) + tmpEndIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								nodes->getInfoVecs().nodeLocX.begin(),
								nodes->getInfoVecs().nodeLocY.begin(),
								nodes->getInfoVecs().nodeIsActive.begin(),
								nodes->getInfoVecs().nodeAdhereIndex.begin(),
								nodes->getInfoVecs().membrIntnlIndex.begin(),
								nodes->getInfoVecs().memNodeType1.begin()   
								)) // the 1 in memNodeType1 is not representing cell number 1 but in the rest it represents
						+ nodeStartIndx);
		cellInfoVecs.activeIntnlNodeCounts[cellRank] =
				divAuxData.tmp1InternalActiveCounts[i];
		// if (cellRank == 31 || cellRank == 86){
		// 	std::cout<<"CellRank = "<<cellRank<<"activeIntnlNodeCounts = "<<cellInfoVecs.activeIntnlNodeCounts[cellRank]<<std::endl;
		// }
		cellInfoVecs.activeMembrNodeCounts[cellRank] =
				divAuxData.tmp1MemActiveCounts[i];
		// if (cellRank == 31 || cellRank == 86){
		// 	std::cout<<"CellRank = "<<cellRank<<"activeMembrNodeCounts = "<<cellInfoVecs.activeMembrNodeCounts[cellRank]<<std::endl;
		// }
		double leftOrRight = cellInfoVecs.centerCoordX[31]*cellInfoVecs.centerCoordY[cellRank] - cellInfoVecs.centerCoordY[31]*cellInfoVecs.centerCoordX[cellRank];
		if (leftOrRight >= 0){
			cellInfoVecs.growthProgress[cellRank] = quiescence1_half; //quiescence1*0.5;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		else{
			cellInfoVecs.growthProgress[cellRank] = quiescence1;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		for (int i = cellRank*allocPara_m.maxAllNodePerCell; i < (cellRank+1)*allocPara_m.maxAllNodePerCell; i++){
			nodes->getInfoVecs().quiescencePerNode[i] = cellInfoVecs.growthProgress[cellRank];
		}
		cellInfoVecs.cellAreaGrowthProgress[cellRank] = 0;
		cellInfoVecs.cellAreaGrowthProgressNonMitotic[cellRank] = 0;
		cellInfoVecs.membrGrowProgress[cellRank] = 0.0;
		cellInfoVecs.isEnteringMitotic[cellRank] = false;
		cellInfoVecs.isRandGrowInited[cellRank] = false;
		cellInfoVecs.lastCheckPoint[cellRank] = 0;
		cellInfoVecs.distFromBasalLoc[cellRank] = 0.0;
		cellInfoVecs.distFromApicalLoc[cellRank] = 0.0;
		cellInfoVecs.daughterCellProduced[cellRank] += 1;
		// std::cout<<"divAuxData.isMotherCellBehind["<<i<<"] = "<<divAuxData.isMotherCellBehind[i]<<std::endl;
		//Ali
		// if (divAuxData.isMotherCellBehind[i]) {
		if (divAuxData.isMotherCellBehind[i]==true) {
			//nodes->getInfoVecs().nodeCellRankBehindNeighb[cellRank] =nodes->getInfoVecs().nodeCellRankBehindNeighb[cellRank] ; //as before so no need to update
	  		nodes->getInfoVecs().nodeCellRankFront[cellRank]  =cellRankDaughter ; 

			int tmpCellRankFront=nodes->getInfoVecs().nodeCellRankFrontOld[cellRank] ;  
	  		nodes->getInfoVecs().nodeCellRankBehind[tmpCellRankFront]  =cellRankDaughter ; 
			// if (cellRank == 31 || cellRank == 86){
			// 	std::cout<<"copy first cell array, cellRank = "<<cellRank<<" , "<<"tmpCellRankFront = "<<tmpCellRankFront<<std::endl;
			// 	std::cout<<"1nodeCellRankBehind["<<tmpCellRankFront<<"] = "<<nodes->getInfoVecs().nodeCellRankBehind[tmpCellRankFront]<<std::endl;
			// 	for (int j = nodeStartIndx; j < nodeStartIndx+(tmpEndIndx-tmpStartIndx); j++){
			// 		std::cout<<nodes->getInfoVecs().nodeLocX[j]<<" "<<nodes->getInfoVecs().nodeLocY[j]<<" "<<nodes->getInfoVecs().memNodeType1[j]<<std::endl;
			// 	}

			// }
		}
		else {
			nodes->getInfoVecs().nodeCellRankBehind[cellRank] =cellRankDaughter ; 
		//	nodes->getInfoVecs().nodeCellRankFrontNeighb[cellRank]  = nodes->getInfoVecs().nodeCellRankFrontNeighb[cellRank]; //as before so no need to update
			
			int tmpCellRankBehind=nodes->getInfoVecs().nodeCellRankBehindOld[cellRank] ;
	  		nodes->getInfoVecs().nodeCellRankFront[tmpCellRankBehind]  =cellRankDaughter ; 
			// if (cellRank == 31 || cellRank == 86){
			// 	std::cout<<"cellRank = "<<cellRank<<" , "<<"tmpCellRankBehind = "<<tmpCellRankBehind<<std::endl;
			// 	std::cout<<"2nodeCellRankFront["<<tmpCellRankBehind<<"] = "<<nodes->getInfoVecs().nodeCellRankFront[tmpCellRankBehind]<<std::endl;
			// 	for (int j = nodeStartIndx; j < nodeStartIndx+(tmpEndIndx-tmpStartIndx); j++){
			// 		std::cout<<nodes->getInfoVecs().nodeLocX[j]<<" "<<nodes->getInfoVecs().nodeLocY[j]<<" "<<nodes->getInfoVecs().memNodeType1[j]<<std::endl;
			// 	}
			// }
		}
		cellInfoVecs.numApicalVec[cellRank] = 0;
		cellInfoVecs.numBasalVec[cellRank] = 0;
		for (int i = cellRank*allocPara_m.maxAllNodePerCell; i < (cellRank+1)*allocPara_m.maxAllNodePerCell; i++){
			if (nodes->getInfoVecs().memNodeType1[i] == apical1){
				cellInfoVecs.numApicalVec[cellRank] += 1;
			}
			if (nodes->getInfoVecs().memNodeType1[i] == basal1){
				cellInfoVecs.numBasalVec[cellRank] += 1;
			}
		}
		cellInfoVecs.cellRankVec[cellRank] = cellRank;
		std::cout<<"Cell["<<cellRank<<"] has "<<cellInfoVecs.numApicalVec[cellRank]<<" apical nodes and "<<cellInfoVecs.numBasalVec[cellRank]<<" basal nodes initially post division"<<std::endl;
		
	// if (cellRank == 31 || cellRank == 86){
	// 	for (int k = 0; k < maxAllNodePerCell; k++){
	// 		std::cout<<"cellRank = "<<cellRank<<std::endl;
	// 		std::cout<<"nodes->getInfoVecs().memNodeType1["<<k+nodeStartIndx<<"] = "<<nodes->getInfoVecs().memNodeType1[k+nodeStartIndx]<<std::endl;
	// 	}
	// }
	}
}

void SceCells::copySecondCellArr_M(double quiescence2, double quiescence2_half) {
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {

		int cellRankMother=divAuxData.tmpCellRank_M[i] ; //Ali 
		divAuxData.cellRankMother = divAuxData.tmpCellRank_M[i];
		// std::cout<<"cellRankMother = "<<cellRankMother<<std::endl;
		uint cellRank = allocPara_m.currentActiveCellCount + i;
		divAuxData.cellRankDaughter = cellRank;
		// std::cout<<"cellRank = "<<cellRank<<std::endl;
		uint nodeStartIndx = cellRank * maxAllNodePerCell
				+ allocPara_m.bdryNodeCount;
		uint tmpStartIndx = i * maxAllNodePerCell;
		uint tmpEndIndx = (i + 1) * maxAllNodePerCell;
		// std::cout<<"tmpStartIndx for tmpNodeType2 = "<<tmpStartIndx<<std::endl;
		// std::cout<<"tmpEndIndx for tmpNodeType2 = "<<tmpEndIndx<<std::endl;
		thrust::constant_iterator<int> noAdhesion(-1), noAdhesion2(-1);
		// for (int index = 0; index < divAuxData.tmpNodeType2.size(); index++){
		// 	std::cout<<"divAuxData.tmpNodeType2["<<index<<"] = "<<divAuxData.tmpNodeType2[index]<<std::endl;
		// }
		thrust::copy(
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos2_M.begin(),
								divAuxData.tmpYPos2_M.begin(),
								divAuxData.tmpIsActive2_M.begin(), noAdhesion,
								noAdhesion2,divAuxData.tmpNodeType2.begin() )) + tmpStartIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos2_M.begin(),
								divAuxData.tmpYPos2_M.begin(),
								divAuxData.tmpIsActive2_M.begin(), noAdhesion,
								noAdhesion2,divAuxData.tmpNodeType2.begin())) + tmpEndIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								nodes->getInfoVecs().nodeLocX.begin(),
								nodes->getInfoVecs().nodeLocY.begin(),
								nodes->getInfoVecs().nodeIsActive.begin(),
								nodes->getInfoVecs().nodeAdhereIndex.begin(),
								nodes->getInfoVecs().membrIntnlIndex.begin(),
								nodes->getInfoVecs().memNodeType1.begin()))  // 1 is not representing cell 1
						+ nodeStartIndx);
		uint cellRankToCopyFrom = cellRankMother * maxAllNodePerCell+ allocPara_m.bdryNodeCount;
		double actomyoMultipMother = cellInfoVecs.contractActomyo_multip_perCell[cellRankMother];
		double actomyoMultipMother_apical = cellInfoVecs.contractActomyo_multip_apical_perCell[cellRankMother];
		cellInfoVecs.contractActomyo_multip_perCell[cellRank] = actomyoMultipMother;
		cellInfoVecs.contractActomyo_multip_apical_perCell[cellRank] = actomyoMultipMother_apical;
		for (int j = nodeStartIndx; j < nodeStartIndx+maxAllNodePerCell; j++){
			nodes->getInfoVecs().contractActomyo_multip[j] = actomyoMultipMother;//nodes->getInfoVecs().contractActomyo_multip[cellRankToCopyFrom];
			nodes->getInfoVecs().contractActomyo_multip_apical[j] = actomyoMultipMother_apical;//nodes->getInfoVecs().contractActomyo_multip_apical[cellRankToCopyFrom];
			if (j==nodeStartIndx){
				std::cout<<"contractActomyo_multip for cellRank = "<<cellRank<<" is "<<nodes->getInfoVecs().contractActomyo_multip[j]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_apical for cellRank = "<<cellRank<<" is "<<nodes->getInfoVecs().contractActomyo_multip_apical[j]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_perCell for cellRank = "<<cellRank<<" is "<<cellInfoVecs.contractActomyo_multip_perCell[cellRank]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_apical_perCell for cellRank = "<<cellRank<<" is "<<cellInfoVecs.contractActomyo_multip_apical_perCell[cellRank]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
			}
			if (j==nodeStartIndx+maxAllNodePerCell-1){
				std::cout<<"contractActomyo_multip for cellRank = "<<cellRank<<" is "<<nodes->getInfoVecs().contractActomyo_multip[j]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_apical for cellRank = "<<cellRank<<" is "<<nodes->getInfoVecs().contractActomyo_multip_apical[j]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_perCell for cellRank = "<<cellRank<<" is "<<cellInfoVecs.contractActomyo_multip_perCell[cellRank]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
				std::cout<<"contractActomyo_multip_apical_perCell for cellRank = "<<cellRank<<" is "<<cellInfoVecs.contractActomyo_multip_apical_perCell[cellRank]<<", inherited from cellRank = "<<cellRankMother<<std::endl;
			}
		}
		cellInfoVecs.activeIntnlNodeCounts[cellRank] =
				divAuxData.tmp2InternalActiveCounts[i];
		// if (cellRank == 31 || cellRank == 86){
		// 	std::cout<<"CellRank = "<<cellRank<<"activeIntnlNodeCounts = "<<cellInfoVecs.activeIntnlNodeCounts[cellRank]<<std::endl;;
		// }
		cellInfoVecs.activeMembrNodeCounts[cellRank] =
				divAuxData.tmp2MemActiveCounts[i];
		// if (cellRank == 31 || cellRank == 86){
		// 	std::cout<<"CellRank = "<<cellRank<<"activeMembrNodeCounts = "<<cellInfoVecs.activeMembrNodeCounts[cellRank]<<std::endl;;
		// }
		double leftOrRight = cellInfoVecs.centerCoordX[31]*cellInfoVecs.centerCoordY[cellRankMother] - cellInfoVecs.centerCoordY[31]*cellInfoVecs.centerCoordX[cellRankMother];
		if (leftOrRight >= 0){
			cellInfoVecs.growthProgress[cellRank] = quiescence2_half;//quiescence2*0.5;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		else{
			cellInfoVecs.growthProgress[cellRank] = quiescence2;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		// cellInfoVecs.growthProgress[cellRank] = quiescence2;
		for (int i = cellRank*allocPara_m.maxAllNodePerCell; i < (cellRank+1)*allocPara_m.maxAllNodePerCell; i++){
			nodes->getInfoVecs().quiescencePerNode[i] = cellInfoVecs.growthProgress[cellRank];
		}
		cellInfoVecs.membrGrowProgress[cellRank] = 0;
		cellInfoVecs.cellAreaGrowthProgress[cellRank] = 0;
		cellInfoVecs.cellAreaGrowthProgressNonMitotic[cellRank] = 0;
		cellInfoVecs.isEnteringMitotic[cellRank] = false;
		cellInfoVecs.isRandGrowInited[cellRank] = false;
		cellInfoVecs.lastCheckPoint[cellRank] = 0;
		cellInfoVecs.cellRoot[cellRank] = cellInfoVecs.cellRoot[cellRankMother]; //Ali 
		cellInfoVecs.eCellTypeV2[cellRank] = cellInfoVecs.eCellTypeV2[cellRankMother]; //Ali
		cellInfoVecs.distFromNucleus_normal[cellRank] = cellInfoVecs.distFromNucleus_normal[cellRankMother];
		cellInfoVecs.distFromNucleus_normal_apical[cellRank] = cellInfoVecs.distFromNucleus_normal_apical[cellRankMother];
		cellInfoVecs.distFromBasalLoc[cellRank] = 0.0;
		cellInfoVecs.distFromApicalLoc[cellRank] = 0.0;
		std::cout<<"New Cell "<<cellRank<<" introduced!"<<std::endl;
		std::cout<<"distFromNucleus_normal : "<<cellInfoVecs.distFromNucleus_normal[cellRank]<<std::endl;
		std::cout<<"distFromNucleus_normal_apical : "<<cellInfoVecs.distFromNucleus_normal_apical[cellRank]<<std::endl;
		std::cout<<"isEnteringMitotic[cellRank] : "<<cellInfoVecs.isEnteringMitotic[cellRank]<<std::endl;
//Ali
		// std::cout<<"divAuxData.isMotherCellBehind["<<i<<"] = "<<divAuxData.isMotherCellBehind[i]<<std::endl;
		if (divAuxData.isMotherCellBehind[i]==true) {
			nodes->getInfoVecs().nodeCellRankBehind[cellRank] =cellRankMother ; 
			nodes->getInfoVecs().nodeCellRankFront[cellRank]  =nodes->getInfoVecs().nodeCellRankFrontOld[cellRankMother]; 
			
		}
		else {
			nodes->getInfoVecs().nodeCellRankBehind[cellRank] =nodes->getInfoVecs().nodeCellRankBehindOld[cellRankMother]; 
			nodes->getInfoVecs().nodeCellRankFront[cellRank]  =cellRankMother ; 
			
		}
		cellInfoVecs.numApicalVec[cellRank] = 0;
		cellInfoVecs.numBasalVec[cellRank] = 0;
		for (int i = cellRank*allocPara_m.maxAllNodePerCell; i < (cellRank+1)*allocPara_m.maxAllNodePerCell; i++){
			if (nodes->getInfoVecs().memNodeType1[i] == apical1){
				cellInfoVecs.numApicalVec[cellRank] += 1;
			}
			if (nodes->getInfoVecs().memNodeType1[i] == basal1){
				cellInfoVecs.numBasalVec[cellRank] += 1;
			}
		}
		cellInfoVecs.cellRankVec[cellRank] = cellRank;
		std::cout<<"Cell["<<cellRank<<"] has "<<cellInfoVecs.numApicalVec[cellRank]<<" apical nodes and "<<cellInfoVecs.numBasalVec[cellRank]<<" basal nodes initially post division"<<std::endl;
		
	}
}

//AAMIRI
/*
void SceCells::removeCellArr_M() {
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		uint cellRank = divAuxData.tmpCellRank_M[i];
		uint nodeStartIndx = cellRank * maxAllNodePerCell
				+ allocPara_m.bdryNodeCount;
		uint tmpStartIndx = i * maxAllNodePerCell;
		uint tmpEndIndx = (i + 1) * maxAllNodePerCell;
		thrust::constant_iterator<int> noAdhesion(-1), noAdhesion2(-1);
		thrust::copy(
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2)) + tmpStartIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2)) + tmpEndIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								nodes->getInfoVecs().nodeLocX.begin(),
								nodes->getInfoVecs().nodeLocY.begin(),
								nodes->getInfoVecs().nodeIsActive.begin(),
								nodes->getInfoVecs().nodeAdhereIndex.begin(),
								nodes->getInfoVecs().membrIntnlIndex.begin()))
						+ nodeStartIndx);
		cellInfoVecs.activeIntnlNodeCounts[cellRank] =
				divAuxData.tmp1InternalActiveCounts[i];
		cellInfoVecs.activeMembrNodeCounts[cellRank] =
				divAuxData.tmp1MemActiveCounts[i];
		cellInfoVecs.growthProgress[cellRank] = 0;
		cellInfoVecs.membrGrowProgress[cellRank] = 0.0;
		cellInfoVecs.isRandGrowInited[cellRank] = false;
		cellInfoVecs.lastCheckPoint[cellRank] = 0;
	}
}

*/

void SceCells::updateActiveCellCount_M() {
	allocPara_m.currentActiveCellCount = allocPara_m.currentActiveCellCount
			+ divAuxData.toBeDivideCount;
	nodes->setActiveCellCount(allocPara_m.currentActiveCellCount);
}

//AAMIRI
/*
void SceCells::updateActiveCellCountAfterRemoval_M() {
	allocPara_m.currentActiveCellCount = allocPara_m.currentActiveCellCount
			+ divAuxData.toBeDivideCount;
	nodes->setActiveCellCount(allocPara_m.currentActiveCellCount);
}

*/

void SceCells::markIsDivideFalse_M() {
	// for (int i = 0; i < cellInfoVecs.isDividing.size(); i++){
		// if (cellInfoVecs.isDividing[i] == true){
		// 	cellInfoVecs.growthProgress[i] = 0.0;
		// 	cellInfoVecs.cellAreaGrowthProgress[i] = 0.0;
		// 	cellInfoVecs.cellAreaGrowthProgressNonMitotic[i] = 0.0;
		// 	cellInfoVecs.isEnteringMitotic[i] = false;
		// }
	// }
	thrust::fill(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin()
					+ allocPara_m.currentActiveCellCount, false);
}

void SceCells::adjustNodeVel_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ allocPara_m.bdryNodeCount + totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			ForceZero());
}

void SceCells::moveNodes_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells + allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
//Ali		SaxpyFunctorDim2(dt));
			SaxpyFunctorDim2_Damp(dt,Damp_Coef));   //Ali
}
//Ali      // This function is written to assigned different damping coefficients to cells, therefore the boundary cells can have more damping

void SceCells::moveNodes_BC_M() {
	thrust::counting_iterator<uint> iBegin2(0); 
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.Cell_Damp.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
                                                        nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.Cell_Damp.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
                                                        nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells + allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			SaxpyFunctorDim2_BC_Damp(dt)); 


}

//Ali



void SceCells::ApplyExtForces()
{ 
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

//for (int i=0 ; i <nodes->getInfoVecs().memNodeType1.size(); i++ ) { 
//	if (nodes->getInfoVecs().memNodeType1[i]==basal1) {
//		cout << "  I am a basal node with id="<< i << " and vx before applying external force is equal to " <<nodes->getInfoVecs().nodeVelX[i] << endl ;  
//	}
//}

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().memNodeType1.begin(),
                            nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().memNodeType1.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeExtForceX.begin(),
							nodes->getInfoVecs().nodeExtForceY.begin())),
			AddExtForces(curTime));
//for (int i=0 ; i <nodes->getInfoVecs().memNodeType1.size(); i++ ) { 
//	if (nodes->getInfoVecs().memNodeType1[i]==basal1) {
//		cout << "  I am a basal node with id="<< i << " and vx is equal to " <<nodes->getInfoVecs().nodeVelX[i]  << endl ; 
//	}
//}

}








void SceCells::applyMemForce_M(bool cellPolar,bool subCellPolar) {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0), iBegin1(0), iBegin2(0)  ; 
 //Ali
        thrust::fill(cellInfoVecs.Cell_Time.begin(),cellInfoVecs.Cell_Time.begin() +allocPara_m.currentActiveCellCount,curTime);
        
       //Ali 
        
         		
        thrust::device_vector<double>::iterator  MinY_Itr_Cell=thrust::min_element(
                                       cellInfoVecs.centerCoordY.begin(),
                                       cellInfoVecs.centerCoordY.begin()+allocPara_m.currentActiveCellCount ) ;

        thrust::device_vector<double>::iterator  MaxY_Itr_Cell=thrust::max_element(
                                       cellInfoVecs.centerCoordY.begin(),
                                       cellInfoVecs.centerCoordY.begin()+allocPara_m.currentActiveCellCount ) ;
        double minY_Cell= *MinY_Itr_Cell ; //This variable doesn't seemed to be used even when passed into functions //Kevin
        double maxY_Cell= *MaxY_Itr_Cell ; //This variable doesn't seemed to be used even when passed into functions //Kevin

		

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));
	int* nodeAdhereIndexAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeAdhereIndex[0])); //assuming that number of boundary nodes are equal to zero
	int* cellRootAddr = thrust::raw_pointer_cast(
			&(cellInfoVecs.cellRoot[0])); // Ali

//	if (curTime>10.05) { 
//		for (int i=0; i<nodes->getInfoVecs().nodeAdhereIndex.size(); i++) {
//			cout<<"node adhere index"<<i+allocPara_m.bdryNodeCount<<" is" <<nodes->getInfoVecs().nodeAdhereIndex[i]<<endl ; 
//		}
//		exit (EXIT_FAILURE) ; 
//	}
	//double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
	//		- growthAuxData.prolifDecay
	//				* (growthAuxData.grthProgrEndCPU
	//						- growthAuxData.grthPrgrCriVal_M_Ori);

	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 

		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.eCellTypeV2.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
                            nodes->getInfoVecs().memNodeType1.begin(),
                            nodes->getInfoVecs().isSubApicalJunction.begin(),
							make_transform_iterator(iBegin2,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin2,
									ModuloFunctor(maxAllNodePerCell)))),
																					

			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.eCellTypeV2.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
                            nodes->getInfoVecs().memNodeType1.begin(),
                            nodes->getInfoVecs().isSubApicalJunction.begin(),
							make_transform_iterator(iBegin2,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin2,
									ModuloFunctor(maxAllNodePerCell))))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeActinLevel.begin(),
			ActinLevelCal(maxAllNodePerCell,nodeIsActiveAddr,cellRootAddr,minY_Cell,maxY_Cell,cellPolar,subCellPolar));
		//double a ; 
	//for(int i=0 ;  i<totalNodeCountForActiveCells ; i++) {
	//	a=static_cast<double>(nodes->getInfoVecs().nodeAdhereIndex[i]-i);  
	//	cout<< "adhere index of node " << i << " is " << nodes->getInfoVecs().nodeAdhereIndex[i] << endl ; 
	//	cout<< "the normalized difference is" <<a/(2.0*680) <<"the difference is " << a << "2 time max node per cell is  " << 2*maxAllNodePerCell << endl ; 
//	}

	
	double* nodeActinLevelAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeActinLevel[0])); //assuming that number of boundary nodes are equal to zero
				

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        nodes->getInfoVecs().nodeAdhereIndex.begin()
									+ allocPara_m.bdryNodeCount,
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount)),

			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        nodes->getInfoVecs().nodeAdhereIndex.begin()
									+ allocPara_m.bdryNodeCount,
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().membrTensionMag.begin(),
							nodes->getInfoVecs().membrTenMagRi.begin(),
							nodes->getInfoVecs().membrLinkRiMidX.begin(),
							nodes->getInfoVecs().membrLinkRiMidY.begin(),
							nodes->getInfoVecs().membrBendLeftX.begin(),
							nodes->getInfoVecs().membrBendLeftY.begin(),
							nodes->getInfoVecs().membrBendRightX.begin(),
							nodes->getInfoVecs().membrBendRightY.begin()))
					+ allocPara_m.bdryNodeCount,
			AddMembrForce(allocPara_m.bdryNodeCount, maxAllNodePerCell,
					nodeLocXAddr, nodeLocYAddr, nodeIsActiveAddr, nodeAdhereIndexAddr,nodeActinLevelAddr, grthPrgrCriVal_M,minY_Cell,maxY_Cell));


thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
							make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                            make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().membrLinSpringEnergy.begin(),
							nodes->getInfoVecs().membrBendSpringEnergy.begin())),
			CalMembrEnergy(maxAllNodePerCell,nodeLocXAddr, nodeLocYAddr, nodeIsActiveAddr,nodeActinLevelAddr, grthPrgrCriVal_M));

energyCell.totalMembrLinSpringEnergyCell=0.5 *(thrust::reduce
 ( nodes->getInfoVecs().membrLinSpringEnergy.begin(),
   nodes->getInfoVecs().membrLinSpringEnergy.begin()+totalNodeCountForActiveCells,
  (double)0.0, thrust::plus<double>() )); 

energyCell.totalMembrBendSpringEnergyCell=thrust::reduce
 ( nodes->getInfoVecs().membrBendSpringEnergy.begin(),
   nodes->getInfoVecs().membrBendSpringEnergy.begin()+totalNodeCountForActiveCells,
  (double)0.0, thrust::plus<double>() );

energyCell.totalNodeIIEnergyCell=0.5*(thrust::reduce
 ( nodes->getInfoVecs().nodeIIEnergy.begin(),
   nodes->getInfoVecs().nodeIIEnergy.begin()+totalNodeCountForActiveCells,
  (double)0.0, thrust::plus<double>() )); 

energyCell.totalNodeIMEnergyCell=0.5*(thrust::reduce
 ( nodes->getInfoVecs().nodeIMEnergy.begin(),
   nodes->getInfoVecs().nodeIMEnergy.begin()+totalNodeCountForActiveCells,
  (double)0.0, thrust::plus<double>() )); 


energyCell.totalNodeEnergyCellOld=energyCell.totalNodeEnergyCell ;  
energyCell.totalNodeEnergyCell=energyCell.totalMembrLinSpringEnergyCell + 
						       energyCell.totalMembrBendSpringEnergyCell + 
						       energyCell.totalNodeIIEnergyCell + 
						       energyCell.totalNodeIMEnergyCell ; 






int timeStep=curTime/dt ; 
if ( (timeStep % 10000)==0 ) {

	string uniqueSymbolOutput =
			globalConfigVars.getConfigValue("UniqueSymbol").toString();
	std::string cSVFileName = "EnergyExportCell_" + uniqueSymbolOutput +  ".CSV";
	ofstream EnergyExportCell ;
	EnergyExportCell.open(cSVFileName.c_str(),ofstream::app);

	EnergyExportCell <<curTime<<","<<energyCell.totalMembrLinSpringEnergyCell << "," <<energyCell.totalMembrBendSpringEnergyCell <<
	"," <<energyCell.totalNodeIIEnergyCell<<"," <<energyCell.totalNodeIMEnergyCell<<", "<< energyCell.totalNodeEnergyCell <<std::endl;
}










	double* bendLeftXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendLeftX[0]));
	double* bendLeftYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendLeftY[0]));
	double* bendRightXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendRightX[0]));
	double* bendRightYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendRightY[0]));

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin1,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin1,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin1,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin1,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			AddMembrBend(maxAllNodePerCell, nodeIsActiveAddr, bendLeftXAddr,
					bendLeftYAddr, bendRightXAddr, bendRightYAddr));
}


//AAMIRI

void SceCells::findTangentAndNormal_M() {

	uint totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;

	thrust::counting_iterator<uint> iBegin(0), iBegin1(0);

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeF_MI_M_x.begin(), //AliE
							nodes->getInfoVecs().nodeF_MI_M_y.begin(), //AliE
 							nodes->getInfoVecs().nodeF_MI_M_T.begin(), //AliE
							nodes->getInfoVecs().nodeF_MI_M_N.begin(), //AliE
							nodes->getInfoVecs().nodeCurvature.begin(),
							nodes->getInfoVecs().nodeInterCellForceX.begin(),
							nodes->getInfoVecs().nodeInterCellForceY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin1,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin1,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeF_MI_M_x.begin(), //AliE
							nodes->getInfoVecs().nodeF_MI_M_y.begin(), //AliE
 							nodes->getInfoVecs().nodeF_MI_M_T.begin(), //AliE
							nodes->getInfoVecs().nodeF_MI_M_N.begin(), //AliE
							nodes->getInfoVecs().nodeCurvature.begin(),
							nodes->getInfoVecs().nodeInterCellForceX.begin(),
							nodes->getInfoVecs().nodeInterCellForceY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeF_MI_M_T.begin(),
							nodes->getInfoVecs().nodeF_MI_M_N.begin(),   //Absoulte value since we know it is always repulsion. only it is used for output data
							nodes->getInfoVecs().nodeCurvature.begin(),
							nodes->getInfoVecs().nodeInterCellForceTangent.begin(),
							nodes->getInfoVecs().nodeInterCellForceNormal.begin(), // Absolute value to be consittent only it is used for output data 
							nodes->getInfoVecs().membrDistToRi.begin())),
			CalCurvatures(maxAllNodePerCell, nodeIsActiveAddr, nodeLocXAddr, nodeLocYAddr));

}





void SceCells::runAblationTest(AblationEvent& ablEvent) {
	for (uint i = 0; i < ablEvent.ablationCells.size(); i++) {
		int cellRank = ablEvent.ablationCells[i].cellNum;
		std::vector<uint> removeSeq = ablEvent.ablationCells[i].nodeNums;
		cellInfoVecs.activeNodeCountOfThisCell[cellRank] =
				cellInfoVecs.activeNodeCountOfThisCell[cellRank]
						- removeSeq.size();
		nodes->removeNodes(cellRank, removeSeq);
	}
}

void SceCells::computeInternalAvgPos_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	//uint totalMembrActiveNodeCount = thrust::reduce(
	//		cellInfoVecs.activeMembrNodeCounts.begin(),
	//		cellInfoVecs.activeMembrNodeCounts.begin()
	//				+ allocPara_m.currentActiveCellCount);
	uint totalIntnlActiveNodeCount = thrust::reduce(
			cellInfoVecs.activeIntnlNodeCounts.begin(),
			cellInfoVecs.activeIntnlNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin()))
					+ allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin())),
			ActiveAndIntnl());

	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalIntnlActiveNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin())),
			thrust::equal_to<uint>(), CVec2Add());
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.activeIntnlNodeCounts.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin())), CVec2Divide());
}

void SceCells::computeInternalMaxMinPos_M() {
	// totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
	// 		* allocPara_m.maxAllNodePerCell;
	// thrust::counting_iterator<uint> iBegin(0);
	// thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	// //uint totalMembrActiveNodeCount = thrust::reduce(
	// //		cellInfoVecs.activeMembrNodeCounts.begin(),
	// //		cellInfoVecs.activeMembrNodeCounts.begin()
	// //				+ allocPara_m.currentActiveCellCount);
	// uint totalIntnlActiveNodeCount = thrust::reduce(
	// 		cellInfoVecs.activeIntnlNodeCounts.begin(),
	// 		cellInfoVecs.activeIntnlNodeCounts.begin()
	// 				+ allocPara_m.currentActiveCellCount);

	// thrust::copy_if(
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(
	// 						make_transform_iterator(iBegin,
	// 								DivideFunctor(
	// 										allocPara_m.maxAllNodePerCell)),
	// 						nodes->getInfoVecs().nodeLocX.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocY.begin()
	// 								+ allocPara_m.bdryNodeCount)),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(
	// 						make_transform_iterator(iBegin,
	// 								DivideFunctor(
	// 										allocPara_m.maxAllNodePerCell)),
	// 						nodes->getInfoVecs().nodeLocX.begin()
	// 								+ allocPara_m.bdryNodeCount,
	// 						nodes->getInfoVecs().nodeLocY.begin()
	// 								+ allocPara_m.bdryNodeCount))
	// 				+ totalNodeCountForActiveCells,
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(
	// 						nodes->getInfoVecs().nodeIsActive.begin(),
	// 						nodes->getInfoVecs().nodeCellType.begin()))
	// 				+ allocPara_m.bdryNodeCount,
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
	// 						cellNodeInfoVecs.activeXPoss.begin(),
	// 						cellNodeInfoVecs.activeYPoss.begin())),
	// 		ActiveAndIntnl());

	// thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
	// 		cellNodeInfoVecs.cellRanks.begin() + totalIntnlActiveNodeCount,
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
	// 						cellNodeInfoVecs.activeYPoss.begin())),
	// 		cellInfoVecs.cellRanksTmpStorage.begin(),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
	// 						cellInfoVecs.InternalAvgY.begin())),
	// 		thrust::equal_to<uint>(), CVec2Add());
	// thrust::transform(
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
	// 						cellInfoVecs.InternalAvgY.begin())),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
	// 						cellInfoVecs.InternalAvgY.begin()))
	// 				+ allocPara_m.currentActiveCellCount,
	// 		cellInfoVecs.activeIntnlNodeCounts.begin(),
	// 		thrust::make_zip_iterator(
	// 				thrust::make_tuple(cellInfoVecs.InternalAvgX.begin(),
	// 						cellInfoVecs.InternalAvgY.begin())), CVec2Divide());
}

void SceCells::applyVolumeConstraint(double timeRatio, double timeRatio_Crit_Division, double volume_Increase_Target_Ratio, double volume_Increase_Scale, double postDivision_restorationRateScale, int cycle) {

	calCellArea();
	// for (int k = 0; k < cellInfoVecs.cellAreaVec.size(); k++){
	// 	if (cellInfoVecs.cellAreaVec[k] < 0){
	// 		cellInfoVecs.cellAreaVec[k] = -1.0*cellInfoVecs.cellAreaVec[k];
	// 	}
	// } // !!!!!NOTE!!!!!! this is currently an ad hoc way to make sure area of each triangle is positive.

	if (timeRatio == timeRatio_Crit_Division || timeRatio == timeRatio_Crit_Division+0.2){
		std::cout<<"Current timeRatio = "<<timeRatio<<std::endl;
		for (int k = 0; k < cellInfoVecs.cellAreaVec.size(); k++){
				std::cout<<"Cell["<<k<<"] area = "<<cellInfoVecs.cellAreaVec[k]<<std::endl;
				std::cout<<"CellCenter = "<<cellInfoVecs.centerCoordX[k]<<" "<<cellInfoVecs.centerCoordY[k]<<std::endl;
		}
	}
	if (timeRatio > timeRatio_Crit_Division && nodes->isCellAreaDisplayed==false){
		for (int k = 0; k < cellInfoVecs.cellAreaVec.size(); k++){
			// if (k == 31 || k == 86){
			// 	std::cout<<"Cell["<<k<<"] area = "<<cellInfoVecs.cellAreaVec[k]<<std::endl;
			// 	std::cout<<"CellCenter = "<<cellInfoVecs.centerCoordX[k]<<" "<<cellInfoVecs.centerCoordY[k]<<std::endl;
			// }
		}
		nodes->isCellAreaDisplayed=true;
	}
	// computeLagrangeForces();

	// if (timeRatio >= timeRatio_Crit_Division){
	// 	double cellAreaDesire=(cellInfoVecs.cellAreaVec[divAuxData.cellRank_division]) + (65 - cellInfoVecs.cellAreaVec[divAuxData.cellRank_division])*postDivision_restorationRateScale*(timeRatio - timeRatio_Crit_Division )/(1.0 - timeRatio_Crit_Division);
	// 	double cellAreaDesire2=(cellInfoVecs.cellAreaVec[divAuxData.cellRank_division2]) + (65 - cellInfoVecs.cellAreaVec[divAuxData.cellRank_division2])*postDivision_restorationRateScale*(timeRatio - timeRatio_Crit_Division )/(1.0 - timeRatio_Crit_Division);
	// 	std::cout<<"cellAreaCurrent["<<divAuxData.cellRank_division<<"] = "<<cellInfoVecs.cellAreaVec[divAuxData.cellRank_division]<<" , "<<"cellAreaCurrent["<<divAuxData.cellRank_division2<<"] = "<<cellInfoVecs.cellAreaVec[divAuxData.cellRank_division2]<<std::endl;
	// 	std::cout<<"cellAreaDesire["<<divAuxData.cellRank_division<<"] = "<<cellAreaDesire<<" , "<<"cellAreaDesire["<<divAuxData.cellRank_division2<<"] = "<<cellAreaDesire2<<std::endl;
	// }
	computeLagrangeForces(timeRatio, volume_Increase_Target_Ratio);

}


void SceCells::computeCenterPos_M2() {

totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	uint totalMembrActiveNodeCount = thrust::reduce(
			cellInfoVecs.activeMembrNodeCounts.begin(),
			cellInfoVecs.activeMembrNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount);
	//uint totalIntnlActiveNodeCount = thrust::reduce(
	//		cellInfoVecs.activeIntnlNodeCounts.begin(),
	//		cellInfoVecs.activeIntnlNodeCounts.begin()
	//				+ allocPara_m.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
											nodes->getInfoVecs().nodeLocX.begin(),
											nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
											nodes->getInfoVecs().nodeLocX.begin(),
											nodes->getInfoVecs().nodeLocY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin()))
					+ allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin())),
			ActiveAndMembr());

	
	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalMembrActiveNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
							           cellNodeInfoVecs.activeYPoss.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
            						   cellInfoVecs.centerCoordY.begin())),
			thrust::equal_to<uint>(), CVec2Add());


	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.activeMembrNodeCounts.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin())), CVec2Divide());
/*
			for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {
				cout << "for cell rank "<<i<< " cell center in X direction is " << cellInfoVecs.centerCoordX[i] << endl ; 
		cout << "for cell rank "<<i<< " cell center in Y direction is " << cellInfoVecs.centerCoordY[i] << endl ; 

		}
*/


}

void SceCells::computeLagrangeForces(double timeRatio,  double volume_Increase_Target_Ratio) {


	uint maxMembrNode = allocPara_m.maxMembrNodePerCell;

	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0) ; 
	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));

	double* cellAreaVecAddr= thrust::raw_pointer_cast(
			&(cellInfoVecs.cellAreaVec[0]));
	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 

	ECellType* eCellTypeV2Addr= thrust::raw_pointer_cast(
			&(cellInfoVecs.eCellTypeV2[0]));
	bool* isEnteringMitotic = thrust::raw_pointer_cast(
			&(cellInfoVecs.isEnteringMitotic[0]));
	double* cellAreaGrowthProgress = thrust::raw_pointer_cast(
			&(cellInfoVecs.cellAreaGrowthProgress[0]));
	double* cellAreaGrowthProgressNonMitotic = thrust::raw_pointer_cast(
			&(cellInfoVecs.cellAreaGrowthProgressNonMitotic[0]));


thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                            thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
									   nodes->getInfoVecs().nodeVelY.begin(),
									   nodes->getInfoVecs().lagrangeFX.begin(),
									   nodes->getInfoVecs().lagrangeFY.begin(),
									   nodes->getInfoVecs().lagrangeFN.begin())),
			AddLagrangeForces(maxAllNodePerCell,nodeLocXAddr, nodeLocYAddr, nodeIsActiveAddr,cellAreaVecAddr,grthPrgrCriVal_M, eCellTypeV2Addr, timeRatio, 
								volume_Increase_Target_Ratio, isEnteringMitotic, cellAreaGrowthProgress, cellAreaGrowthProgressNonMitotic));
			
			uint maxNPerCell = allocPara_m.maxAllNodePerCell;
			totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
			thrust::counting_iterator<uint> iBegin2(0);


			thrust::reduce_by_key(
									make_transform_iterator(iBegin2, DivideFunctor(maxNPerCell)),
									make_transform_iterator(iBegin2, DivideFunctor(maxNPerCell))+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().lagrangeFX.begin(),
							           nodes->getInfoVecs().lagrangeFY.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.sumLagrangeFPerCellX.begin(),
							           cellInfoVecs.sumLagrangeFPerCellY.begin())),
			thrust::equal_to<uint>(), CVec2Add());
/*
			for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {
				cout << "for cell rank "<<i<< " the summation of lagrangian force in X direction is " << cellInfoVecs.sumLagrangeFPerCellX[i] << endl ; 
		cout << "for cell rank "<<i<< " the summation of lagrangian force in Y direction is " << cellInfoVecs.sumLagrangeFPerCellY[i] << endl ; 

		}
*/

}

void SceCells::computeContractileRingForces() {


	uint maxMembrNode = allocPara_m.maxMembrNodePerCell;

	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0) ; 
	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));

	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 



thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
														thrust::make_permutation_iterator(
									cellInfoVecs.ringApicalId.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                            thrust::make_permutation_iterator(
									cellInfoVecs.ringBasalId.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
														thrust::make_permutation_iterator(
									cellInfoVecs.ringApicalId.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        thrust::make_permutation_iterator(
									cellInfoVecs.ringBasalId.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
									   nodes->getInfoVecs().nodeVelY.begin())),
			AddContractileRingForces(maxAllNodePerCell,nodeLocXAddr, nodeLocYAddr, grthPrgrCriVal_M));



}


void SceCells::BC_Imp_M() {

        thrust::device_vector<double>::iterator  MinX_Itr=thrust::min_element(
                                       cellInfoVecs.centerCoordX.begin(),
                                       cellInfoVecs.centerCoordX.begin()+allocPara_m.currentActiveCellCount ) ;

        thrust::device_vector<double>::iterator  MaxX_Itr=thrust::max_element(
                                       cellInfoVecs.centerCoordX.begin(),
                                       cellInfoVecs.centerCoordX.begin()+allocPara_m.currentActiveCellCount ) ;

        thrust::device_vector<double>::iterator  MinY_Itr=thrust::min_element(
                                       cellInfoVecs.centerCoordY.begin(),
                                       cellInfoVecs.centerCoordY.begin()+allocPara_m.currentActiveCellCount ) ;

        thrust::device_vector<double>::iterator  MaxY_Itr=thrust::max_element(
                                       cellInfoVecs.centerCoordY.begin(),
                                       cellInfoVecs.centerCoordY.begin()+allocPara_m.currentActiveCellCount ) ;
        double MinX= *MinX_Itr ; 
        double MaxX= *MaxX_Itr ; 
        double MinY= *MinY_Itr ; 
        double MaxY= *MaxY_Itr ;
  
/**	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							   cellInfoVecs.centerCoordY.begin())
						           ),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							   cellInfoVecs.centerCoordY.begin())) + allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(   
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
                                                           cellInfoVecs.centerCoordY.begin())),
			BC_Tissue_Damp(Damp_Coef)) ; 


**/
        int  NumActCells=allocPara_m.currentActiveCellCount ; 

        //Ali 
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							   cellInfoVecs.centerCoordY.begin(),
                                                           cellInfoVecs.Cell_Damp.begin())
						           ),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							   cellInfoVecs.centerCoordY.begin(),
							   cellInfoVecs.Cell_Damp.begin())) + allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(   
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
                                                           cellInfoVecs.Cell_Damp.begin())),
			BC_Tissue_Damp(MinX,MaxX,MinY,MaxY,Damp_Coef,NumActCells)) ; 


/**void SceCells::randomizeGrowth() {
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)) + allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin())),
			AssignRandIfNotInit(growthAuxData.randomGrowthSpeedMin,
					growthAuxData.randomGrowthSpeedMax,
					allocPara.currentActiveCellCount,
					growthAuxData.randGenAuxPara));
}



**/


}


void SceCells::assignMemNodeType() {

	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint>  iBegin2(0)  ; 

	thrust::transform(
			thrust::make_zip_iterator(
				     thrust::make_tuple(nodes->getInfoVecs().nodeIsActive.begin(),	
										nodes->getInfoVecs().memNodeType1.begin(),
									    make_transform_iterator(iBegin2,ModuloFunctor(maxAllNodePerCell)),
									    thrust::make_permutation_iterator(
									                                     cellInfoVecs.activeMembrNodeCounts.begin(),
									                                     make_transform_iterator(iBegin2,
											                             DivideFunctor(maxAllNodePerCell))))),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeIsActive.begin(),	
									   nodes->getInfoVecs().memNodeType1.begin(),
									   make_transform_iterator(iBegin2,ModuloFunctor(maxAllNodePerCell)),
									   thrust::make_permutation_iterator(
									                                     cellInfoVecs.activeMembrNodeCounts.begin(),
									                                     make_transform_iterator(iBegin2,
											                             DivideFunctor(maxAllNodePerCell)))))								
									   + totalNodeCountForActiveCells,
										thrust::make_zip_iterator(
											thrust::make_tuple(
									   			nodes->getInfoVecs().memNodeType1.begin(),
												nodes->getInfoVecs().nodeIsApicalMem.begin(),
												nodes->getInfoVecs().nodeIsBasalMem.begin()))
												,AssignMemNodeType());

}

// This function is written with the assumption that there is at least one basal point for each cell.
void SceCells::computeBasalLoc() {

	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	int* basalNodeCountAddr = thrust::raw_pointer_cast(
			&(cellInfoVecs.basalNodeCount[0]));

	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeIsBasalMem.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.basalNodeCount.begin(),
			thrust::equal_to<uint>(), thrust::plus<int>());


	uint totalBasalNodeCount = thrust::reduce(
			cellInfoVecs.basalNodeCount.begin(),
			cellInfoVecs.basalNodeCount.begin()
					+ allocPara_m.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().memNodeType1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeLocXBasal.begin(),
							cellNodeInfoVecs.activeLocYBasal.begin())),
			ActiveAndBasal());

	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalBasalNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeLocXBasal.begin(),
									   cellNodeInfoVecs.activeLocYBasal.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.basalLocX.begin(),
					            	   cellInfoVecs.basalLocY.begin())),
			thrust::equal_to<uint>(), CVec2Add());
	// up to here basaLocX and basalLocY are the summation. We divide them //
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.basalLocX.begin(),
							           cellInfoVecs.basalLocY.begin(),
			                           cellInfoVecs.cellRanksTmpStorage.begin()
									   )),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.basalLocX.begin(),
							           cellInfoVecs.basalLocY.begin(),
			                           cellInfoVecs.cellRanksTmpStorage.begin()
									   ))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.basalLocX.begin(),
							      cellInfoVecs.basalLocY.begin())), BasalLocCal(basalNodeCountAddr));
	
}


void SceCells::computeApicalLoc(double timeRatio, double timeRatio_Crit_Division) {

	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
		// std::cout<<"AHAHA 1"<<std::endl;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	// int* apicalNodeCountAddr = thrust::raw_pointer_cast(
	// 		&(cellInfoVecs.apicalNodeCount[0]));
		// std::cout<<"AHAHA 2"<<std::endl;
		
	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeIsApicalMem.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.apicalNodeCount.begin(),
			thrust::equal_to<uint>(), thrust::plus<int>());
			
	int* apicalNodeCountAddr = thrust::raw_pointer_cast(
			&(cellInfoVecs.apicalNodeCount[0]));
// std::cout<<"AHAHA 2.5"<<std::endl;
	int sizeApical=cellInfoVecs.apicalNodeCount.size() ; 

		// std::cout<<"AHAHA 3"<<std::endl;

	uint totalApicalNodeCount = thrust::reduce(
			cellInfoVecs.apicalNodeCount.begin(),
			cellInfoVecs.apicalNodeCount.begin()
					+ allocPara_m.currentActiveCellCount);
		// std::cout<<"AHAHA 4"<<std::endl;
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().memNodeType1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeLocXApical.begin(),
							cellNodeInfoVecs.activeLocYApical.begin())),
			ActiveAndApical());


	//for (int i=sizeApical-40 ; i<sizeApical ; i++) {
	//	cout << " the location of apical node " << i << " is "<<cellNodeInfoVecs.activeLocXApical[i] << " and " << cellNodeInfoVecs.activeLocYApical[i] << endl ;   
	//}

		// std::cout<<"AHAHA 5"<<std::endl;
	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalApicalNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeLocXApical.begin(),
									   cellNodeInfoVecs.activeLocYApical.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.apicalLocX.begin(),
					            	   cellInfoVecs.apicalLocY.begin())),
			thrust::equal_to<uint>(), CVec2Add());
	// up to here apicalLocX and apicalLocY are the summation. We divide them if at lease one apical node exist.
	// 0,0 location for apical node indicates that there is no apical node.
	/* // I comment this section since for now all the cells have apical node //
	// special consideration for the cells with no apical nodes
	int  NumCellsWithApicalNode=0 ; 
	for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {
		if (cellInfoVecs.apicalNodeCount[i]!=0) {
			NumCellsWithApicalNode=NumCellsWithApicalNode +1; 
		}
	}
	*/
	//finish commenting speical consideration for the cells with no apical node 
	//simply these two are equal
	int NumCellsWithApicalNode=allocPara_m.currentActiveCellCount ; 
	//
		// std::cout<<"AHAHA 6"<<std::endl;
	//cout << "num of cells with apical node is " << NumCellsWithApicalNode << endl ; 
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
			                           cellInfoVecs.cellRanksTmpStorage.begin()
									   )),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
			                           cellInfoVecs.cellRanksTmpStorage.begin()
									   ))
					+ NumCellsWithApicalNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin())), ApicalLocCal(apicalNodeCountAddr));
	
	/* I comment this section since for this simulation all the cells have apical node
       // start special consideration for the cells which have no apical node
	   //reargment to also include the cell which have not apical cells and assign the location for them as 0,0
		for (int i=0 ; i<allocPara_m.currentActiveCellCount-1 ; i++) {  // if the cell with 0 apical node is at the end, we are fine.
			if (cellInfoVecs.apicalNodeCount[i]==0) {
				cout << " I am inside complicated loop" << endl ; 
				for (int j=allocPara_m.currentActiveCellCount-2 ; j>=i ; j--) {
					cellInfoVecs.apicalLocX[j+1]=cellInfoVecs.apicalLocX[j] ;
					cellInfoVecs.apicalLocY[j+1]=cellInfoVecs.apicalLocY[j] ;
				}
				cellInfoVecs.apicalLocX[i]=0 ;
				cellInfoVecs.apicalLocY[i]=0 ; 
			}
		}

		if (cellInfoVecs.apicalNodeCount[allocPara_m.currentActiveCellCount-1]==0) { // if the cell with 0 apical node is at the end, no rearrngment is required
			cellInfoVecs.apicalLocX[allocPara_m.currentActiveCellCount-1]=0 ;
			cellInfoVecs.apicalLocY[allocPara_m.currentActiveCellCount-1]=0 ; 
		}
	// finish special consideration for the cells that have not apical nodes 
	*/
	// if (timeRatio == timeRatio_Crit_Division){
	// 	std::cout<<"totalNodeCountForActiveCells = "<<totalNodeCountForActiveCells<<std::endl;
	// 	std::cout<<"nodes->getInfoVecs().nodeIsApicalMem = "<<nodes->getInfoVecs().nodeIsApicalMem.size()<<std::endl;
	// 	for (int i = 0; i < nodes->getInfoVecs().nodeIsApicalMem.size(); i++){
	// 		if (i == 0){
	// 			std::cout<<"i = "<<i<<", "<<nodes->getInfoVecs().nodeIsApicalMem[i]<<std::endl;
	// 		}
	// 		else if (i == nodes->getInfoVecs().nodeIsApicalMem.size()-1){
	// 			std::cout<<"i = "<<i<<", "<<nodes->getInfoVecs().nodeIsApicalMem[i]<<std::endl;
	// 		}
	// 		else{continue;}
	// 	}
	// 	std::cout<<"cellInfoVecs.cellRanksTmpStorage = "<<cellInfoVecs.cellRanksTmpStorage.size()<<std::endl;
	// 	for (int i = 0; i < cellInfoVecs.cellRanksTmpStorage.size(); i++){
	// 		if (i == 0){
	// 			std::cout<<"i = "<<i<<", "<<cellInfoVecs.cellRanksTmpStorage[i]<<std::endl;
	// 		}
	// 		else if (i == cellInfoVecs.cellRanksTmpStorage.size()-1){
	// 			std::cout<<"i = "<<i<<", "<<cellInfoVecs.cellRanksTmpStorage[i]<<std::endl;
	// 		}
	// 		else{continue;}
	// 	}
	// 	std::cout<<"cellInfoVecs.apicalNodeCount = "<<cellInfoVecs.apicalNodeCount.size()<<std::endl;
	// 	for (int i = 0; i < cellInfoVecs.apicalNodeCount.size(); i++){
	// 		if (i == 0){
	// 			std::cout<<"i = "<<i<<", "<<cellInfoVecs.apicalNodeCount[i]<<std::endl;
	// 		}
	// 		else if (i == cellInfoVecs.apicalNodeCount.size()-1){
	// 			std::cout<<"i = "<<i<<", "<<cellInfoVecs.apicalNodeCount[i]<<std::endl;
	// 		}
	// 		else{continue;}
	// 	}
	// }
}

// this function is not currently active. It is useful when the level of growth needs to be related to nucleus location.
void SceCells::computeNucleusLoc() {

		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin()))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusLocX.begin(),
							           cellInfoVecs.nucleusLocY.begin())), CalNucleusLoc());
//for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {

//	cout << "for cell rank "<< i << " Cell progress is " << cellInfoVecs.growthProgress[i] << endl ; 
//	cout << "for cell rank "<< i << " Nucleus location in X direction is " << cellInfoVecs.nucleusLocX[i] <<" in Y direction is " << cellInfoVecs.nucleusLocY[i] << endl ; 
//	cout << "for cell rank "<< i << " apical  location in X direction is " << cellInfoVecs.apicalLocX[i] <<" in Y direction is " << cellInfoVecs.apicalLocY[i] << endl ; 
//	cout << "for cell rank "<< i << " center  location in X direction is " << cellInfoVecs.centerCoordX[i] <<" in Y direction is " << cellInfoVecs.centerCoordY[i] << endl ; 
//}

}

void SceCells::computeIndividualCellHeight(double distFromNucleus_normalMax1,double distFromNucleus_normalMax2,double distFromNucleus_normalMax3,
											 double distFromNucleus_normalMax_apical1, double distFromNucleus_normalMax_apical2, double distFromNucleus_normalMax_apical3){
	// double* individualCellHeightPreMitotic = thrust::raw_pointer_cast(&(cellInfoVecs.individualCellHeight[0]));
	thrust::counting_iterator<uint> iBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
									   cellInfoVecs.basalLocX.begin(),
									   cellInfoVecs.basalLocY.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.individualCellHeight.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
									   cellInfoVecs.basalLocX.begin(),
									   cellInfoVecs.basalLocY.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.individualCellHeight.begin()))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.individualCellHeight.begin(),
							           cellInfoVecs.distFromNucleus_normal.begin(),
									   cellInfoVecs.distFromNucleus_normal_apical.begin())), CalCellHeight( distFromNucleus_normalMax1, distFromNucleus_normalMax2, distFromNucleus_normalMax3,
											 distFromNucleus_normalMax_apical1, distFromNucleus_normalMax_apical2, distFromNucleus_normalMax_apical3));

}

void SceCells::computeIndividualCellHeight_Ver2(){
	// double* individualCellHeightPreMitotic = thrust::raw_pointer_cast(&(cellInfoVecs.individualCellHeight[0]));
	thrust::counting_iterator<uint> iBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
									   cellInfoVecs.basalLocX.begin(),
									   cellInfoVecs.basalLocY.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.individualCellHeight.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin(),
									   cellInfoVecs.basalLocX.begin(),
									   cellInfoVecs.basalLocY.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.individualCellHeight.begin()))
					+ allocPara_m.currentActiveCellCount,
					cellInfoVecs.individualCellHeight.begin(), CalCellHeight_Ver2());

}


void SceCells::computeNucleusDesireLoc() {

		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusLocPercent.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusLocPercent.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin()))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusDesireLocX.begin(),
							           cellInfoVecs.nucleusDesireLocY.begin(),
									   cellInfoVecs.nucDesireDistApical.begin())), CalNucleusDesireLoc());
}

void SceCells::computeCellCenterPerturbedLoc() {

		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusLocPercent.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.nucleusLocPercent.begin(),
									   cellInfoVecs.centerCoordX.begin(),
									   cellInfoVecs.centerCoordY.begin(),
									   cellInfoVecs.apicalLocX.begin(),
							           cellInfoVecs.apicalLocY.begin()))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.cellCenterPerturbedLocLocX.begin(),
							           cellInfoVecs.cellCenterPerturbedLocLocY.begin(),
									   cellInfoVecs.cellCenterDesireDistApical.begin())), CalCellCenterPerturbedLoc());
}




void SceCells::computeNucleusIniLocPercent() {

		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgY.begin(),
									   cellInfoVecs.centerCoordY.begin(),
							           cellInfoVecs.apicalLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.InternalAvgY.begin(),
									   cellInfoVecs.centerCoordY.begin(),
							           cellInfoVecs.apicalLocY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.nucleusLocPercent.begin(), CalNucleusIniLocPercent());
/*
for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {

	cout << "for cell rank "<< i << " nucleus cell percent is " <<         cellInfoVecs.nucleusLocPercent[i] << endl ; 
	cout << "for cell rank "<< i << " cell center in Y direction is "  <<  cellInfoVecs.centerCoordY[i] << endl ; 
	cout << "for cell rank "<< i << " apical location  in Y direction is " << cellInfoVecs.apicalLocY[i] << endl ; 
	cout << "for cell rank "<< i << " Internal average in Y direction is " << cellInfoVecs.InternalAvgY[i] << endl ; 
}

*/

}

// this function is not currently active. It is used when 1) internal nodes are used to represent the nucleus 2) we wanted to force the internal nodes to be at desired location. The problem with this method is that it will create net unphysical force on the cell.
void SceCells::updateInternalAvgPosByNucleusLoc_M () {

	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);


		thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
									   cellInfoVecs.InternalAvgX.begin(),
									   cellInfoVecs.InternalAvgY.begin(),
									   cellInfoVecs.nucleusLocX.begin(),
							           cellInfoVecs.nucleusLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
									   cellInfoVecs.InternalAvgX.begin(),
									   cellInfoVecs.InternalAvgY.begin(),
									   cellInfoVecs.nucleusLocX.begin(),
							           cellInfoVecs.nucleusLocY.begin()))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.tmpShiftVecX.begin(),
							           cellInfoVecs.tmpShiftVecY.begin())), CalShiftVec());


thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.tmpShiftVecX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.tmpShiftVecY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.tmpShiftVecX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.tmpShiftVecY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
									make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(), 
							   		   nodes->getInfoVecs().nodeLocY.begin())),
			AdjustInternalNodesLoc(maxMemNodePerCell));

}




// void SceCells::growAtRandom_M(double dt) {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;
// 	cout << "dt inside growAtRandom_M is: "<< dt << endl ; 
// 	randomizeGrowth_M();

// 	updateGrowthProgress_M();

// 	decideIsScheduleToGrow_M();

// 	//computeCellTargetLength_M();

// 	//computeDistToCellCenter_M();

// 	//findMinAndMaxDistToCenter_M();

// 	//computeLenDiffExpCur_M();

// 	//stretchCellGivenLenDiff_M();

// 	addPointIfScheduledToGrow_M();

// 	//decideIsScheduleToShrink_M();// AAMIRI May5

// 	//delPointIfScheduledToGrow_M();//AAMIRI - commented out on June20

// 	int currentActiveCellCount = allocPara_m.currentActiveCellCount ; 
// 	thrust::device_vector<double>::iterator  minCellProgress_Itr=thrust::min_element(cellInfoVecs.growthProgress.begin(),
//                                               cellInfoVecs.growthProgress.begin()+ currentActiveCellCount) ;

//     double minCell_Progress= *minCellProgress_Itr ; 
//     if (minCell_Progress > 0 ) {   // to not intefer with initialization with negative progress and no cell should divide before every one is positive.
// 		adjustGrowthInfo_M(); // 
// 	}
// }
// void SceCells::growAtRandom_M(double dt) {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;
// 	// cout << "dt inside growAtRandom_M is: "<< dt << endl ; 
// 	// randomizeGrowth_M();

// 	updateGrowthProgress_M();

// 	// decideIsScheduleToGrow_M();

// 	// int currentActiveCellCount = allocPara_m.currentActiveCellCount ; 
// 	// thrust::device_vector<double>::iterator  minCellProgress_Itr=thrust::min_element(cellInfoVecs.growthProgress.begin(),
//     //                                           cellInfoVecs.growthProgress.begin()+ currentActiveCellCount) ;

//     // double minCell_Progress= *minCellProgress_Itr ; 
//     // if (minCell_Progress > 0 ) {   // to not intefer with initialization with negative progress and no cell should divide before every one is positive.
// 	// 	adjustGrowthInfo_M(); // 
// 	// }
// }
void SceCells::growAtRandom_M(double growthProgressSpeed) {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	updateGrowthProgress_M(growthProgressSpeed);
}

void SceCells::growAtRandom_M_Ver2(double growthProgressSpeed, double mitoticThreshold) {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	// cellRankFront = thrust::raw_pointer_cast(&(nodes->getInfoVecs().nodeCellRankFront[0]));
	// cellRankBehind = thrust::raw_pointer_cast(&(nodes->getInfoVecs().nodeCellRankBehind[0]));
	updateGrowthProgress_M_Ver2(growthProgressSpeed, mitoticThreshold);
}

//Ali
void SceCells::enterMitoticCheckForDivAxisCal(double mitoticThreshold) {

    bool isEnteringMitotic = decideIfAnyCellEnteringMitotic(mitoticThreshold) ; //A&A

        
        //A&A
	// if (isEnteringMitotic){
    //     std::cout<< "I am in EnteringMitotic"<< std::endl; 
	// 	copyCellsEnterMitotic();
	// 	// findHertwigAxis();
	// 	findHertwigAxis_useBasalApicalLoc();
	// }
}

void SceCells::divide2D_M(double volume_Increase_Target_Ratio, double timeRatio, double thresholdToIntroduceNewCell) {
	bool isDivisionPresent = decideIfGoingToDivide_M(volume_Increase_Target_Ratio);
	// std::cout<<"isDivisionPresent = "<<isDivisionPresent<<", number of cells undergoing division = "<<divAuxData.toBeDivideCount<<std::endl;
			// std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
			// std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
			// std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
	if (isDivisionPresent == false) {
		// std::cout<<"No cell division is present!"<<std::endl;
		return;
	}
	
	uint seed = time(NULL);
	thrust::default_random_engine rng(seed);
	rng.discard(allocPara_m.currentActiveCellCount);
		thrust::uniform_real_distribution<double> distribution(0.0, 1.0);
		thrust::uniform_real_distribution<double> distribution_half(0.0, 0.5);
		// thrust::uniform_real_distribution<double> distribution_half(0.0, 0.25);
	
	// if (timeRatio >= 0.5){
	// 	thrust::uniform_real_distribution<double> distribution(2.0, 4.0);
	// }
	double quiescence1, quiescence1_half;
	double quiescence2, quiescence2_half;
	double isNewCellIntroduced;
	// double thresholdToIntroduceNewCell = 0.3;//0.15;
	std::cout<<"Chance of new cell introduction = "<<thresholdToIntroduceNewCell<<std::endl;

	for (int i = 0; i < cellInfoVecs.isDividing.size(); i++){
		if (cellInfoVecs.isDividing[i] == true){
			if (cellInfoVecs.cellAreaVec[i] < (65*volume_Increase_Target_Ratio*0.9)){
				std::cout<<"Cell Division occurs without reaching at least 90% of the target cell volume for division! Need to readjust parameters."<<std::endl;
			}
		}
	}

	isNewCellIntroduced = distribution(rng);//(distribution(rng) - 1.0)/2.0;
	// if (1 > 0){
	if (isNewCellIntroduced < thresholdToIntroduceNewCell){
		quiescence1 = -1.0*distribution(rng);
		// quiescence1_half = -1.0*distribution_half(rng);
		quiescence1_half = 1.0*distribution_half(rng);
		quiescence2 = -1.0*distribution(rng);
		// quiescence2_half = -1.0*distribution_half(rng);
		quiescence2_half = 1.0*distribution_half(rng);
	
			// std::cout<<"cellArea[10] = "<<cellInfoVecs.cellAreaVec[10]<<std::endl;
			// std::cout<<"cellArea[19] = "<<cellInfoVecs.cellAreaVec[19]<<std::endl;
			// std::cout<<"cellArea[28] = "<<cellInfoVecs.cellAreaVec[28]<<std::endl;
		copyCellsEnterDivision();
		// std::cout<<"ERROR HERE 0?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		findHertwigAxis_useBasalApicalLoc();
		//aniDebug = true;
		// std::cout<<"ERROR HERE 1?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		copyCellsPreDivision_M(); 
		// std::cout<<"ERROR HERE 2?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		createTwoNewCellArr_M(); // main function which plays with position of internal nodes and membrane new created nodes.
		// std::cout<<"ERROR HERE 3?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		copyFirstCellArr_M(quiescence1, quiescence1_half); // copy the first cell information to GPU level and initilize values such as cell prgoress and cell rank .. 
		// std::cout<<"ERROR HERE 4?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		copySecondCellArr_M(quiescence2, quiescence2_half);// copy the second cell information to GPU level and initilize values such as cell prgoress and cell rank ..
		// std::cout<<"ERROR HERE 5?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		updateActiveCellCount_M();
		// std::cout<<"ERROR HERE 6?"<<std::endl;
		// 		std::cout<<"isDividing["<<10<<"] = "<<cellInfoVecs.isDividing[10]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[10]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[10]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[10]<<std::endl;
		// 		std::cout<<"isDividing["<<19<<"] = "<<cellInfoVecs.isDividing[19]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[19]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[19]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[19]<<std::endl;
		// 		std::cout<<"isDividing["<<28<<"] = "<<cellInfoVecs.isDividing[28]<<", isEnteringMitotic = "<<cellInfoVecs.isEnteringMitotic[28]<<", cellAreaGrowthProgress = "<<cellInfoVecs.cellAreaGrowthProgress[28]<<", cellAreaGrowthProgressNonMitotic = "<<cellInfoVecs.cellAreaGrowthProgressNonMitotic[28]<<std::endl;
		markIsDivideFalse_M();
		// std::cout<<"ERROR HERE 7?"<<std::endl;
		//divDebug();
		//Ali 
		nodes->isMemNodeTypeAssigned_postCellDivision = false;
		cellInfoVecs.isPostDivision = true;
		// for (int k = 0; k < nodes->getInfoVecs().nodeCellRankFront.size(); k++){
		uint cellRank_print = 0;
		for (int k = 0; k < allocPara_m.currentActiveCellCount; k++){
			if (k == 0){
				std::cout<<"Printing out order of cells in the tissue post growth"<<std::endl;
				std::cout<<k<<" "<<nodes->getInfoVecs().nodeCellRankFront[cellRank_print]<<" ";
			}
			// if (k == allocPara_m.currentActiveCellCount){
			// 	std::cout<<nodes->getInfoVecs().nodeCellRankFront[k]<<" "<<std::endl;;
			// }
			else{
				std::cout<<nodes->getInfoVecs().nodeCellRankFront[nodes->getInfoVecs().nodeCellRankFront[cellRank_print]]<<" ";
			}
			cellRank_print = nodes->getInfoVecs().nodeCellRankFront[cellRank_print];
			
		}
		// for (int k = 0; k < nodes->getInfoVecs().nodeCellRankFront.size(); k++){
		for (int k = 0; k < allocPara_m.currentActiveCellCount; k++){
			if (k == 0){
				std::cout<<"Printing out associated multip of cells in the tissue post growth"<<std::endl;
				std::cout<<nodes->getInfoVecs().contractActomyo_multip[k*allocPara_m.maxAllNodePerCell]<<" "<<nodes->getInfoVecs().contractActomyo_multip[nodes->getInfoVecs().nodeCellRankFront[k]*allocPara_m.maxAllNodePerCell]<<" ";
			}
			// if (k == allocPara_m.currentActiveCellCount){
			// 	std::cout<<nodes->getInfoVecs().nodeCellRankFront[k]<<" "<<std::endl;;
			// }
			else{
				// std::cout<<nodes->getInfoVecs().nodeCellRankFront[nodes->getInfoVecs().nodeCellRankFront[k-1]]<<" ";
				std::cout<< nodes->getInfoVecs().contractActomyo_multip[nodes->getInfoVecs().nodeCellRankFront[k-1]*allocPara_m.maxAllNodePerCell]<<" ";
			}
			
		}
		for (int k = 0; k < allocPara_m.currentActiveCellCount; k++){
			if (k == 0){
				std::cout<<"Printing out number of cell daughter procued"<<std::endl;
			}
			if (k == allocPara_m.currentActiveCellCount-1){
				std::cout<<"["<<k<<", "<<cellInfoVecs.daughterCellProduced[k]<<"] "<<std::endl;;
			}
			else{
				std::cout<<"["<<k<<", "<<cellInfoVecs.daughterCellProduced[k]<<"] ";
			}
			
		}
	}
	else{
		double quiescence3 = -1.0*distribution(rng);
		// double quiescence3_half = -1.0*distribution_half(rng);
		double quiescence3_half = 1.0*distribution_half(rng);
		copyCellsEnterDivision();
		uint cellRank = divAuxData.tmpCellRank_M[0];
		double leftOrRight = cellInfoVecs.centerCoordX[31]*cellInfoVecs.centerCoordY[cellRank] - cellInfoVecs.centerCoordY[31]*cellInfoVecs.centerCoordX[cellRank];
		if (leftOrRight >= 0){
			cellInfoVecs.growthProgress[cellRank] = quiescence3_half;//quiescence3*0.5;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		else{
			cellInfoVecs.growthProgress[cellRank] = quiescence3;
			std::cout<<"leftOrRight : "<<leftOrRight<<", quiescence : "<<cellInfoVecs.growthProgress[cellRank]<<std::endl;
		}
		for (int i = cellRank*allocPara_m.maxAllNodePerCell; i < (cellRank+1)*allocPara_m.maxAllNodePerCell; i++){
			nodes->getInfoVecs().quiescencePerNode[i] = cellInfoVecs.growthProgress[cellRank];//quiescence3;
		}
		// cellInfoVecs.growthProgress[cellRank] = quiescence3;
		cellInfoVecs.cellAreaGrowthProgress[cellRank] = 0;
		cellInfoVecs.cellAreaGrowthProgressNonMitotic[cellRank] = 0;
		cellInfoVecs.isEnteringMitotic[cellRank] = false;
		cellInfoVecs.isDividing[cellRank] = false;
		cellInfoVecs.distFromBasalLoc[cellRank] = 0.0;
		cellInfoVecs.distFromApicalLoc[cellRank] = 0.0;
		cellInfoVecs.daughterCellProduced[cellRank] += 1;
		for (int k = 0; k < allocPara_m.currentActiveCellCount; k++){
			if (k == 0){
				std::cout<<"Printing out number of cell daughter procued"<<std::endl;
			}
			if (k == allocPara_m.currentActiveCellCount){
				std::cout<<"["<<k<<", "<<cellInfoVecs.daughterCellProduced[k]<<" "<<std::endl;;
			}
			else{
				std::cout<<"["<<k<<", "<<cellInfoVecs.daughterCellProduced[k]<<"] ";
			}
			
		}
		int startIndex = cellRank*allocPara_m.maxAllNodePerCell + allocPara_m.maxMembrNodePerCell;
		int numOfNucleus = cellInfoVecs.activeIntnlNodeCounts[cellRank];
		std::cout<<"Current number of nucleus labeled as 'active' = "<<numOfNucleus<<", before attempting to reduce the number due to out-of-plane cell division"<<std::endl;
		// if (cellInfoVecs.activeIntnlNodeCounts[cellRank] < 48){
		// 	std::cout<<"Reducing nucleus size in out-of-plane cell division. But the current nucleus account is less than the default value (default = 24) ! SOMETHING IS WRONG!"<<std::endl;
		// }
		double activeNucNodeX = 0.0, activeNucNodeY = 0.0;
		cellInfoVecs.activeIntnlNodeCounts[cellRank] = 0;
		for (int i = 0; i < numOfNucleus; i++){
			if (i < 24){
				nodes->getInfoVecs().nodeIsActive[i + startIndex] = true;
				activeNucNodeX += nodes->getInfoVecs().nodeLocX[i + startIndex];
				activeNucNodeY += nodes->getInfoVecs().nodeLocY[i + startIndex];
				// std::cout<<"Saved nuc["<<i<<"] : "<<nodes->getInfoVecs().nodeLocX[i + startIndex]<<", "<<nodes->getInfoVecs().nodeLocY[i + startIndex]<<std::endl;
				cellInfoVecs.activeIntnlNodeCounts[cellRank] += 1;
			}
			else{
				nodes->getInfoVecs().nodeIsActive[i + startIndex] = false;
				nodes->getInfoVecs().nodeLocX[i + startIndex] = 0.0;
				nodes->getInfoVecs().nodeLocY[i + startIndex] = 0.0;
			}
		}
		activeNucNodeX = activeNucNodeX/24.0;
		activeNucNodeY = activeNucNodeY/24.0;
		// std::cout<<"activeNucNodeX, Y: "<<activeNucNodeX<<", "<<activeNucNodeY<<std::endl;
		for (int i = 0; i < 24; i++){
			if (nodes->getInfoVecs().nodeIsActive[i + startIndex] == true){
			// std::cout<<"To be moved nuc["<<i<<"] : "<<nodes->getInfoVecs().nodeLocX[i + startIndex]<<", "<<nodes->getInfoVecs().nodeLocY[i + startIndex]<<std::endl;
			// std::cout<<"shifted by X, Y ["<<i<<"]: "<<(activeNucNodeX - nodes->getInfoVecs().nodeLocX[i + startIndex])*0.9<<", "<<(activeNucNodeY - nodes->getInfoVecs().nodeLocY[i + startIndex])*0.9<<std::endl;
				double distFromNucToNucCenter = sqrt((-activeNucNodeX + nodes->getInfoVecs().nodeLocX[i + startIndex])*(-activeNucNodeX + nodes->getInfoVecs().nodeLocX[i + startIndex]) +
													(-activeNucNodeY + nodes->getInfoVecs().nodeLocY[i + startIndex])*(-activeNucNodeY + nodes->getInfoVecs().nodeLocY[i + startIndex]));
				nodes->getInfoVecs().nodeLocX[i + startIndex] = (-activeNucNodeX + nodes->getInfoVecs().nodeLocX[i + startIndex])/distFromNucToNucCenter + activeNucNodeX;//nodes->getInfoVecs().nodeLocX[i + startIndex];
				nodes->getInfoVecs().nodeLocY[i + startIndex] = (-activeNucNodeY + nodes->getInfoVecs().nodeLocY[i + startIndex])/distFromNucToNucCenter + activeNucNodeY;//nodes->getInfoVecs().nodeLocY[i + startIndex];
			// std::cout<<"Resulting nuc["<<i<<"] : "<<nodes->getInfoVecs().nodeLocX[i + startIndex]<<", "<<nodes->getInfoVecs().nodeLocY[i + startIndex]<<std::endl;
			}
		}
		// std::cout<<"Current number of nucleus labeled as 'active' = "<<numOfNucleus<<", after reducing the number for out-of-plane cell division"<<std::endl;
	}
}

void SceCells::eCMCellInteraction(bool cellPolar,bool subCellPolar, bool isInitPhase, double timeRatio, double timeRatio_Crit_ECM, double timeRatio_Crit_Division, int relaxCount, double mitoticThreshold) {
	int totalNodeCountForActiveCellsECM = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	// eCMPointerCells->ApplyECMConstrain(allocPara_m.currentActiveCellCount,totalNodeCountForActiveCellsECM,curTime,dt,Damp_Coef,cellPolar,subCellPolar,isInitPhase);//, timeRatio, timeRatio_Crit_ECM, timeRatio_Crit_Division, relaxCount);
	eCMPointerCells->ApplyECMConstrain(allocPara_m.currentActiveCellCount,totalNodeCountForActiveCellsECM,curTime,dt,Damp_Coef,cellPolar,subCellPolar,isInitPhase, mitoticThreshold);
   }




void SceCells::distributeCellGrowthProgress_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::copy(
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingEnd,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			nodes->getInfoVecs().nodeGrowPro.begin()
					+ allocPara_m.bdryNodeCount);
			if (curTime <= InitTimeStage+dt)//AAMIRI   /A & A 
				thrust::copy(
					cellInfoVecs.growthProgress.begin(),
					cellInfoVecs.growthProgress.end(),
					cellInfoVecs.lastCheckPoint.begin()
				);
}

void SceCells::allComponentsMove_M() {
	//moveNodes_M();  //Ali 
        moveNodes_BC_M();      //Ali
}

void SceCells::randomizeGrowth_M() {

thrust::device_vector<double>::iterator  MinY_Itr=thrust::min_element(nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        thrust::device_vector<double>::iterator  MaxY_Itr=thrust::max_element(nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        double minY_Tisu= *MinY_Itr ; 
        double maxY_Tisu= *MaxY_Itr ;  


	uint seed = time(NULL);
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.eCellTypeV2.begin(),
							cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.eCellTypeV2.begin(),
							cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.isRandGrowInited.begin())),
			RandomizeGrow_M(minY_Tisu,maxY_Tisu,growthAuxData.randomGrowthSpeedMin,
					growthAuxData.randomGrowthSpeedMax, seed));
	for (int i=0 ; i<1 ;  i++) {
	cout << "cell growth speed for rank " <<i << " is " << cellInfoVecs.growthSpeed [i] << endl ; 
	}
	cout << "the min growth speed is " << growthAuxData.randomGrowthSpeedMin << endl ; 
	cout << "the max growth speed is " << growthAuxData.randomGrowthSpeedMax << endl ; 


}

// void SceCells::updateGrowthProgress_M() {


// 	thrust::counting_iterator<uint> iBegin(0);
// 	thrust::counting_iterator<uint> iEnd(allocPara_m.currentActiveCellCount);

//         thrust::copy(cellInfoVecs.growthProgress.begin(),
// 			cellInfoVecs.growthProgress.begin()
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.growthProgressOld.begin());

        
// //	thrust::transform(cellInfoVecs.growthSpeed.begin(),
//   //            			cellInfoVecs.growthSpeed.begin()
// //					+ allocPara_m.currentActiveCellCount,
// //			cellInfoVecs.growthProgress.begin(),
// //			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));

// thrust::transform(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
// 							   cellInfoVecs.growthSpeed.begin(),
// 							   iBegin)),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							cellInfoVecs.growthProgress.begin()+ allocPara_m.currentActiveCellCount,
// 							cellInfoVecs.growthSpeed.begin()   + allocPara_m.currentActiveCellCount,
// 							iEnd)),
// 					cellInfoVecs.growthProgress.begin(),
// 			progress_BCImp(dt));
// }

void SceCells::updateGrowthProgress_M(double growthProgressSpeed) {


	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(allocPara_m.currentActiveCellCount);
	int* daughterCellProduced = thrust::raw_pointer_cast(
			&(cellInfoVecs.daughterCellProduced[0]));
	

        thrust::copy(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.growthProgressOld.begin());

        
//	thrust::transform(cellInfoVecs.growthSpeed.begin(),
  //            			cellInfoVecs.growthSpeed.begin()
//					+ allocPara_m.currentActiveCellCount,
//			cellInfoVecs.growthProgress.begin(),
//			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));

thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							   cellInfoVecs.growthSpeed.begin(),
							   iBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.growthProgress.begin()+ allocPara_m.currentActiveCellCount,
							cellInfoVecs.growthSpeed.begin()   + allocPara_m.currentActiveCellCount,
							iEnd)),
					cellInfoVecs.growthProgress.begin(),
			progress_BCImp(growthProgressSpeed, daughterCellProduced));
}

void SceCells::updateGrowthProgress_M_Ver2(double growthProgressSpeed, double mitoticThreshold) {


	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(allocPara_m.currentActiveCellCount);
	int* cellRankFront = thrust::raw_pointer_cast(&(nodes->getInfoVecs().nodeCellRankFront[0]));
	int* cellRankBehind = thrust::raw_pointer_cast(&(nodes->getInfoVecs().nodeCellRankBehind[0]));
	double* currentGrowthProgress = thrust::raw_pointer_cast(&cellInfoVecs.growthProgress[0]);

	int* daughterCellProduced = thrust::raw_pointer_cast(
			&(cellInfoVecs.daughterCellProduced[0]));
	

        thrust::copy(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.growthProgressOld.begin());

        
//	thrust::transform(cellInfoVecs.growthSpeed.begin(),
  //            			cellInfoVecs.growthSpeed.begin()
//					+ allocPara_m.currentActiveCellCount,
//			cellInfoVecs.growthProgress.begin(),
//			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));

thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							   cellInfoVecs.growthSpeed.begin(),
							   iBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.growthProgress.begin()+ allocPara_m.currentActiveCellCount,
							cellInfoVecs.growthSpeed.begin()   + allocPara_m.currentActiveCellCount,
							iEnd)),
					cellInfoVecs.growthProgress.begin(),
			progress_BCImp_Ver2(growthProgressSpeed, daughterCellProduced, cellRankFront, cellRankBehind, currentGrowthProgress, mitoticThreshold));
}

void SceCells::decideIsScheduleToGrow_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			PtCondiOp(miscPara.growThreshold));
}

//AAMIRI May5
void SceCells::decideIsScheduleToShrink_M() {

	double laserCenterX = 26.0;
	double laserCenterY = 25.0;
	double laserRadius = 4.0; 

	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(allocPara_m.currentActiveCellCount);


	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, 
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), 
							cellInfoVecs.isScheduledToShrink.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iEnd, 
							cellInfoVecs.centerCoordX.begin()+allocPara_m.currentActiveCellCount,
							cellInfoVecs.centerCoordY.begin()+allocPara_m.currentActiveCellCount,
							cellInfoVecs.isScheduledToShrink.begin()+allocPara_m.currentActiveCellCount)),
			cellInfoVecs.isScheduledToShrink.begin(),
			isDelOp(laserCenterX, laserCenterY, laserRadius));
}



void SceCells::computeCellTargetLength_M() {
	thrust::transform(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.expectedLength.begin(),
			CompuTarLen(bioPara.cellInitLength, bioPara.cellFinalLength));
}

void SceCells::computeDistToCellCenter_M() {
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(totalNodeCountForActiveCells);
	uint endIndx = allocPara_m.bdryNodeCount + totalNodeCountForActiveCells;
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							nodes->getInfoVecs().nodeLocX.begin() + endIndx,
							nodes->getInfoVecs().nodeLocY.begin() + endIndx,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ endIndx)),
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(), CompuDist());
}

void SceCells::findMinAndMaxDistToCenter_M() {
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.smallestDistance.begin(), thrust::equal_to<uint>(),
			thrust::minimum<double>());

// for nodes of each cell, find the maximum distance from the node to the corresponding
// cell center along the pre-defined growth direction.

	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.biggestDistance.begin(), thrust::equal_to<uint>(),
			thrust::maximum<double>());
}

void SceCells::computeLenDiffExpCur_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.lengthDifference.begin(), CompuDiff());
}

void SceCells::stretchCellGivenLenDiff_M() {
	uint count = allocPara_m.maxAllNodePerCell;
	uint bdry = allocPara_m.bdryNodeCount;
	uint actCount = totalNodeCountForActiveCells;
	uint all = bdry + actCount;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(actCount);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							nodes->getInfoVecs().nodeVelX.begin() + bdry,
							nodes->getInfoVecs().nodeVelY.begin() + bdry,
							make_transform_iterator(iBegin,
									ModuloFunctor(count)))),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin()
									+ actCount,
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							nodes->getInfoVecs().nodeVelX.begin() + all,
							nodes->getInfoVecs().nodeVelY.begin() + all,
							make_transform_iterator(iEnd,
									ModuloFunctor(count)))),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin() + bdry,
							nodes->getInfoVecs().nodeVelY.begin() + bdry)),
			ApplyStretchForce_M(bioPara.elongationCoefficient,
					allocPara_m.maxMembrNodePerCell));
}

void SceCells::addPointIfScheduledToGrow_M() {
	uint seed = time(NULL);
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(activeCellCount);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.InternalAvgX.begin(),
							cellInfoVecs.InternalAvgY.begin(), iBegin,
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.isScheduledToGrow.begin()
									+ activeCellCount,
							cellInfoVecs.activeIntnlNodeCounts.begin()
									+ activeCellCount,
							cellInfoVecs.InternalAvgX.begin() + activeCellCount,
							cellInfoVecs.InternalAvgY.begin() + activeCellCount,
							iEnd,
							cellInfoVecs.lastCheckPoint.begin()
									+ activeCellCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lastCheckPoint.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			AddPtOp_M(seed, miscPara.addNodeDistance, miscPara.growThreshold,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress,
					growthAuxData.nodeIsActiveAddress));
}


//AAMIRI
void SceCells::delPointIfScheduledToGrow_M() {
	uint seed = time(NULL);
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(activeCellCount);

	int timeStep = curTime/dt;

	if (curTime>70000.0 && curTime<70000.1){

	decideIsScheduleToShrink_M();// AAMIRI
	}

 
	if (curTime > 70000.0)
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToShrink.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), iBegin,
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.isCellActive.begin(),
							cellInfoVecs.growthSpeed.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.isScheduledToShrink.begin()
									+ activeCellCount,
							cellInfoVecs.activeIntnlNodeCounts.begin()
									+ activeCellCount,
							cellInfoVecs.centerCoordX.begin() + activeCellCount,
							cellInfoVecs.centerCoordY.begin() + activeCellCount,
							iEnd,
							cellInfoVecs.activeMembrNodeCounts.begin()
									+ activeCellCount,
							cellInfoVecs.isCellActive.begin()
									+ activeCellCount,
							cellInfoVecs.growthSpeed.begin()
									+ activeCellCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin(),												 cellInfoVecs.isCellActive.begin(),
							cellInfoVecs.growthSpeed.begin())),
			DelPtOp_M(seed, timeStep, growthAuxData.adhIndxAddr,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress,
					growthAuxData.nodeIsActiveAddress));
}



bool SceCells::decideIfGoingToDivide_M(double volume_Increase_Target_Ratio) {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.cellAreaVec.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.cellAreaVec.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			CompuIsDivide_M(65.0*volume_Increase_Target_Ratio*0.9));
	// sum all bool values which indicate whether the cell is going to divide.
	// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toBeDivideCount = thrust::reduce(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin()
					+ allocPara_m.currentActiveCellCount, (uint) (0));
	if (divAuxData.toBeDivideCount > 1){
		int howManyCellsAreDividing = 0;
		for (int i = 0; i < cellInfoVecs.isDividing.size(); i++){
			if (cellInfoVecs.isDividing[i] == true){
				howManyCellsAreDividing += 1;
				if (howManyCellsAreDividing > 1){
					cellInfoVecs.isDividing[i] = false;
				}
			}
		}
	divAuxData.toBeDivideCount = thrust::reduce(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin()
					+ allocPara_m.currentActiveCellCount, (uint) (0));
	std::cout<<"More than one cell is undergoing division. The number is trimmed down to "<<divAuxData.toBeDivideCount<<" to avoid model instability!"<<std::endl;
	}
	
	if (divAuxData.toBeDivideCount > 0) {
		return true;
	} else {
		return false;
	}
}
//A&A
// bool SceCells::decideIfAnyCellEnteringMitotic() {

//     //    double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
// 	//		- growthAuxData.prolifDecay
// 	//				* (growthAuxData.grthProgrEndCPU
// 	//						- growthAuxData.grthPrgrCriVal_M_Ori);

// 	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 
// 	thrust::transform(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
// 							cellInfoVecs.growthProgressOld.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
// 							cellInfoVecs.growthProgressOld.begin()))
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.isEnteringMitotic.begin(),
// 			CompuIsEnteringMitotic_M(grthPrgrCriVal_M));
// 			//CompuIsEnteringMitotic_M(0.98)); // Ali for cross section modeling 
// 	// sum all bool values which indicate whether the cell is going to divide.
// 	// toBeDivideCount is the total number of cells going to divide.
// 	// divAuxData.toEnterMitoticCount = thrust::reduce(cellInfoVecs.isEnteringMitotic.begin(),
// 	// 		cellInfoVecs.isEnteringMitotic.begin()
// 	// 				+ allocPara_m.currentActiveCellCount, (uint) (0));
// 	// if (cycle == 0){
// 		std::fill(cellInfoVecs.isEnteringMitotic.begin(), cellInfoVecs.isEnteringMitotic.end(),false);
// 		cellInfoVecs.isEnteringMitotic[divAuxData.cellRank_division] = true;
// 	// }
// 	// else if (cycle == 1){
// 	// 	std::fill(cellInfoVecs.isEnteringMitotic.begin(), cellInfoVecs.isEnteringMitotic.end(),false);
// 	// 	cellInfoVecs.isEnteringMitotic[divAuxData.cellRank_division] = true;
// 	// }
// 	divAuxData.toEnterMitoticCount = 1;
// 	if (divAuxData.toEnterMitoticCount > 0) {
// 		return true;
// 	} else {
// 		return false;
// 	}
// }

bool SceCells::decideIfAnyCellEnteringMitotic(double grthPrgrCriVal_M) {

    //    double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
	//		- growthAuxData.prolifDecay
	//				* (growthAuxData.grthProgrEndCPU
	//						- growthAuxData.grthPrgrCriVal_M_Ori);

	// double grthPrgrCriVal_M = 0.8973;//divAuxData.targetCellDividingArea; 
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.growthProgressOld.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.growthProgressOld.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isEnteringMitotic.begin(),
			CompuIsEnteringMitotic_M(grthPrgrCriVal_M));
			//CompuIsEnteringMitotic_M(0.98)); // Ali for cross section modeling 
	// sum all bool values which indicate whether the cell is going to divide.
	// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toEnterMitoticCount = thrust::reduce(cellInfoVecs.isEnteringMitotic.begin(),
			cellInfoVecs.isEnteringMitotic.begin()
					+ allocPara_m.currentActiveCellCount, (uint) (0));
	
	if (divAuxData.toEnterMitoticCount > 0) {
		return true;
	} else {
		return false;
	}
}


//AAMIRI
/*
bool SceCells::decideIfGoingToRemove_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isRemoving.begin(),
			CompuIsRemoving_M(allocPara_m.maxIntnlNodePerCell));
	// sum all bool values which indicate whether the cell is going to divide.
	// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toBeRemovingCount = thrust::reduce(cellInfoVecs.isRemoving.begin(),
			cellInfoVecs.isRemoving.begin()
					+ allocPara_m.currentActiveCellCount, (uint) (0));
	if (divAuxData.toBeRemovingCount > 0) {
		return true;
	} else {
		return false;
	}
}

*/

AniRawData SceCells::obtainAniRawData(AnimationCriteria& aniCri) {
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	uint beginIndx = allocPara_m.bdryNodeCount;

	AniRawData rawAniData;
	//cout << "size of potential pairs = " << pairs.size() << endl;

// unordered_map is more efficient than map, but it is a c++ 11 feature
// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

	uint maxActiveNode = activeCellCount * maxNodePerCell;
	thrust::host_vector<double> hostTmpVectorLocX(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorLocY(maxActiveNode);
	thrust::host_vector<bool> hostIsActiveVec(maxActiveNode);
	thrust::host_vector<int> hostBondVec(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorTenMag(maxActiveNode);

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin()))
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpVectorLocX.begin(),
							hostTmpVectorLocY.begin(),
							hostIsActiveVec.begin(),
							hostBondVec.begin(), hostTmpVectorTenMag.begin())));

	thrust::host_vector<uint> curActiveMemNodeCounts =
			cellInfoVecs.activeMembrNodeCounts;

	CVector tmpPos;
	uint index1;
	int index2;
	std::vector<BondInfo> bondInfoVec;

	double node1X, node1Y;
	double node2X, node2Y;
	double aniVal;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < maxMemNodePerCell; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (hostIsActiveVec[index1] == true) {
				index2 = hostBondVec[index1];
				if (index2 > index1 && index2 != -1) {
					BondInfo bond;
					bond.cellRank1 = i;
					bond.pos1 = CVector(hostTmpVectorLocX[index1],
							hostTmpVectorLocY[index1], 0);
					bond.cellRank2 = (index2 - beginIndx) / maxNodePerCell;
					bond.pos2 = CVector(hostTmpVectorLocX[index2],
							hostTmpVectorLocY[index2], 0);
					bondInfoVec.push_back(bond);
				}
			}
		}
	}

	rawAniData.bondsArr = bondInfoVec;

	uint curIndex = 0;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (j == curActiveMemNodeCounts[i] - 1) {
				index2 = beginIndx + i * maxNodePerCell;
			} else {
				index2 = beginIndx + i * maxNodePerCell + j + 1;
			}

			if (hostIsActiveVec[index1] == true
					&& hostIsActiveVec[index2] == true) {
				node1X = hostTmpVectorLocX[index1];
				node1Y = hostTmpVectorLocY[index1];
				node2X = hostTmpVectorLocX[index2];
				node2Y = hostTmpVectorLocY[index2];
				IndexMap::iterator it = locIndexToAniIndexMap.find(index1);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index1, curIndex));
					curIndex++;
					tmpPos = CVector(node1X, node1Y, 0);
					//aniVal = hostTmpVectorNodeType[index1];
					aniVal = hostTmpVectorTenMag[index1];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}
				it = locIndexToAniIndexMap.find(index2);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index2, curIndex));
					curIndex++;
					tmpPos = CVector(node2X, node2Y, 0);
					//aniVal = hostTmpVectorNodeType[index2];
					aniVal = hostTmpVectorTenMag[index2];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}

				it = locIndexToAniIndexMap.find(index1);
				uint aniIndex1 = it->second;
				it = locIndexToAniIndexMap.find(index2);
				uint aniIndex2 = it->second;

				LinkAniData linkData;
				linkData.node1Index = aniIndex1;
				linkData.node2Index = aniIndex2;
				rawAniData.memLinks.push_back(linkData);
			}
		}
	}

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < allocPara_m.maxIntnlNodePerCell; j++) {
			for (uint k = j + 1; k < allocPara_m.maxIntnlNodePerCell; k++) {
				index1 = i * maxNodePerCell + maxMemNodePerCell + j;
				index2 = i * maxNodePerCell + maxMemNodePerCell + k;
				if (hostIsActiveVec[index1] && hostIsActiveVec[index2]) {
					node1X = hostTmpVectorLocX[index1];
					node1Y = hostTmpVectorLocY[index1];
					node2X = hostTmpVectorLocX[index2];
					node2Y = hostTmpVectorLocY[index2];
					if (aniCri.isPairQualify_M(node1X, node1Y, node2X,
							node2Y)) {
						IndexMap::iterator it = locIndexToAniIndexMap.find(
								index1);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index1, curIndex));
							curIndex++;
							tmpPos = CVector(node1X, node1Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = -1;
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}
						it = locIndexToAniIndexMap.find(index2);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index2, curIndex));
							curIndex++;
							tmpPos = CVector(node2X, node2Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = -1;
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}

						it = locIndexToAniIndexMap.find(index1);
						uint aniIndex1 = it->second;
						it = locIndexToAniIndexMap.find(index2);
						uint aniIndex2 = it->second;

						LinkAniData linkData;
						linkData.node1Index = aniIndex1;
						linkData.node2Index = aniIndex2;
						rawAniData.internalLinks.push_back(linkData);
					}
				}
			}
		}
	}
	return rawAniData;
}

AniRawData SceCells::obtainAniRawDataGivenCellColor(vector<double>& cellColors,
		AnimationCriteria& aniCri, vector<double>& cellsPerimeter) {   //AliE

	uint activeCellCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	uint beginIndx = allocPara_m.bdryNodeCount;

	assert(cellColors.size() >= activeCellCount);
	assert(cellsPerimeter.size() == activeCellCount); //AliE

	AniRawData rawAniData;
	//cout << "size of potential pairs = " << pairs.size() << endl;

	// unordered_map is more efficient than map, but it is a c++ 11 feature
	// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

	uint maxActiveNode = activeCellCount * maxNodePerCell;

	thrust::host_vector<double> hostTmpVectorLocX(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorLocY(maxActiveNode);
	thrust::host_vector<bool> hostIsActiveVec(maxActiveNode);
	thrust::host_vector<int> hostBondVec(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorTenMag(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorF_MI_M_x(maxActiveNode);//AAMIRI //AliE
	thrust::host_vector<double> hostTmpVectorF_MI_M_y(maxActiveNode);//AAMIRI //AliE
	//thrust::host_vector<double> hostTmpVectorF_MI_M_T(maxActiveNode); //AliE
	//thrust::host_vector<double> hostTmpVectorF_MI_M_N(maxActiveNode);//AliE
	thrust::host_vector<double> hostTmpVectorNodeCurvature(maxActiveNode);//AAMIRI
	thrust::host_vector<double> hostTmpVectorNodeActinLevel(maxActiveNode);//Ali
	thrust::host_vector<double> hostTmpVectorInterCellForceTangent(maxActiveNode);//AAMIRI
	thrust::host_vector<double> hostTmpVectorInterCellForceNormal(maxActiveNode);//AAMIRI
	thrust::host_vector<int> hostTmpContractPair(maxActiveNode);
	// thrust::host_vector<double> hostTmpNodeVelX(maxActiveNode);
	// thrust::host_vector<double> hostTmpNodeVelY(maxActiveNode);
	thrust::host_vector<double> hostTmpNodeContrApi(maxActiveNode);



	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeF_MI_M_x.begin(),//AAMIRI //AliE
							nodes->getInfoVecs().nodeF_MI_M_y.begin(),//AAMIRI //AliE
							nodes->getInfoVecs().nodeCurvature.begin(),//AAMIRI
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin(),
							nodes->getInfoVecs().nodeInterCellForceTangent.begin(),//AAMIRI
							nodes->getInfoVecs().nodeInterCellForceNormal.begin())),//AAMIRI
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeF_MI_M_x.begin(),//AAMIRI //AliE
							nodes->getInfoVecs().nodeF_MI_M_y.begin(),//AAMIRI //AliE
							nodes->getInfoVecs().nodeCurvature.begin(),//AAMIRI
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin(),
							nodes->getInfoVecs().nodeInterCellForceTangent.begin(),//AAMIRI
							nodes->getInfoVecs().nodeInterCellForceNormal.begin()))//AAMIRI
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpVectorLocX.begin(),
							hostTmpVectorLocY.begin(), 
							hostTmpVectorF_MI_M_x.begin(), hostTmpVectorF_MI_M_y.begin(),//AAMIRI
							hostTmpVectorNodeCurvature.begin(), //AAMIRI
							hostIsActiveVec.begin(),
							hostBondVec.begin(), hostTmpVectorTenMag.begin(),
							hostTmpVectorInterCellForceTangent.begin(), hostTmpVectorInterCellForceNormal.begin())));//AAMIRI

//Copy more than 10 elements is not allowed so, I separate it
/*
	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeF_MI_M_T.begin(), //Ali
							nodes->getInfoVecs().nodeF_MI_M_N.begin(), //Ali
							nodes->getInfoVecs().nodeActinLevel.begin() //Ali
							)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeF_MI_M_T.begin(),//AliE
							nodes->getInfoVecs().nodeF_MI_M_N.begin(), //AliE
							nodes->getInfoVecs().nodeActinLevel.begin() //Ali
							))
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							hostTmpVectorF_MI_M_T.begin(), hostTmpVectorF_MI_M_N.begin(),hostTmpVectorNodeActinLevel.begin()
							)));
*/
	thrust::copy(nodes->getInfoVecs().nodeActinLevel.begin(),nodes->getInfoVecs().nodeActinLevel.begin()+ maxActiveNode,hostTmpVectorNodeActinLevel.begin()); //Ali 
	thrust::copy(nodes->getInfoVecs().basalContractPair.begin()  ,nodes->getInfoVecs().basalContractPair.begin()  + maxActiveNode,hostTmpContractPair.begin()); //Ali 
	// thrust::copy(nodes->getInfoVecs().nodeVelX.begin(), nodes->getInfoVecs().nodeVelX.begin() + maxActiveNode, hostTmpNodeVelX.begin());
	// thrust::copy(nodes->getInfoVecs().nodeVelY.begin(), nodes->getInfoVecs().nodeVelY.begin() + maxActiveNode, hostTmpNodeVelY.begin());
	thrust::copy(nodes->getInfoVecs().contractActomyo_multip_apical.begin(), nodes->getInfoVecs().contractActomyo_multip_apical.begin() + maxActiveNode, hostTmpNodeContrApi.begin());





	thrust::host_vector<uint> curActiveMemNodeCounts =
			cellInfoVecs.activeMembrNodeCounts;

	thrust::host_vector<uint> curActiveIntnlNodeCounts =
				cellInfoVecs.activeIntnlNodeCounts;


	CVector tmpPos;
	CVector tmpF_MI_M ;//AAmiri
	CVector tmpInterCellForce;//AAMIRI
	double tmpCurv;
	double tmpMembTen ; 
	double tmpActinLevel ; 
	uint index1;
	int index2;
	std::vector<BondInfo> bondInfoVec;

	double node1X, node1Y;
	double node2X, node2Y;
	double node1F_MI_M_x, node1F_MI_M_y;//AAMIRI //AliE
	double nodeInterCellForceT, nodeInterCellForceN;//AAMIRI 
	double aniVal;
	double tmpNodeVel_Mag;
	double tmpNodeContrApi;
        //double tmpF_MI_M_MagN_Int[activeCellCount-1] ; //AliE

         //This is how the VTK file is intended to be written. First the memmbraen nodes are going to be written and then internal nodes.
        //loop on membrane nodes
	for (uint i = 0; i < activeCellCount; i++) {
		//tmpF_MI_M_MagN_Int[i]=0.0   ;   
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if ( hostIsActiveVec[index1]==true) {
				tmpCurv = hostTmpVectorNodeCurvature[index1];//AAMIRI
				rawAniData.aniNodeCurvature.push_back(tmpCurv);//AAMIRI
				tmpMembTen = hostTmpVectorTenMag[index1];//Ali
				rawAniData.aniNodeMembTension.push_back(tmpMembTen);//Ali
				tmpActinLevel = hostTmpVectorNodeActinLevel[index1];//Ali
				rawAniData.aniNodeActinLevel.push_back(tmpActinLevel);//Ali

				node1F_MI_M_x= hostTmpVectorF_MI_M_x[index1]; //AliE
				node1F_MI_M_y= hostTmpVectorF_MI_M_y[index1]; //AliE
				tmpF_MI_M= CVector(node1F_MI_M_x, node1F_MI_M_y, 0.0); //AliE
				rawAniData.aniNodeF_MI_M.push_back(tmpF_MI_M); //AliE
                               // tmpF_MI_M_MagN_Int[i]=tmpF_MI_M_MagN_Int[i]+sqrt(pow(hostTmpVectorF_MI_M_x[index1],2)+pow(hostTmpVectorF_MI_M_y[index1],2)) ; //AliE
                //tmpF_MI_M_MagN_Int[i]=tmpF_MI_M_MagN_Int[i]+hostTmpVectorF_MI_M_N[index1] ; //AliE

				nodeInterCellForceT = hostTmpVectorInterCellForceTangent[index1];//AAMIRI
				nodeInterCellForceN = hostTmpVectorInterCellForceNormal[index1];//AAMIRI
				tmpInterCellForce = CVector(nodeInterCellForceT, nodeInterCellForceN, 0.0);//AAMIRI
				rawAniData.aniNodeInterCellForceArr.push_back(tmpInterCellForce);

				// tmpNodeVel_Mag = sqrt(hostTmpNodeVelX[index1]*hostTmpNodeVelX[index1] + hostTmpNodeVelY[index1]*hostTmpNodeVelY[index1]);
				// rawAniData.aniNodeVel_Mag.push_back(tmpNodeVel_Mag);
				tmpNodeContrApi = hostTmpNodeContrApi[index1];
				rawAniData.aniNodeContrApi.push_back(tmpNodeContrApi);

				rawAniData.aniNodeRank.push_back(i);//AAMIRI

				}
			
			}

	}
        //loop on internal nodes
	for (uint i=0; i<activeCellCount; i++){
			for (uint j = maxMemNodePerCell; j < maxNodePerCell; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if ( hostIsActiveVec[index1]==true ) {
				tmpCurv = hostTmpVectorNodeCurvature[index1];//AAMIRI
				rawAniData.aniNodeCurvature.push_back(tmpCurv);//AAMIRI
				tmpMembTen = hostTmpVectorTenMag[index1];//Ali
				rawAniData.aniNodeMembTension.push_back(tmpMembTen);//Ali
				tmpActinLevel = hostTmpVectorNodeActinLevel[index1];//Ali
				rawAniData.aniNodeActinLevel.push_back(tmpActinLevel);//Ali

				node1F_MI_M_x= hostTmpVectorF_MI_M_x[index1]; //AliE
				node1F_MI_M_y= hostTmpVectorF_MI_M_y[index1]; //AliE
				tmpF_MI_M= CVector(node1F_MI_M_x, node1F_MI_M_y, 0.0); //AliE
				rawAniData.aniNodeF_MI_M.push_back(tmpF_MI_M);

				nodeInterCellForceT = hostTmpVectorInterCellForceTangent[index1];//AAMIRI
				nodeInterCellForceN = hostTmpVectorInterCellForceNormal[index1];//AAMIRI
				tmpInterCellForce = CVector(nodeInterCellForceT, nodeInterCellForceN, 0.0);//AAMIRI

				// tmpNodeVel_Mag = sqrt(hostTmpNodeVelX[index1]*hostTmpNodeVelX[index1] + hostTmpNodeVelY[index1]*hostTmpNodeVelY[index1]);
				// rawAniData.aniNodeVel_Mag.push_back(tmpNodeVel_Mag);
				tmpNodeContrApi = hostTmpNodeContrApi[index1];
				rawAniData.aniNodeContrApi.push_back(tmpNodeContrApi);
				
				rawAniData.aniNodeInterCellForceArr.push_back(tmpInterCellForce);
				rawAniData.aniNodeRank.push_back(i);//AAMIRI
				}
			
			}

	}


// for adhesion pair
	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < maxMemNodePerCell; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (hostIsActiveVec[index1] == true) {
				index2 = hostBondVec[index1];
				if (index2 > index1 && index2 != -1) {
					BondInfo bond;
					bond.cellRank1 = i;
					bond.pos1 = CVector(hostTmpVectorLocX[index1],
							hostTmpVectorLocY[index1], 0);
					bond.cellRank2 = (index2 - beginIndx) / maxNodePerCell;
					bond.pos2 = CVector(hostTmpVectorLocX[index2],
							hostTmpVectorLocY[index2], 0);
					bondInfoVec.push_back(bond);
				}
			}
		}
	}

	rawAniData.bondsArr = bondInfoVec;

	uint curIndex = 0;
        //loop on membrane nodes
	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (j == curActiveMemNodeCounts[i] - 1) {
				index2 = beginIndx + i * maxNodePerCell;
			} else {
				index2 = beginIndx + i * maxNodePerCell + j + 1;
			}

			if (hostIsActiveVec[index1] == true
					&& hostIsActiveVec[index2] == true) {
				node1X = hostTmpVectorLocX[index1];
				node1Y = hostTmpVectorLocY[index1];
				node2X = hostTmpVectorLocX[index2];
				node2Y = hostTmpVectorLocY[index2];

				IndexMap::iterator it = locIndexToAniIndexMap.find(index1);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index1, curIndex));
					curIndex++;
					tmpPos = CVector(node1X, node1Y, 0);
					//aniVal = hostTmpVectorNodeType[index1];
					aniVal = cellColors[i];
                    //rawAniData.aniNodeF_MI_M_MagN_Int.push_back(tmpF_MI_M_MagN_Int[i]/cellsPerimeter[i]) ; //Ali added 
                    rawAniData.aniNodeF_MI_M_MagN_Int.push_back(cellInfoVecs.cellPressure[i]) ; //Ali added 
					
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);

				}
				it = locIndexToAniIndexMap.find(index2);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index2, curIndex));
					curIndex++;
					tmpPos = CVector(node2X, node2Y, 0);
					//aniVal = hostTmpVectorNodeType[index2];
					aniVal = cellColors[i];
                    //rawAniData.aniNodeF_MI_M_MagN_Int.push_back(tmpF_MI_M_MagN_Int[i]/cellsPerimeter[i]) ; //Ali Added
                    rawAniData.aniNodeF_MI_M_MagN_Int.push_back(cellInfoVecs.cellPressure[i]) ; //Ali Added
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);

				}

				it = locIndexToAniIndexMap.find(index1);
				uint aniIndex1 = it->second;
				it = locIndexToAniIndexMap.find(index2);
				uint aniIndex2 = it->second;

				LinkAniData linkData;
				linkData.node1Index = aniIndex1;
				linkData.node2Index = aniIndex2;
		//		rawAniData.memLinks.push_back(linkData); I don't want this type of membrane nodes links be shown.
			}
		}
	}

// loop for links for basal contraction. Since the links are between membrane nodes, no new map needs to be created.
	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			index2 = hostTmpContractPair[index1];
			if (index2 == -1) {
				continue; 
			}
			IndexMap::iterator it = locIndexToAniIndexMap.find(index1);
			uint aniIndex1 = it->second;
			it = locIndexToAniIndexMap.find(index2);
			uint aniIndex2 = it->second;

			LinkAniData linkData;
			linkData.node1Index = aniIndex1;
			linkData.node2Index = aniIndex2;
			rawAniData.memLinks.push_back(linkData);
		}
	} 


        //loop on internal nodes
	for (uint i = 0; i < activeCellCount; i++) {
	//	for (uint j = 0; j < allocPara_m.maxAllNodePerCell; j++) {
		for (uint j = 0; j < allocPara_m.maxIntnlNodePerCell; j++) {
			for (uint k = 0; k < allocPara_m.maxAllNodePerCell; k++) {   //Ali
			//for (uint k = j + 1; k < allocPara_m.maxIntnlNodePerCell; k++) {  //Ali comment 
				index1 = i * maxNodePerCell + maxMemNodePerCell + j;
				index2 = i * maxNodePerCell  + k;         //Ali
			//	index2 = i * maxNodePerCell + maxMemNodePerCell + k;   //Ali comment
			//	if (hostIsActiveVec[index1] && hostIsActiveVec[index2]) {
				if (hostIsActiveVec[index1] && hostIsActiveVec[index2]&& index1 !=index2 ) {
					node1X = hostTmpVectorLocX[index1];
					node1Y = hostTmpVectorLocY[index1];
					node2X = hostTmpVectorLocX[index2];
					node2Y = hostTmpVectorLocY[index2];
					if (aniCri.isPairQualify_M(node1X, node1Y, node2X,
							node2Y)) {
						IndexMap::iterator it = locIndexToAniIndexMap.find(
								index1);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index1, curIndex));
							curIndex++;
							tmpPos = CVector(node1X, node1Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = cellColors[i];
                            //rawAniData.aniNodeF_MI_M_MagN_Int.push_back(tmpF_MI_M_MagN_Int[i]/cellsPerimeter[i]) ; //Ali Added 
                            rawAniData.aniNodeF_MI_M_MagN_Int.push_back(cellInfoVecs.cellPressure[i]) ; //Ali Added 
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}
						it = locIndexToAniIndexMap.find(index2);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index2, curIndex));
							curIndex++;
							tmpPos = CVector(node2X, node2Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = cellColors[i];
                            //rawAniData.aniNodeF_MI_M_MagN_Int.push_back(tmpF_MI_M_MagN_Int[i]/cellsPerimeter[i]) ; //Ali Added
                            rawAniData.aniNodeF_MI_M_MagN_Int.push_back(cellInfoVecs.cellPressure[i]) ; //Ali Added
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}

						it = locIndexToAniIndexMap.find(index1);
						uint aniIndex1 = it->second;
						it = locIndexToAniIndexMap.find(index2);
						uint aniIndex2 = it->second;

						LinkAniData linkData;
						linkData.node1Index = aniIndex1;
						linkData.node2Index = aniIndex2;
					//	rawAniData.internalLinks.push_back(linkData); I don't want internal node links be shown.
					}
				}
			}
		}
	}
	return rawAniData;
}

vector<AniResumeData> SceCells::obtainResumeData() {   //AliE
	
	//Copy from GPU to CPU node properties
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	
	uint maxActiveNode = activeCellCount * maxNodePerCell;

	thrust::host_vector<double> hostTmpNodeLocX(maxActiveNode);
	thrust::host_vector<double> hostTmpNodeLocY(maxActiveNode);
	thrust::host_vector<double> hostTmpDppLevel(maxActiveNode);
	thrust::host_vector<bool>   hostTmpNodeIsActive(maxActiveNode);
	thrust::host_vector<MembraneType1>   hostTmpMemNodeType(maxActiveNode);

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(   
									   nodes->getInfoVecs().dppLevel.begin(),
									   nodes->getInfoVecs().nodeIsActive.begin(),
									   nodes->getInfoVecs().nodeLocX.begin(),
									   nodes->getInfoVecs().nodeLocY.begin(),
									   nodes->getInfoVecs().memNodeType1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
								       nodes->getInfoVecs().dppLevel.begin(),
								       nodes->getInfoVecs().nodeIsActive.begin(),
								       nodes->getInfoVecs().nodeLocX.begin(),
									   nodes->getInfoVecs().nodeLocY.begin(),
									   nodes->getInfoVecs().memNodeType1.begin()))
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpDppLevel.begin(),
									   hostTmpNodeIsActive.begin(),
									   hostTmpNodeLocX.begin(),
							           hostTmpNodeLocY.begin(),
									   hostTmpMemNodeType.begin())));  

	// Copy from GPU to CPU cell properties. Since cell vectors are small copy with thrust function seems unnecessary
	thrust::host_vector<uint> 		hostTmpActiveMemNodeCounts   =cellInfoVecs.activeMembrNodeCounts;
	thrust::host_vector<ECellType>	hostTmpCellType				 =cellInfoVecs.eCellTypeV2  ; 
	thrust::host_vector<double>hostTmpCellCntrX                  =cellInfoVecs.centerCoordX ; 
	thrust::host_vector<double>hostTmpCellCntrY 			     =cellInfoVecs.centerCoordY ;

	// Write it nicely in CPU vectorial form that can be easily wirtten in an output file.
	vector <AniResumeData> aniResumeDatas ; 
	AniResumeData membraneResumeData;
	AniResumeData internalResumeData;
	AniResumeData cellResumeData;
	
	CVector tmpPos;
	uint index1;

    //loop on membrane nodes
	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < hostTmpActiveMemNodeCounts[i]; j++) {
			index1 = i * maxNodePerCell + j;
			if ( hostTmpNodeIsActive[index1]==true) {
				membraneResumeData.cellRank.push_back(i);  // it is cell rank
				membraneResumeData.nodeType.push_back(hostTmpMemNodeType[index1]);
				membraneResumeData.signalLevel.push_back(hostTmpDppLevel[index1]); 
				
				tmpPos=CVector(hostTmpNodeLocX[index1],hostTmpNodeLocY[index1],0)  ; 
				membraneResumeData.nodePosArr.push_back(tmpPos) ;  
			}
		}
	}
	aniResumeDatas.push_back(membraneResumeData) ; 
    
	//loop on internal nodes
	for (uint i=0; i<activeCellCount; i++){
		for (uint j = maxMemNodePerCell; j < maxNodePerCell; j++) {
			index1 = i * maxNodePerCell + j;
			if ( hostTmpNodeIsActive[index1]==true ) {
				internalResumeData.cellRank.push_back(i);  // it is cell rank
				
				tmpPos=CVector(hostTmpNodeLocX[index1],hostTmpNodeLocY[index1],0)  ; 
				internalResumeData.nodePosArr.push_back(tmpPos) ;  
			}
		}
	}
	aniResumeDatas.push_back(internalResumeData) ; 
	// loop for cells 
	for (uint i=0; i<activeCellCount; i++){
		cellResumeData.cellRank.push_back(i);
		cellResumeData.cellType.push_back(hostTmpCellType[i]);
		
		tmpPos=CVector(hostTmpCellCntrX[i],hostTmpCellCntrY[i],0)  ; 
		cellResumeData.nodePosArr.push_back(tmpPos) ;  

	}
	aniResumeDatas.push_back(cellResumeData) ; 
   return aniResumeDatas;
}


void SceCells::copyInitActiveNodeCount_M(
		std::vector<uint>& initMembrActiveNodeCounts,
		std::vector<uint>& initIntnlActiveNodeCounts,
		std::vector<double> &initGrowProgVec,
		std::vector<ECellType> &eCellTypeV1) {
	assert(
			initMembrActiveNodeCounts.size()
					== initIntnlActiveNodeCounts.size());
	totalNodeCountForActiveCells = initMembrActiveNodeCounts.size()
			* allocPara_m.maxAllNodePerCell;

	thrust::copy(initMembrActiveNodeCounts.begin(),
			initMembrActiveNodeCounts.end(),
			cellInfoVecs.activeMembrNodeCounts.begin());
	thrust::copy(initIntnlActiveNodeCounts.begin(),
			initIntnlActiveNodeCounts.end(),
			cellInfoVecs.activeIntnlNodeCounts.begin());
	thrust::copy(initGrowProgVec.begin(), initGrowProgVec.end(),
			cellInfoVecs.growthProgress.begin());
	thrust::copy(eCellTypeV1.begin(), eCellTypeV1.end(),
			cellInfoVecs.eCellTypeV2.begin());   // v2 might be bigger
	//for (int i=0 ; i<eCellTypeV1.size() ; i++ ) {
	//	cout << "fourth check for cell type" << cellInfoVecs.eCellTypeV2[i] << endl ; 
//	}
}

void SceCells::myDebugFunction() {

	uint maxActiveNodeCount = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxActiveCellCount = allocPara_m.currentActiveCellCount;
	std::cout << "totalNodeCountforActiveCells: "
			<< totalNodeCountForActiveCells << std::endl;
	std::cout << "maxAllNodePerCell: " << allocPara_m.maxAllNodePerCell
			<< std::endl;
	std::cout << "maxActiveCellCount: " << maxActiveCellCount << std::endl;
	std::cout << "bdryNodeCount: " << allocPara_m.bdryNodeCount << std::endl;

	std::cout << "grow threshold: " << miscPara.growThreshold << std::endl;

	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthProgress[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.isScheduledToGrow[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.lastCheckPoint[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeCount; i++) {
		if (nodes->getInfoVecs().nodeIsActive[i]
				&& nodes->getInfoVecs().nodeCellType[i] == CellIntnl) {
			std::cout << nodes->getInfoVecs().nodeVelX[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.activeIntnlNodeCounts[i] << " ";
	}
	std::cout << std::endl;

	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.expectedLength[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.smallestDistance[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.biggestDistance[i] << " ";
	}
	std::cout << std::endl;

	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.lengthDifference[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.centerCoordX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.centerCoordY[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthXDir[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthYDir[i] << " ";
	}
	std::cout << std::endl;

	int jj;
	std::cin >> jj;
}

void SceCells::divDebug() {

	std::cout << "tmpIsActive_M: ";
	for (uint i = 0; i < divAuxData.tmpIsActive_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpNodePosX_M: ";
	for (uint i = 0; i < divAuxData.tmpNodePosX_M.size(); i++) {
		std::cout << divAuxData.tmpNodePosX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpNodePosY_M : ";
	for (uint i = 0; i < divAuxData.tmpNodePosY_M.size(); i++) {
		std::cout << divAuxData.tmpNodePosY_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpCellRank_M : ";
	for (uint i = 0; i < divAuxData.tmpCellRank_M.size(); i++) {
		std::cout << divAuxData.tmpCellRank_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpDivDirX_M : ";
	for (uint i = 0; i < divAuxData.tmpDivDirX_M.size(); i++) {
		std::cout << divAuxData.tmpDivDirX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpDivDirY_M : ";
	for (uint i = 0; i < divAuxData.tmpDivDirY_M.size(); i++) {
		std::cout << divAuxData.tmpDivDirY_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpCenterPosX_M : ";
	for (uint i = 0; i < divAuxData.tmpCenterPosX_M.size(); i++) {
		std::cout << divAuxData.tmpCenterPosX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpCenterPosY_M : ";
	for (uint i = 0; i < divAuxData.tmpCenterPosY_M.size(); i++) {
		std::cout << divAuxData.tmpCenterPosY_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpIsActive1_M : ";
	for (uint i = 0; i < divAuxData.tmpIsActive1_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive1_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpXPos1_M : ";
	for (uint i = 0; i < divAuxData.tmpXPos1_M.size(); i++) {
		std::cout << divAuxData.tmpXPos1_M[i] << " ";
		if (i > 0 && i < allocPara_m.maxMembrNodePerCell
				&& divAuxData.tmpIsActive1_M[i]
				&& divAuxData.tmpIsActive1_M[i - 1]
				&& fabs(divAuxData.tmpXPos1_M[i] - divAuxData.tmpXPos1_M[i - 1])
						> 0.1) {
			std::cout << "11111111111111111111111, " << i << std::endl;
			int jj;
			cin >> jj;
		}
	}
	std::cout << std::endl;
	std::cout << "XPos1_onDevice : ";
	for (uint i = 0; i < divAuxData.tmpCellRank_M.size(); i++) {
		for (uint j = 0; j < allocPara_m.maxAllNodePerCell; j++) {
			uint index = divAuxData.tmpCellRank_M[i]
					* allocPara_m.maxAllNodePerCell + j;
			std::cout << nodes->getInfoVecs().nodeLocX[index] << " ";
		}
	}
	std::cout << std::endl;

	std::cout << "tmpYPos1_M : ";
	for (uint i = 0; i < divAuxData.tmpYPos1_M.size(); i++) {
		std::cout << divAuxData.tmpYPos1_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpIsActive2_M: ";
	for (uint i = 0; i < divAuxData.tmpIsActive2_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive2_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpXPos2_M : ";
	for (uint i = 0; i < divAuxData.tmpXPos2_M.size(); i++) {
		std::cout << divAuxData.tmpXPos2_M[i] << " ";
		if (i > 0 && i < allocPara_m.maxMembrNodePerCell
				&& divAuxData.tmpIsActive2_M[i]
				&& divAuxData.tmpIsActive2_M[i - 1]
				&& fabs(divAuxData.tmpXPos2_M[i] - divAuxData.tmpXPos2_M[i - 1])
						> 0.1) {
			std::cout << "2222222222222222222, " << i << std::endl;
			int jj;
			cin >> jj;
		}
	}
	std::cout << std::endl;
	std::cout << "tmpYPos2_M : ";
	for (uint i = 0; i < divAuxData.tmpYPos2_M.size(); i++) {
		std::cout << divAuxData.tmpYPos2_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp1InternalActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp1InternalActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp1InternalActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp2InternalActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp2InternalActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp2InternalActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp1MemActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp1MemActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp1MemActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp2MemActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp2MemActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp2MemActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	int jj;
	std::cin >> jj;
}

void SceCells::adjustGrowthInfo_M() {
	uint halfMax = allocPara_m.maxIntnlNodePerCell / 2;
	thrust::transform_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			AdjustGrowth(halfMax), thrust::identity<bool>());
}

VtkAnimationData SceCells::outputVtkData(AniRawData& rawAniData,
		AnimationCriteria& aniCri) {
	VtkAnimationData vtkData;
	for (uint i = 0; i < rawAniData.aniNodePosArr.size(); i++) {
		PointAniData ptAniData;
		ptAniData.pos = rawAniData.aniNodePosArr[i];
		ptAniData.F_MI_M_MagN_Int= rawAniData.aniNodeF_MI_M_MagN_Int[i]; //AliE
		ptAniData.F_MI_M = rawAniData.aniNodeF_MI_M[i];//AAMIRI
		ptAniData.colorScale = rawAniData.aniNodeVal[i];
		ptAniData.colorScale2 = rawAniData.aniNodeCurvature[i];//AAMIRI
		ptAniData.colorScale3 = rawAniData.aniNodeMembTension[i];//Ali 
		//ptAniData.colorScale4 = rawAniData.aniNodeVel_Mag[i];//rawAniData.aniNodeActinLevel[i];//Ali 
		ptAniData.colorScale4 = rawAniData.aniNodeContrApi[i];
		ptAniData.rankScale = rawAniData.aniNodeRank[i];//AAMIRI
		ptAniData.intercellForce = rawAniData.aniNodeInterCellForceArr[i];//AAMIRI
		vtkData.pointsAniData.push_back(ptAniData);
	}
	for (uint i = 0; i < rawAniData.internalLinks.size(); i++) {
		LinkAniData linkData = rawAniData.internalLinks[i];
		vtkData.linksAniData.push_back(linkData);
	}
	for (uint i = 0; i < rawAniData.memLinks.size(); i++) {
		LinkAniData linkData = rawAniData.memLinks[i];
		vtkData.linksAniData.push_back(linkData);
	}
	vtkData.isArrowIncluded = false;
	return vtkData;
}

void SceCells::copyToGPUConstMem() {
	double pI_CPU = acos(-1.0);
	double minLengthCPU =
			globalConfigVars.getConfigValue("MinLength").toDouble();
	cudaMemcpyToSymbol(minLength, &minLengthCPU, sizeof(double));
	double minDivisorCPU =
			globalConfigVars.getConfigValue("MinDivisor").toDouble();
	cudaMemcpyToSymbol(minDivisor, &minDivisorCPU, sizeof(double));
	cudaMemcpyToSymbol(membrEquLen, &membrPara.membrEquLenCPU, sizeof(double));
	cudaMemcpyToSymbol(membrStiff, &membrPara.membrStiffCPU, sizeof(double));
	cudaMemcpyToSymbol(membrStiff_Mitotic, &membrPara.membrStiff_Mitotic, sizeof(double)); // Ali June 30
	cudaMemcpyToSymbol(kContractMemb, &membrPara.kContractMemb, sizeof(double)); 
	cudaMemcpyToSymbol(pI, &pI_CPU, sizeof(double));

	cudaMemcpyToSymbol(bendCoeff, &membrPara.membrBendCoeff, sizeof(double));

	cudaMemcpyToSymbol(bendCoeff_Mitotic, &membrPara.membrBendCoeff_Mitotic, sizeof(double));//AAMIRI

	cudaMemcpyToSymbol(F_Ext_Incline_M2, &membrPara.F_Ext_Incline, sizeof(double)); //Ali
      
	uint maxAllNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxAllNodeCountPerCell").toInt();
	uint maxMembrNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxMembrNodeCountPerCell").toInt();
	uint maxIntnlNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();

	cudaMemcpyToSymbol(maxAllNodePerCell, &maxAllNodePerCellCPU, sizeof(uint));
	cudaMemcpyToSymbol(maxMembrPerCell, &maxMembrNodePerCellCPU, sizeof(uint));
	cudaMemcpyToSymbol(maxIntnlPerCell, &maxIntnlNodePerCellCPU, sizeof(uint));

	double sceIntnlBParaCPU_M[5];
	double sceIntraParaCPU_M[5];
	double sceIntraParaDivCPU_M[5];

	double U0_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_U0").toDouble();
	double V0_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_V0").toDouble();
	double k1_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_k1").toDouble();
	double k2_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_k2").toDouble();
	double intnlBEffectiveRange = globalConfigVars.getConfigValue(
			"IntnlBEffectRange").toDouble();
	sceIntnlBParaCPU_M[0] = U0_IntnlB;
	sceIntnlBParaCPU_M[1] = V0_IntnlB;
	sceIntnlBParaCPU_M[2] = k1_IntnlB;
	sceIntnlBParaCPU_M[3] = k2_IntnlB;
	sceIntnlBParaCPU_M[4] = intnlBEffectiveRange;


        
 
	//////////////////////
	//// Block 3 /////////
	//////////////////////
	double U0_Intra =
			globalConfigVars.getConfigValue("IntraCell_U0").toDouble();
	double V0_Intra =
			globalConfigVars.getConfigValue("IntraCell_V0").toDouble();
	double k1_Intra =
			globalConfigVars.getConfigValue("IntraCell_k1").toDouble();
	double k2_Intra =
			globalConfigVars.getConfigValue("IntraCell_k2").toDouble();
	double intraLinkEffectiveRange = globalConfigVars.getConfigValue(
			"IntraEffectRange").toDouble();
	sceIntraParaCPU_M[0] = U0_Intra;
	sceIntraParaCPU_M[1] = V0_Intra;
	sceIntraParaCPU_M[2] = k1_Intra;
	sceIntraParaCPU_M[3] = k2_Intra;
	sceIntraParaCPU_M[4] = intraLinkEffectiveRange;

	//////////////////////
	//// Block 4 /////////
	//////////////////////
	double U0_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_U0_Div").toDouble();
	double V0_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_V0_Div").toDouble();
	double k1_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_k1_Div").toDouble();
	double k2_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_k2_Div").toDouble();
	double intraDivEffectiveRange = globalConfigVars.getConfigValue(
			"IntraDivEffectRange").toDouble();
	sceIntraParaDivCPU_M[0] = U0_Intra_Div;
	sceIntraParaDivCPU_M[1] = V0_Intra_Div;
	sceIntraParaDivCPU_M[2] = k1_Intra_Div;
	sceIntraParaDivCPU_M[3] = k2_Intra_Div;
	sceIntraParaDivCPU_M[4] = intraDivEffectiveRange;

	cudaMemcpyToSymbol(grthPrgrCriEnd_M, &growthAuxData.grthProgrEndCPU,
			sizeof(double));
	//cudaMemcpyToSymbol(grthPrgrCriVal_M, &growthPrgrCriVal, sizeof(double));
	cudaMemcpyToSymbol(sceIB_M, sceIntnlBParaCPU_M, 5 * sizeof(double));
	cudaMemcpyToSymbol(sceII_M, sceIntraParaCPU_M, 5 * sizeof(double));
	cudaMemcpyToSymbol(sceIIDiv_M, sceIntraParaDivCPU_M, 5 * sizeof(double));

	
	double IBDivHost[5];
	IBDivHost[0] =
			globalConfigVars.getConfigValue("SceIntnlB_U0_Div").toDouble();
	IBDivHost[1] =
			globalConfigVars.getConfigValue("SceIntnlB_V0_Div").toDouble();
	IBDivHost[2] =
			globalConfigVars.getConfigValue("SceIntnlB_k1_Div").toDouble();
	IBDivHost[3] =
			globalConfigVars.getConfigValue("SceIntnlB_k2_Div").toDouble();
	IBDivHost[4] =
			globalConfigVars.getConfigValue("IntnlBDivEffectRange").toDouble();
	cudaMemcpyToSymbol(sceIBDiv_M, IBDivHost, 5 * sizeof(double));


     //////////////////////
	//// Block Nucleus  /////////
	//////////////////////


	double sceNucleusParaCPU_M[5];

	double U0_Nucleus =
			globalConfigVars.getConfigValue("NucleusCell_U0").toDouble();
	double V0_Nucleus =
			globalConfigVars.getConfigValue("NucleusCell_V0").toDouble();
	double k1_Nucleus =
			globalConfigVars.getConfigValue("NucleusCell_k1").toDouble();
	double k2_Nucleus =
			globalConfigVars.getConfigValue("NucleusCell_k2").toDouble();
	double nucleusLinkEffectiveRange = globalConfigVars.getConfigValue(
			"NucleusEffectRange").toDouble();
	sceNucleusParaCPU_M[0] = U0_Nucleus;
	sceNucleusParaCPU_M[1] = V0_Nucleus;
	sceNucleusParaCPU_M[2] = k1_Nucleus;
	sceNucleusParaCPU_M[3] = k2_Nucleus;
	sceNucleusParaCPU_M[4] = nucleusLinkEffectiveRange;

	//////////////////////
	//// Block Nucleus Division /////////
	//////////////////////


	double sceNucleusParaDivCPU_M[5];

	double U0_Nucleus_Div =
			globalConfigVars.getConfigValue("NucleusCell_U0_Div").toDouble();
	double V0_Nucleus_Div =
			globalConfigVars.getConfigValue("NucleusCell_V0_Div").toDouble();
	double k1_Nucleus_Div =
			globalConfigVars.getConfigValue("NucleusCell_k1_Div").toDouble();
	double k2_Nucleus_Div =
			globalConfigVars.getConfigValue("NucleusCell_k2_Div").toDouble();
	double nucleusDivEffectiveRange = globalConfigVars.getConfigValue(
			"NucleusDivEffectRange").toDouble();
	sceNucleusParaDivCPU_M[0] = U0_Nucleus_Div;
	sceNucleusParaDivCPU_M[1] = V0_Nucleus_Div;
	sceNucleusParaDivCPU_M[2] = k1_Nucleus_Div;
	sceNucleusParaDivCPU_M[3] = k2_Nucleus_Div;
	sceNucleusParaDivCPU_M[4] = nucleusDivEffectiveRange;


	cudaMemcpyToSymbol(sceN_M,    sceNucleusParaCPU_M,    5 * sizeof(double));  //Ali 
	cudaMemcpyToSymbol(sceNDiv_M, sceNucleusParaDivCPU_M, 5 * sizeof(double)); //Ali
	

}

void SceCells::updateMembrGrowthProgress_M() {

	// figure out membr growth speed
	calMembrGrowSpeed_M();  //Ali: to my understanding it doesn't do anything right now. it will be override by adjustMembrGrowSpeed_M 
	// figure out which cells will add new point and which cell needs to delete node.

	// adjustMembrGrowSpeed_M(); // for now just a constant speed to give some relaxation before adding another node.

	// returning a bool and progress for each cell. if bool is true (a node sould be added) progress will be reset to give relaxation time after adding the node. Otherwise growth prgoress will be incremented
// add membr nodes  // In each time step either adding mechanism is active or deleting mechanism. It is an unneccessary complication to manage memory for both operations at one time step.

	// uint curActCellCt = allocPara_m.currentActiveCellCount;
	
	// thrust::transform(cellInfoVecs.membrGrowSpeed.begin(),
	// 		cellInfoVecs.membrGrowSpeed.begin() + curActCellCt,
	// 		cellInfoVecs.membrGrowProgress.begin(),
	// 		cellInfoVecs.membrGrowProgress.begin(), SaxpyFunctor(dt));





}
void SceCells::handleMembrGrowth_M(int maxApicalBasalNodeNum, double maxLengthToAddMemNodes) {
	
if (1>0){//addNode) {

		decideIfAddMembrNode_M(maxApicalBasalNodeNum, maxLengthToAddMemNodes); 
		addMembrNodes_M();
		// addNode=false  ; 
		// cout << " I am in add membrane node " << endl ; 
	}	
    else  {

		// decideIfDelMembrNode_M(); //Ali 
		// delMembrNodes_M();
		// addNode=true ; 
		// cout << " I am in del membrane node " << endl ; 
		}
	//membrDebug();
}





// void SceCells::calMembrGrowSpeed_M() {
// 	membrPara.membrGrowCoeff = growthAuxData.prolifDecay
// 			* membrPara.membrGrowCoeff_Ori;
// 	membrPara.membrGrowLimit = growthAuxData.prolifDecay
// 			* membrPara.membrGrowLimit_Ori;
// // reduce_by_key, find value of max tension and their index
// 	thrust::counting_iterator<uint> iBegin(0);
// 	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	
// 	thrust::reduce_by_key(
// 			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
// 			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
// 					+ totalNodeCountForActiveCells,
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().membrTenMagRi.begin(),
// 							make_transform_iterator(iBegin,
// 									ModuloFunctor(maxNPerCell)),
// 							nodes->getInfoVecs().membrLinkRiMidX.begin(),
// 							nodes->getInfoVecs().membrLinkRiMidY.begin(),
// 							nodes->getInfoVecs().membrDistToRi.begin())),
// 			cellInfoVecs.cellRanksTmpStorage.begin(),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(cellInfoVecs.maxTenRiVec.begin(),
// 							cellInfoVecs.maxTenIndxVec.begin(),
// 							cellInfoVecs.maxTenRiMidXVec.begin(),
// 							cellInfoVecs.maxTenRiMidYVec.begin(),
// 							cellInfoVecs.maxDistToRiVec.begin())),
// 			thrust::equal_to<uint>(), MaxWInfo());

// //	for (int i=0 ; i<cellInfoVecs.maxDistToRiVec.size() ; i++) {
// //		cout << "the max distance in cell" << i << " is "<<cellInfoVecs.maxDistToRiVec[i] << endl ; 
// //	}

// 	//Ali for min Distance

// 	thrust::counting_iterator<uint> iBegin_min(0);
// thrust::reduce_by_key(
// 			make_transform_iterator(iBegin_min, DivideFunctor(maxNPerCell)), // begin of the key 
// 			make_transform_iterator(iBegin_min, DivideFunctor(maxNPerCell))  // end of the key 
// 					+ totalNodeCountForActiveCells,
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							nodes->getInfoVecs().membrDistToRi.begin(),
// 							make_transform_iterator(iBegin_min,   // values to reduce by key 
// 									ModuloFunctor(maxNPerCell))  
// 							)),
// 			cellInfoVecs.cellRanksTmpStorage1.begin(),  // to Store reduced version of key 
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							cellInfoVecs.minDistToRiVec.begin(),
// 							cellInfoVecs.minTenIndxVec.begin()  // to sotred the reduce verision of values 
// 							)), 
// 			thrust::equal_to<uint>(), MinWInfo());  // how to sort the keys & how to reduce the parameters assigned to based on each key
// // equal_to mean how we set the beans to reduce. For example here we are saying if they are equal in Int we compare them and would peroform the reduction.

// //	for (int i=0 ; i<cellInfoVecs.minDistToRiVec.size() ; i++) {
// //		cout << "the min distance in cell" << i << " is "<<cellInfoVecs.minDistToRiVec[i] << endl ; 
// //		cout << "the min tension index vec" << i << " is "<<cellInfoVecs.minTenIndxVec[i] << endl ; 
// //	}


// 	thrust::reduce_by_key(
// 			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
// 			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
// 					+ totalNodeCountForActiveCells,
// 			nodes->getInfoVecs().membrTensionMag.begin(),
// 			cellInfoVecs.cellRanksTmpStorage.begin(),
// 			cellInfoVecs.aveTension.begin(), thrust::equal_to<uint>(),
// 			thrust::plus<double>());

// 	thrust::transform(cellInfoVecs.aveTension.begin(),
// 			cellInfoVecs.aveTension.begin()
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.activeMembrNodeCounts.begin(),
// 			cellInfoVecs.aveTension.begin(), thrust::divides<double>());

// 	// linear relationship with highest tension; capped by a given value
// 	thrust::transform(cellInfoVecs.aveTension.begin(),
// 			cellInfoVecs.aveTension.begin()
// 					+ allocPara_m.currentActiveCellCount,
// 			cellInfoVecs.membrGrowSpeed.begin(),
// 			MultiWithLimit(membrPara.membrGrowCoeff, membrPara.membrGrowLimit));
// }

void SceCells::calMembrGrowSpeed_M() {
	membrPara.membrGrowCoeff = growthAuxData.prolifDecay
			* membrPara.membrGrowCoeff_Ori;
	membrPara.membrGrowLimit = growthAuxData.prolifDecay
			* membrPara.membrGrowLimit_Ori;
// reduce_by_key, find value of max tension and their index
	thrust::counting_iterator<uint> iBegin(0);
	
	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	
	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().membrTenMagRi.begin(),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxNPerCell)),
							nodes->getInfoVecs().membrLinkRiMidX.begin(),
							nodes->getInfoVecs().membrLinkRiMidY.begin(),
							nodes->getInfoVecs().membrDistToRi.begin(),
							nodes->getInfoVecs().memNodeType1.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.maxTenRiVec.begin(),
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin(),
							cellInfoVecs.maxDistToRiVec.begin(),
							cellInfoVecs.maxTenIndxTypeVec.begin())),
			thrust::equal_to<uint>(), MaxWInfo());

	// for (int i=0 ; i<cellInfoVecs.maxDistToRiVec.size() ; i++) {
	// 	cout << "the max distance in cell" << i << " is "<<cellInfoVecs.maxDistToRiVec[i] << endl ;
	// 	cout << "At index "<<cellInfoVecs.maxTenIndxVec[i]<<std::endl; 
	// }

	//Ali for min Distance

	thrust::counting_iterator<uint> iBegin_min(0);
thrust::reduce_by_key(
			make_transform_iterator(iBegin_min, DivideFunctor(maxNPerCell)), // begin of the key 
			make_transform_iterator(iBegin_min, DivideFunctor(maxNPerCell))  // end of the key 
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().membrDistToRi.begin(),
							make_transform_iterator(iBegin_min,   // values to reduce by key 
									ModuloFunctor(maxNPerCell))  
							)),
			cellInfoVecs.cellRanksTmpStorage1.begin(),  // to Store reduced version of key 
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.minDistToRiVec.begin(),
							cellInfoVecs.minTenIndxVec.begin()  // to sotred the reduce verision of values 
							)), 
			thrust::equal_to<uint>(), MinWInfo());  // how to sort the keys & how to reduce the parameters assigned to based on each key
// equal_to mean how we set the beans to reduce. For example here we are saying if they are equal in Int we compare them and would peroform the reduction.

//	for (int i=0 ; i<cellInfoVecs.minDistToRiVec.size() ; i++) {
//		cout << "the min distance in cell" << i << " is "<<cellInfoVecs.minDistToRiVec[i] << endl ; 
//		cout << "the min tension index vec" << i << " is "<<cellInfoVecs.minTenIndxVec[i] << endl ; 
//	}


	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().membrTensionMag.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.aveTension.begin(), thrust::equal_to<uint>(),
			thrust::plus<double>());

	thrust::transform(cellInfoVecs.aveTension.begin(),
			cellInfoVecs.aveTension.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.activeMembrNodeCounts.begin(),
			cellInfoVecs.aveTension.begin(), thrust::divides<double>());

	// linear relationship with highest tension; capped by a given value
	thrust::transform(cellInfoVecs.aveTension.begin(),
			cellInfoVecs.aveTension.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.membrGrowSpeed.begin(),
			MultiWithLimit(membrPara.membrGrowCoeff, membrPara.membrGrowLimit));
}

void SceCells::adjustMembrGrowSpeed_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.membrGrowSpeed.begin(),
			AdjustMembrGrow(membrPara.growthConst_N, membrPara.initMembrCt_N,
					membrPara.initIntnlCt_N));
}

void SceCells::decideIfAddMembrNode_M(int maxApicalBasalNodeNum, double maxLengthToAddMemNodes) {
// decide if add membrane node given current active node count and
// membr growth progresss
	uint curActCellCt = allocPara_m.currentActiveCellCount;
	uint maxMembrNode = allocPara_m.maxMembrNodePerCell;
	bool isInitPhase= nodes->isInitPhase ; 
	/*
	thrust::transform(cellInfoVecs.membrGrowSpeed.begin(),
			cellInfoVecs.membrGrowSpeed.begin() + curActCellCt,
			cellInfoVecs.membrGrowProgress.begin(),
			cellInfoVecs.membrGrowProgress.begin(), SaxpyFunctor(dt));
*/

/*
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.membrGrowProgress.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.membrGrowProgress.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin()))
					+ curActCellCt,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isMembrAddingNode.begin(),
							cellInfoVecs.membrGrowProgress.begin())),
			MemGrowFunc(maxMembrNode));
*/
         thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							           cellInfoVecs.membrGrowProgress.begin(),
									   cellInfoVecs.maxDistToRiVec.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.numApicalVec.begin(),
									   cellInfoVecs.numBasalVec.begin(),
									   cellInfoVecs.maxTenIndxTypeVec.begin(),
									   cellInfoVecs.cellRankVec.begin()
									   )),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							           cellInfoVecs.membrGrowProgress.begin(),
									   cellInfoVecs.maxDistToRiVec.begin(),
									   cellInfoVecs.isEnteringMitotic.begin(),
									   cellInfoVecs.numApicalVec.begin(),
									   cellInfoVecs.numBasalVec.begin(),
									   cellInfoVecs.maxTenIndxTypeVec.begin(),
									   cellInfoVecs.cellRankVec.begin()
									   ))
					+ curActCellCt,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isMembrAddingNode.begin(),
							cellInfoVecs.membrGrowProgress.begin(),
							cellInfoVecs.numApicalVec.begin(),
							cellInfoVecs.numBasalVec.begin())),
			MemGrowFunc(maxMembrNode,isInitPhase, maxApicalBasalNodeNum, maxLengthToAddMemNodes));

}
//Ali
void SceCells::decideIfDelMembrNode_M() {
	uint curActCellCt = allocPara_m.currentActiveCellCount;
		uint maxMembrNode = allocPara_m.maxMembrNodePerCell;

	bool isInitPhase= nodes->isInitPhase ; 
         thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							   		   cellInfoVecs.membrGrowProgress.begin(),
									   cellInfoVecs.minDistToRiVec.begin(),
									  cellInfoVecs.growthProgress.begin()
									  )),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							   		   cellInfoVecs.membrGrowProgress.begin(),
									   cellInfoVecs.minDistToRiVec.begin(),
									  cellInfoVecs.growthProgress.begin()
									  ))
					+ curActCellCt,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isMembrRemovingNode.begin(),
							cellInfoVecs.membrGrowProgress.begin())),
			MemDelFunc(maxMembrNode, isInitPhase));

}


/**
 * Add new membrane elements to cells.
 * This operation is relatively expensive because of memory rearrangement.
 */
// void SceCells::addMembrNodes_M() {
// 	thrust::counting_iterator<uint> iBegin(0);
// 	uint curAcCCount = allocPara_m.currentActiveCellCount;
// 	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
// 	bool* nodeIsActiveAddress = thrust::raw_pointer_cast(
// 					&(nodes->getInfoVecs().nodeIsActive[0]));
// 	double* nodeXPosAddress = thrust::raw_pointer_cast(
// 					&(nodes->getInfoVecs().nodeLocX[0]));
// 	double* nodeYPosAddress = thrust::raw_pointer_cast(
// 					&(nodes->getInfoVecs().nodeLocY[0]));
// 	int* adhIndxAddr = thrust::raw_pointer_cast(
// 					&(nodes->getInfoVecs().nodeAdhereIndex[0]));
// 	MembraneType1* memNodeType1 = thrust::raw_pointer_cast(
// 					&(nodes->getInfoVecs().memNodeType1[0]));
// 	thrust::transform_if(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin,
// 							cellInfoVecs.maxTenIndxVec.begin(),
// 							cellInfoVecs.activeMembrNodeCounts.begin(),
// 							cellInfoVecs.maxTenRiMidXVec.begin(),
// 							cellInfoVecs.maxTenRiMidYVec.begin(),
// 							cellInfoVecs.ringApicalId.begin(),
// 							cellInfoVecs.ringBasalId.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(iBegin,
// 							cellInfoVecs.maxTenIndxVec.begin(),
// 							cellInfoVecs.activeMembrNodeCounts.begin(),
// 							cellInfoVecs.maxTenRiMidXVec.begin(),
// 							cellInfoVecs.maxTenRiMidYVec.begin(),
// 							cellInfoVecs.ringApicalId.begin(),
// 							cellInfoVecs.ringBasalId.begin()))
// 					+ curAcCCount, cellInfoVecs.isMembrAddingNode.begin(),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							cellInfoVecs.activeMembrNodeCounts.begin(),
// 							cellInfoVecs.ringApicalId.begin(),
// 							cellInfoVecs.ringBasalId.begin())),
// 			AddMemNode(maxNodePerCell, 
// 					nodeIsActiveAddress,
// 					nodeXPosAddress,
// 					nodeYPosAddress, 
// 					adhIndxAddr, 
// 					memNodeType1),
// 			thrust::identity<bool>());
// }

void SceCells::addMembrNodes_M() {
	thrust::counting_iterator<uint> iBegin(0);
	uint curAcCCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	bool* nodeIsActiveAddress = thrust::raw_pointer_cast(
					&(nodes->getInfoVecs().nodeIsActive[0]));
	double* nodeXPosAddress = thrust::raw_pointer_cast(
					&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeYPosAddress = thrust::raw_pointer_cast(
					&(nodes->getInfoVecs().nodeLocY[0]));
	int* adhIndxAddr = thrust::raw_pointer_cast(
					&(nodes->getInfoVecs().nodeAdhereIndex[0]));
	MembraneType1* memNodeType1 = thrust::raw_pointer_cast(
					&(nodes->getInfoVecs().memNodeType1[0]));
	thrust::transform_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin(),
							cellInfoVecs.ringApicalId.begin(),
							cellInfoVecs.ringBasalId.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin(),
							cellInfoVecs.ringApicalId.begin(),
							cellInfoVecs.ringBasalId.begin()))
					+ curAcCCount, cellInfoVecs.isMembrAddingNode.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
			AddMemNode(maxNodePerCell, 
					nodeIsActiveAddress,
					nodeXPosAddress,
					nodeYPosAddress, 
					adhIndxAddr, 
					memNodeType1),
			thrust::identity<bool>());
			for (int z = 0; z < allocPara_m.currentActiveCellCount; z++){
				if (cellInfoVecs.isMembrAddingNode[z] == true){
					nodes->isMemNodeTypeAssigned_postAddNode=false;
					cellInfoVecs.isPostAddMembrNodes = true;
					break;
				}
			}
}

//Ali
void SceCells::delMembrNodes_M() {
	thrust::counting_iterator<uint> iBegin(0);
	uint curAcCCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::transform_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.minTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.ringApicalId.begin(),
							cellInfoVecs.ringBasalId.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.ringApicalId.begin(),
							cellInfoVecs.ringBasalId.begin()))
					+ curAcCCount, cellInfoVecs.isMembrRemovingNode.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.ringApicalId.begin(),
							cellInfoVecs.ringBasalId.begin())),
			DelMemNode(maxNodePerCell, growthAuxData.nodeIsActiveAddress,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress, growthAuxData.adhIndxAddr,growthAuxData.memNodeType1Address),
			thrust::identity<bool>());
}

void SceCells::membrDebug() {
	uint curAcCCount = allocPara_m.currentActiveCellCount;
	uint maxActiveNodeC = curAcCCount * allocPara_m.maxAllNodePerCell;
	uint maxNodePC = allocPara_m.maxAllNodePerCell;
	//uint tmp = 0;
	//for (uint i = 0; i < curAcCCount; i++) {
	//	tmp += cellInfoVecs.isMembrAddingNode[i];
	//}
	//if (tmp != 0) {
	//	tmpDebug = true;
	//}
	//if (!tmpDebug) {
	//	return;
	//}
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrTensionMag[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrTenMagRi[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrLinkRiMidX[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrLinkRiMidY[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendLeftX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendLeftY[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendRightX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendRightX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < curAcCCount; i++) {
		std::cout << "(" << cellInfoVecs.maxTenIndxVec[i] << ","
				<< cellInfoVecs.activeMembrNodeCounts[i] << ","
				<< cellInfoVecs.maxTenRiMidXVec[i] << ","
				<< cellInfoVecs.maxTenRiMidYVec[i] << ")" << std::endl;
	}
	int jj;
	std::cin >> jj;
}

void SceCells::assembleVecForTwoCells(uint i) {
	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	for (uint j = 0; j < membThreshold; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < divAuxData.tmp1VecMem.size()) {
			divAuxData.tmpXPos1_M[index] = divAuxData.tmp1VecMem[j].x;
			divAuxData.tmpYPos1_M[index] = divAuxData.tmp1VecMem[j].y;
			divAuxData.tmpNodeType1[index] = divAuxData.tmp1VecMemNodeType[j] ; //Ali 
				// std::cout<<"divAuxData.tmpPos1_M["<<index<<"] = "<<divAuxData.tmpXPos1_M[index]<<" "<<divAuxData.tmpYPos1_M[index]<<", tmpNodeType1 = "<<divAuxData.tmpNodeType1[index]<<std::endl;
			divAuxData.tmpIsActive1_M[index] = true;
		} else {
			divAuxData.tmpIsActive1_M[index] = false;
		}
	}
	for (uint j = 0; j < membThreshold; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < divAuxData.tmp2VecMem.size()) {
			divAuxData.tmpXPos2_M[index] = divAuxData.tmp2VecMem[j].x;
			divAuxData.tmpYPos2_M[index] = divAuxData.tmp2VecMem[j].y;
			divAuxData.tmpNodeType2[index] = divAuxData.tmp2VecMemNodeType[j] ; //Ali 
				// std::cout<<"divAuxData.tmpPos2_M["<<index<<"] = "<<divAuxData.tmpXPos2_M[index]<<" "<<divAuxData.tmpYPos2_M[index]<<", tmpNodeType2 = "<<divAuxData.tmpNodeType2[index]<<std::endl;
			divAuxData.tmpIsActive2_M[index] = true;
		} else {
			divAuxData.tmpIsActive2_M[index] = false;
		}
	}

	divAuxData.tmp1MemActiveCounts.push_back(divAuxData.tmp1VecMem.size());
	std::cout<<"divAuxData.tmp1MemActiveCounts size = "<<divAuxData.tmp1MemActiveCounts[0]<<std::endl;
	divAuxData.tmp2MemActiveCounts.push_back(divAuxData.tmp2VecMem.size());
	std::cout<<"divAuxData.tmp2MemActiveCounts size = "<<divAuxData.tmp2MemActiveCounts[0]<<std::endl;

	for (uint j = membThreshold; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		uint shift_j = j - membThreshold;
		if (shift_j < divAuxData.tmp1IntnlVec.size()) {
			divAuxData.tmpXPos1_M[index] = divAuxData.tmp1IntnlVec[shift_j].x;
			divAuxData.tmpYPos1_M[index] = divAuxData.tmp1IntnlVec[shift_j].y;
			divAuxData.tmpNodeType1[index] = notAssigned1 ; //Ali 
			divAuxData.tmpIsActive1_M[index] = true;
		} else {
			divAuxData.tmpIsActive1_M[index] = false;
		}
		if (shift_j < divAuxData.tmp2IntnlVec.size()) {
			divAuxData.tmpXPos2_M[index] = divAuxData.tmp2IntnlVec[shift_j].x;
			divAuxData.tmpYPos2_M[index] = divAuxData.tmp2IntnlVec[shift_j].y;
			divAuxData.tmpNodeType2[index] = notAssigned1 ; //Ali 
			divAuxData.tmpIsActive2_M[index] = true;
		} else {
			divAuxData.tmpIsActive2_M[index] = false;
		}
	}
	divAuxData.tmp1InternalActiveCounts.push_back(
			divAuxData.tmp1IntnlVec.size());
	divAuxData.tmp2InternalActiveCounts.push_back(
			divAuxData.tmp2IntnlVec.size());
}


// we have two new center of internal node positions. 
// we already shrinked the internal nodes around their old internal nodes center
// here we shift the internal nodes of each cell around the new internal node position
void SceCells::shiftIntnlNodesByCellCenter(CVector intCell1Center,
		CVector intCell2Center) {
	CVector tmpCell1Center(0, 0, 0);
	for (uint j = 0; j < divAuxData.tmp1IntnlVec.size(); j++) {
		tmpCell1Center = tmpCell1Center + divAuxData.tmp1IntnlVec[j];
	}
	tmpCell1Center = tmpCell1Center / divAuxData.tmp1IntnlVec.size();
	CVector shiftVec1 = intCell1Center - tmpCell1Center; // it should be new nucleus center for cell1 
	for (uint j = 0; j < divAuxData.tmp1IntnlVec.size(); j++) {
		divAuxData.tmp1IntnlVec[j] = divAuxData.tmp1IntnlVec[j] + shiftVec1;
	}

	CVector tmpCell2Center(0, 0, 0);
	for (uint j = 0; j < divAuxData.tmp2IntnlVec.size(); j++) {
		tmpCell2Center = tmpCell2Center + divAuxData.tmp2IntnlVec[j];
	}
	tmpCell2Center = tmpCell2Center / divAuxData.tmp2IntnlVec.size();
	CVector shiftVec2 = intCell2Center - tmpCell2Center; // it should be new nucleus center for cell 2 
	for (uint j = 0; j < divAuxData.tmp2IntnlVec.size(); j++) {
		divAuxData.tmp2IntnlVec[j] = divAuxData.tmp2IntnlVec[j] + shiftVec2;
	}
}

void SceCells::processMemVec(uint i, std::vector<VecValT>& tmp1,
		std::vector<VecValT>& tmp2, CVector oldNucleusCenter) {
	divAuxData.tmp1VecMem.clear();
	divAuxData.tmp2VecMem.clear();
	divAuxData.tmp1VecMemNodeType.clear(); //Ali
	divAuxData.tmp2VecMemNodeType.clear(); //Ali

	uint membThreshold = allocPara_m.maxMembrNodePerCell;

	std::sort(tmp1.begin(), tmp1.end());
	std::sort(tmp2.begin(), tmp2.end());
	// Here we perform a cross-product computation to check the node order orientation. If it is clockwise, reverse it so it is counter-clockwise.
	double vec1x = (tmp1[0].vec.x - oldNucleusCenter.x);
	double vec2x = (tmp1[1].vec.x - oldNucleusCenter.x);
	double vec1y = (tmp1[0].vec.y - oldNucleusCenter.y);
	double vec2y = (tmp1[1].vec.y - oldNucleusCenter.y);
	double orientation = vec1x*vec2y - vec2x*vec1y;
	if (orientation < 0){
		std::cout<<"tmp1 orientation incorrect. Rearranging tmp1."<<std::endl;
		std::vector<VecValT> tmp_tmp1 = tmp1;
		tmp1.clear();
		for (int i = 0; i < tmp_tmp1.size(); i++){
			tmp1.push_back(tmp_tmp1[(tmp_tmp1.size()-1) - i]);
		}
	}
	vec1x = (tmp2[0].vec.x - oldNucleusCenter.x);
	vec2x = (tmp2[1].vec.x - oldNucleusCenter.x);
	vec1y = (tmp2[0].vec.y - oldNucleusCenter.y);
	vec2y = (tmp2[1].vec.y - oldNucleusCenter.y);
	orientation = vec1x*vec2y - vec2x*vec1y;
	if (orientation < 0){
		std::cout<<"tmp2 orientation incorrect. Rearranging tmp2."<<std::endl;
		std::vector<VecValT> tmp_tmp2 = tmp2;
		tmp2.clear();
		for (int i = 0; i < tmp_tmp2.size(); i++){
			tmp2.push_back(tmp_tmp2[(tmp_tmp2.size()-1) - i]);
		}
	}

	// std::vector<VecValT> sorted_tmp1 = tmp1;
	// std::vector<VecValT> sorted_tmp2 = tmp2;
	// std::sort(sorted_tmp1.begin(), sorted_tmp1.end());
	// std::sort(sorted_tmp2.begin(), sorted_tmp2.end());
	// for (int j = 0; j < tmp1.size(); j++){
	// 	std::cout<<"tmp1["<<j<<"].val = "<<tmp1[j].val<<", vec = "<<tmp1[j].vec.x<<" "<<tmp1[j].vec.y<<std::endl;
	// }
	// for (int j = 0; j < tmp2.size(); j++){
	// 	std::cout<<"tmp2["<<j<<"].val = "<<tmp2[j].val<<", vec = "<<tmp2[j].vec.x<<" "<<tmp2[j].vec.y<<std::endl;
	// }

	//assert(tmp1.size() < allocPara_m.maxMembrNodePerCell);
	//assert(tmp2.size() < allocPara_m.maxMembrNodePerCell);

	uint maxDivMembrNodeCount1 = allocPara_m.maxMembrNodePerCell - tmp1.size();
	uint maxDivMembrNodeCount2 = allocPara_m.maxMembrNodePerCell - tmp2.size();

	std::vector<CVector> ptsBetween1, ptsBetween2;

	uint Num_of_NodeType0_MotherCell = 0, Num_of_NodeType0_DaughterCell = 0, Num_of_NodeType1_DaughterCell = 0, Num_of_NodeType1_MotherCell = 0;
	if (divAuxData.isMotherCellBehind[i] == true){
		// Recall that tmp1 will always be the mother cell data //Kevin
		for (int j = 0; j < tmp1.size(); j++){
			if (tmp1[j].type == lateralB){
				Num_of_NodeType0_MotherCell+=1;
			}
		}
		if (Num_of_NodeType0_MotherCell > maxDivMembrNodeCount1){
			std::cout<<"Too many new nodes are needed to be introduced for mother cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		for (int j = 0; j < tmp2.size(); j++){
			if (tmp2[j].type == lateralA){
				Num_of_NodeType1_DaughterCell+=1;
			}
		}
		if (Num_of_NodeType1_DaughterCell > maxDivMembrNodeCount2){
			std::cout<<"Too many new nodes are needed to be introducedj for daughter cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		// Num_of_NodeType0_MotherCell = ;
		// Num_of_NodeType1_DaughterCell = ;
	}
	else{
		for (int j = 0; j < tmp1.size(); j++){
			if (tmp1[j].type == lateralA){
				Num_of_NodeType1_MotherCell+=1;
			}
		}
		if (Num_of_NodeType1_MotherCell > maxDivMembrNodeCount1){
			std::cout<<"Too many new nodes are needed to be introduced for mother cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		for (int j = 0; j < tmp2.size(); j++){
			if (tmp2[j].type == lateralB){
				Num_of_NodeType0_DaughterCell+=1;
			}
		}
		if (Num_of_NodeType0_DaughterCell > maxDivMembrNodeCount2){
			std::cout<<"Too many new nodes are needed to be introduced for daughter cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		
	}
	std::cout<<"Num_of_NodeType0_MotherCell = "<<Num_of_NodeType0_MotherCell<<std::endl;
	std::cout<<"Num_of_NodeType1_MotherCell = "<<Num_of_NodeType1_MotherCell<<std::endl;
	std::cout<<"Num_of_NodeType0_DaughterCell = "<<Num_of_NodeType0_DaughterCell<<std::endl;
	std::cout<<"Num_of_NodeType1_DaughterCell = "<<Num_of_NodeType1_DaughterCell<<std::endl;

	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp1.size() >= 1) {
		if (divAuxData.isMotherCellBehind[i] == true){
			ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
				Num_of_NodeType0_MotherCell, maxDivMembrNodeCount1);
			
		}
		else{
			ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
				Num_of_NodeType1_MotherCell, maxDivMembrNodeCount1);
			
		}
		// std::cout<<"tmp1[tmp1.size() - 1] = "<<tmp1[tmp1.size() - 1].vec.x<<" "<<tmp1[tmp1.size()-1].vec.y<<std::endl;
		// std::cout<<"tmp1[0] = "<<tmp1[0].vec.x<<" "<<tmp1[0].vec.y<<std::endl;
		// ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount1);
		// std::cout<<"sorted_tmp1[tmp1.size() - 1].val = "<<sorted_tmp1[sorted_tmp1.size()-1].val<<", .vec = "<< sorted_tmp1[sorted_tmp1.size() - 1].vec.x<<" "<<sorted_tmp1[sorted_tmp1.size()-1].vec.y<<std::endl;
		// std::cout<<"sorted_tmp1[0].val = "<<sorted_tmp1[0].val<<", .vec = "<<sorted_tmp1[0].vec.x<<" "<<sorted_tmp1[0].vec.y<<std::endl;
		// ptsBetween1 = obtainPtsBetween(sorted_tmp1[sorted_tmp1.size() - 1].vec, sorted_tmp1[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount1);
	}
	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp2.size() >= 1) {
		if (divAuxData.isMotherCellBehind[i] == true){
			ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
				Num_of_NodeType1_DaughterCell, maxDivMembrNodeCount2);
			
		}
		else{
			ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
				Num_of_NodeType0_DaughterCell, maxDivMembrNodeCount2);
			
		}
		// std::cout<<"tmp2[tmp1.size() - 1] = "<<tmp2[tmp2.size() - 1].vec.x<<" "<<tmp2[tmp2.size()-1].vec.y<<std::endl;
		// std::cout<<"tmp2[0] = "<<tmp2[0].vec.x<<" "<<tmp2[0].vec.y<<std::endl;
		// ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount2);
		// std::cout<<"sorted_tmp2[sorted_tmp2.size() - 1].val = "<<sorted_tmp2[sorted_tmp2.size()-1].val<<", .vec = "<<sorted_tmp2[sorted_tmp2.size() - 1].vec.x<<" "<<sorted_tmp2[sorted_tmp2.size()-1].vec.y<<std::endl;
		// std::cout<<"sorted_tmp2[0].val = "<<sorted_tmp2[0].val<<", vec = "<<sorted_tmp2[0].vec.x<<" "<<sorted_tmp2[0].vec.y<<std::endl;
		// ptsBetween2 = obtainPtsBetween(sorted_tmp2[sorted_tmp2.size() - 1].vec, sorted_tmp2[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount2);
	}

	for (uint j = 0; j < tmp1.size(); j++) {
		divAuxData.tmp1VecMem.push_back(tmp1[j].vec);
		divAuxData.tmp1VecMemNodeType.push_back(tmp1[j].type);
	}
	for (uint j = 0; j < tmp2.size(); j++) {
		divAuxData.tmp2VecMem.push_back(tmp2[j].vec);
		divAuxData.tmp2VecMemNodeType.push_back(tmp2[j].type);
	}
	std::cout<<"ptsBetween1 size = "<<ptsBetween1.size()<<std::endl;
	if (divAuxData.isMotherCellBehind[i] == false){
		for (uint j = 0; j < ptsBetween1.size(); j++) {
			divAuxData.tmp1VecMem.push_back(ptsBetween1[j]);
			// std::cout<<"ptsBtween1 "<<ptsBetween1[j].x<<" "<<ptsBetween1[j].y<<std::endl;
			divAuxData.tmp1VecMemNodeType.push_back(lateralB); 
		}
		// std::cout<<"size of tmp1VecMemNodeType = "<<divAuxData.tmp1VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp1VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp1VecMem["<<j<<"] = "<<divAuxData.tmp1VecMem[j].x<<" "<<divAuxData.tmp1VecMem[j].y<<" "<<divAuxData.tmp1VecMem[j].z<<", Type = "<<divAuxData.tmp1VecMemNodeType[j]<<std::endl;
		// }
	}
	else if (divAuxData.isMotherCellBehind[i] == true){
		for (uint j = 0; j < ptsBetween1.size(); j++) {
			divAuxData.tmp1VecMem.push_back(ptsBetween1[j]);
			// std::cout<<"ptsBtween1 "<<ptsBetween1[j].x<<" "<<ptsBetween1[j].y<<std::endl;
			divAuxData.tmp1VecMemNodeType.push_back(lateralA); 
		}
	// 	std::cout<<"size of tmp1VecMemNodeType = "<<divAuxData.tmp1VecMemNodeType.size()<<std::endl;
	// 	for (uint j = 0; j < divAuxData.tmp1VecMemNodeType.size();j++){
	// 		std::cout<<"divAuxData.tmp1VecMem["<<j<<"] = "<<divAuxData.tmp1VecMem[j].x<<" "<<divAuxData.tmp1VecMem[j].y<<" "<<divAuxData.tmp1VecMem[j].z<<", Type = "<<divAuxData.tmp1VecMemNodeType[j]<<std::endl;
	// 	}
	}
	std::cout<<"ptsBetween2 size = "<<ptsBetween2.size()<<std::endl;
	if (divAuxData.isMotherCellBehind[i] == false){
		for (uint j = 0; j < ptsBetween2.size(); j++) {
			divAuxData.tmp2VecMem.push_back(ptsBetween2[j]);
			// std::cout<<"ptsBtween2 "<<ptsBetween2[j].x<<" "<<ptsBetween2[j].y<<std::endl;
			divAuxData.tmp2VecMemNodeType.push_back(lateralA);
		}
		// std::cout<<"size of tmp2VecMemNodeType = "<<divAuxData.tmp2VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp2VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp2VecMem["<<j<<"] = "<<divAuxData.tmp2VecMem[j].x<<" "<<divAuxData.tmp2VecMem[j].y<<" "<<divAuxData.tmp2VecMem[j].z<<", Type = "<<divAuxData.tmp2VecMemNodeType[j]<<std::endl;
		// }
	}
	else if (divAuxData.isMotherCellBehind[i] == true){
		for (uint j = 0; j < ptsBetween2.size(); j++) {
			divAuxData.tmp2VecMem.push_back(ptsBetween2[j]);
			// std::cout<<"ptsBtween2 "<<ptsBetween2[j].x<<" "<<ptsBetween2[j].y<<std::endl;
			divAuxData.tmp2VecMemNodeType.push_back(lateralB);
		}
	// 	std::cout<<"size of tmp2VecMemNodeType = "<<divAuxData.tmp2VecMemNodeType.size()<<std::endl;
	// 	for (uint j = 0; j < divAuxData.tmp2VecMemNodeType.size();j++){
	// 		std::cout<<"divAuxData.tmp2VecMem["<<j<<"] = "<<divAuxData.tmp2VecMem[j].x<<" "<<divAuxData.tmp2VecMem[j].y<<" "<<divAuxData.tmp2VecMem[j].z<<", Type = "<<divAuxData.tmp2VecMemNodeType[j]<<std::endl;
	// 	}
	}

	//Here we will try to adjust the positions of newly added nodes (from ptsBetween1 and ptsBetween2) according to the 'center line' of the cell.

	assert(divAuxData.tmp1VecMem.size() <= membThreshold);
	assert(divAuxData.tmp2VecMem.size() <= membThreshold);
}

void SceCells::processMemVec_Ver2(uint i, std::vector<VecValT>& tmp1,
		std::vector<VecValT>& tmp2, CVector oldNucleusCenter,
		std::vector<CVector>& cellCenterLine_Basal2Apical,
		std::vector<CVector>& cellCenterLine_Apical2Basal,
		std::vector<CVector>& cellCenterLine_Basal2Apical_leftShift,
		std::vector<CVector>& cellCenterLine_Basal2Apical_rightShift,
		std::vector<CVector>& cellCenterLine_Apical2Basal_leftShift,
		std::vector<CVector>& cellCenterLine_Apical2Basal_rightShift,
		std::vector<double>& cellCenterLine_MirrorLength_Basal2Apical,
		std::vector<double>& cellCenterLine_MirrorLength_Apical2Basal) {
	divAuxData.tmp1VecMem.clear();
	divAuxData.tmp2VecMem.clear();
	divAuxData.tmp1VecMemNodeType.clear(); //Ali
	divAuxData.tmp2VecMemNodeType.clear(); //Ali

	uint membThreshold = allocPara_m.maxMembrNodePerCell;

	// std::sort(tmp1.begin(), tmp1.end());
	// std::sort(tmp2.begin(), tmp2.end());
	// Here we perform a cross-product computation to check the node order orientation. If it is clockwise, reverse it so it is counter-clockwise.
	double vec1x = (tmp1[0].vec.x - oldNucleusCenter.x);
	double vec2x = (tmp1[1].vec.x - oldNucleusCenter.x);
	double vec1y = (tmp1[0].vec.y - oldNucleusCenter.y);
	double vec2y = (tmp1[1].vec.y - oldNucleusCenter.y);
	double orientation = vec1x*vec2y - vec2x*vec1y;
	std::cout<<"tmp1 orientation = "<<orientation<<std::endl;
	if (orientation < 0){
		std::cout<<"tmp1 orientation incorrect. Rearranging tmp1."<<std::endl;
		std::vector<VecValT> tmp_tmp1 = tmp1;
		for (int i = 0; i < tmp1.size(); i++){
			// std::cout<<"tmp1 "<<tmp1[i].vec.x<<" "<<tmp1[i].vec.y<<" "<<tmp1[i].val<<" "<<tmp1[i].type<<std::endl;
		}
		tmp1.clear();
		for (int i = 0; i < tmp_tmp1.size(); i++){
			tmp1.push_back(tmp_tmp1[(tmp_tmp1.size()-1) - i]);
		}
		// for (int i = 0; i < tmp_tmp1.size(); i++){
		// 	std::cout<<"tmp_tmp1 "<<tmp_tmp1[i].vec.x<<" "<<tmp_tmp1[i].vec.y<<" "<<tmp_tmp1[i].val<<" "<<tmp_tmp1[i].type<<std::endl;
		// }
	}
	vec1x = (tmp2[0].vec.x - oldNucleusCenter.x);
	vec2x = (tmp2[1].vec.x - oldNucleusCenter.x);
	vec1y = (tmp2[0].vec.y - oldNucleusCenter.y);
	vec2y = (tmp2[1].vec.y - oldNucleusCenter.y);
	orientation = vec1x*vec2y - vec2x*vec1y;
	std::cout<<"tmp2 orientation = "<<orientation<<std::endl;
	if (orientation < 0){
		std::cout<<"tmp2 orientation incorrect. Rearranging tmp2."<<std::endl;
		std::vector<VecValT> tmp_tmp2 = tmp2;
		for (int i = 0; i < tmp2.size(); i++){
			// std::cout<<"tmp2 "<<tmp2[i].vec.x<<" "<<tmp2[i].vec.y<<" "<<tmp2[i].val<<" "<<tmp2[i].type<<std::endl;
		}
		tmp2.clear();
		for (int i = 0; i < tmp_tmp2.size(); i++){
			tmp2.push_back(tmp_tmp2[(tmp_tmp2.size()-1) - i]);
		}
		// for (int i = 0; i < tmp_tmp2.size(); i++){
		// 	std::cout<<"tmp_tmp2 "<<tmp_tmp2[i].vec.x<<" "<<tmp_tmp2[i].vec.y<<" "<<tmp_tmp2[i].val<<" "<<tmp_tmp2[i].type<<std::endl;
		// }
	}

	// std::vector<VecValT> sorted_tmp1 = tmp1;
	// std::vector<VecValT> sorted_tmp2 = tmp2;
	// std::sort(sorted_tmp1.begin(), sorted_tmp1.end());
	// std::sort(sorted_tmp2.begin(), sorted_tmp2.end());
	// for (int j = 0; j < sorted_tmp1.size(); j++){
	// 	std::cout<<"sorted_tmp1["<<j<<"].val = "<<sorted_tmp1[j].val<<", vec = "<<sorted_tmp1[j].vec.x<<" "<<sorted_tmp1[j].vec.y<<std::endl;
	// }
	// for (int j = 0; j < sorted_tmp2.size(); j++){
	// 	std::cout<<"sorted_tmp2["<<j<<"].val = "<<sorted_tmp2[j].val<<", vec = "<<sorted_tmp2[j].vec.x<<" "<<sorted_tmp2[j].vec.y<<std::endl;
	// }

	// for (int j = 0; j < tmp1.size(); j++){
	// 	std::cout<<"tmp1["<<j<<"].val = "<<tmp1[j].val<<", vec = "<<tmp1[j].vec.x<<" "<<tmp1[j].vec.y<<" , type = "<<tmp2[j].type<<std::endl;
	// }
	// for (int j = 0; j < tmp2.size(); j++){
	// 	std::cout<<"tmp2["<<j<<"].val = "<<tmp2[j].val<<", vec = "<<tmp2[j].vec.x<<" "<<tmp2[j].vec.y<<" , type = "<<tmp2[j].type<<std::endl;
	// }
	// //assert(tmp1.size() < allocPara_m.maxMembrNodePerCell);
	// //assert(tmp2.size() < allocPara_m.maxMembrNodePerCell);

	uint maxDivMembrNodeCount1 = allocPara_m.maxMembrNodePerCell - tmp1.size();
	uint maxDivMembrNodeCount2 = allocPara_m.maxMembrNodePerCell - tmp2.size();

	std::vector<CVector> ptsBetween1, ptsBetween2;

	uint Num_of_NodeType0_MotherCell = 0, Num_of_NodeType0_DaughterCell = 0, Num_of_NodeType1_DaughterCell = 0, Num_of_NodeType1_MotherCell = 0;
	if (divAuxData.isMotherCellBehind[i] == true){
		// Recall that tmp1 will always be the mother cell data //Kevin
		for (int j = 0; j < tmp1.size(); j++){
			if (tmp1[j].type == lateralB){
				Num_of_NodeType0_MotherCell+=1;
			}
		}
		if (Num_of_NodeType0_MotherCell > maxDivMembrNodeCount1){
			std::cout<<"Too many new nodes are needed to be introduced for mother cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		for (int j = 0; j < tmp2.size(); j++){
			if (tmp2[j].type == lateralA){
				Num_of_NodeType1_DaughterCell+=1;
			}
		}
		if (Num_of_NodeType1_DaughterCell > maxDivMembrNodeCount2){
			std::cout<<"Too many new nodes are needed to be introducedj for daughter cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		// Num_of_NodeType0_MotherCell = ;
		// Num_of_NodeType1_DaughterCell = ;
	}
	else{
		for (int j = 0; j < tmp1.size(); j++){
			if (tmp1[j].type == lateralA){
				Num_of_NodeType1_MotherCell+=1;
			}
		}
		if (Num_of_NodeType1_MotherCell > maxDivMembrNodeCount1){
			std::cout<<"Too many new nodes are needed to be introduced for mother cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		for (int j = 0; j < tmp2.size(); j++){
			if (tmp2[j].type == lateralB){
				Num_of_NodeType0_DaughterCell+=1;
			}
		}
		if (Num_of_NodeType0_DaughterCell > maxDivMembrNodeCount2){
			std::cout<<"Too many new nodes are needed to be introduced for daughter cell! Change the max number of mem nodes allowed!"<<std::endl;
		}
		
	}
	std::cout<<"Num_of_NodeType0_MotherCell = "<<Num_of_NodeType0_MotherCell<<std::endl;
	std::cout<<"Num_of_NodeType1_MotherCell = "<<Num_of_NodeType1_MotherCell<<std::endl;
	std::cout<<"Num_of_NodeType0_DaughterCell = "<<Num_of_NodeType0_DaughterCell<<std::endl;
	std::cout<<"Num_of_NodeType1_DaughterCell = "<<Num_of_NodeType1_DaughterCell<<std::endl;

	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp1.size() >= 1) {
		if (divAuxData.isMotherCellBehind[i] == true){
			// ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
			// 	Num_of_NodeType0_MotherCell, maxDivMembrNodeCount1);
			ptsBetween1 = obtainPtsBetween_cellCenterLine(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
				cellCenterLine_Basal2Apical,
				cellCenterLine_Basal2Apical_leftShift,
				cellCenterLine_MirrorLength_Basal2Apical);
		}
		else{
			// ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
			// 	Num_of_NodeType1_MotherCell, maxDivMembrNodeCount1);
			ptsBetween1 = obtainPtsBetween_cellCenterLine(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
				cellCenterLine_Apical2Basal,
				cellCenterLine_Apical2Basal_rightShift,
				cellCenterLine_MirrorLength_Apical2Basal);
		}
		// std::cout<<"tmp1[tmp1.size() - 1] = "<<tmp1[tmp1.size() - 1].vec.x<<" "<<tmp1[tmp1.size()-1].vec.y<<std::endl;
		// std::cout<<"tmp1[0] = "<<tmp1[0].vec.x<<" "<<tmp1[0].vec.y<<std::endl;
		// ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount1);
		// std::cout<<"sorted_tmp1[tmp1.size() - 1].val = "<<sorted_tmp1[sorted_tmp1.size()-1].val<<", .vec = "<< sorted_tmp1[sorted_tmp1.size() - 1].vec.x<<" "<<sorted_tmp1[sorted_tmp1.size()-1].vec.y<<std::endl;
		// std::cout<<"sorted_tmp1[0].val = "<<sorted_tmp1[0].val<<", .vec = "<<sorted_tmp1[0].vec.x<<" "<<sorted_tmp1[0].vec.y<<std::endl;
		// ptsBetween1 = obtainPtsBetween(sorted_tmp1[sorted_tmp1.size() - 1].vec, sorted_tmp1[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount1);
	}
	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp2.size() >= 1) {
		if (divAuxData.isMotherCellBehind[i] == true){
			// ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
			// 	Num_of_NodeType1_DaughterCell, maxDivMembrNodeCount2);
			ptsBetween2 = obtainPtsBetween_cellCenterLine(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
				cellCenterLine_Apical2Basal,
				cellCenterLine_Apical2Basal_rightShift,
				cellCenterLine_MirrorLength_Apical2Basal);
		}
		else{
			// ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
			// 	Num_of_NodeType0_DaughterCell, maxDivMembrNodeCount2);
			ptsBetween2 = obtainPtsBetween_cellCenterLine(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
				cellCenterLine_Basal2Apical,
				cellCenterLine_Basal2Apical_leftShift,
				cellCenterLine_MirrorLength_Basal2Apical);
		}
		// std::cout<<"tmp2[tmp1.size() - 1] = "<<tmp2[tmp2.size() - 1].vec.x<<" "<<tmp2[tmp2.size()-1].vec.y<<std::endl;
		// std::cout<<"tmp2[0] = "<<tmp2[0].vec.x<<" "<<tmp2[0].vec.y<<std::endl;
		// ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount2);
		// std::cout<<"sorted_tmp2[sorted_tmp2.size() - 1].val = "<<sorted_tmp2[sorted_tmp2.size()-1].val<<", .vec = "<<sorted_tmp2[sorted_tmp2.size() - 1].vec.x<<" "<<sorted_tmp2[sorted_tmp2.size()-1].vec.y<<std::endl;
		// std::cout<<"sorted_tmp2[0].val = "<<sorted_tmp2[0].val<<", vec = "<<sorted_tmp2[0].vec.x<<" "<<sorted_tmp2[0].vec.y<<std::endl;
		// ptsBetween2 = obtainPtsBetween(sorted_tmp2[sorted_tmp2.size() - 1].vec, sorted_tmp2[0].vec,
		// 		memNewSpacing, maxDivMembrNodeCount2);
	}

	for (uint j = 0; j < tmp1.size(); j++) {
		divAuxData.tmp1VecMem.push_back(tmp1[j].vec);
		divAuxData.tmp1VecMemNodeType.push_back(tmp1[j].type);
	}
	for (uint j = 0; j < tmp2.size(); j++) {
		divAuxData.tmp2VecMem.push_back(tmp2[j].vec);
		divAuxData.tmp2VecMemNodeType.push_back(tmp2[j].type);
	}
	std::cout<<"ptsBetween1 size = "<<ptsBetween1.size()<<std::endl;
	if (divAuxData.isMotherCellBehind[i] == false){
		for (uint j = 0; j < ptsBetween1.size(); j++) {
			divAuxData.tmp1VecMem.push_back(ptsBetween1[j]);
			// std::cout<<"ptsBtween1 "<<ptsBetween1[j].x<<" "<<ptsBetween1[j].y<<std::endl;
			divAuxData.tmp1VecMemNodeType.push_back(lateralB); 
		}
		// std::cout<<"size of tmp1VecMemNodeType = "<<divAuxData.tmp1VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp1VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp1VecMem["<<j<<"] = "<<divAuxData.tmp1VecMem[j].x<<" "<<divAuxData.tmp1VecMem[j].y<<" "<<divAuxData.tmp1VecMem[j].z<<", Type = "<<divAuxData.tmp1VecMemNodeType[j]<<std::endl;
		// }
	}
	else if (divAuxData.isMotherCellBehind[i] == true){
		for (uint j = 0; j < ptsBetween1.size(); j++) {
			divAuxData.tmp1VecMem.push_back(ptsBetween1[j]);
			// std::cout<<"ptsBtween1 "<<ptsBetween1[j].x<<" "<<ptsBetween1[j].y<<std::endl;
			divAuxData.tmp1VecMemNodeType.push_back(lateralA); 
		}
		std::cout<<"size of tmp1VecMemNodeType = "<<divAuxData.tmp1VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp1VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp1VecMem["<<j<<"] = "<<divAuxData.tmp1VecMem[j].x<<" "<<divAuxData.tmp1VecMem[j].y<<" "<<divAuxData.tmp1VecMem[j].z<<", Type = "<<divAuxData.tmp1VecMemNodeType[j]<<std::endl;
		// }
	}
	std::cout<<"ptsBetween2 size = "<<ptsBetween2.size()<<std::endl;
	if (divAuxData.isMotherCellBehind[i] == false){
		for (uint j = 0; j < ptsBetween2.size(); j++) {
			divAuxData.tmp2VecMem.push_back(ptsBetween2[j]);
			// std::cout<<"ptsBtween2 "<<ptsBetween2[j].x<<" "<<ptsBetween2[j].y<<std::endl;
			divAuxData.tmp2VecMemNodeType.push_back(lateralA);
		}
		// std::cout<<"size of tmp2VecMemNodeType = "<<divAuxData.tmp2VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp2VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp2VecMem["<<j<<"] = "<<divAuxData.tmp2VecMem[j].x<<" "<<divAuxData.tmp2VecMem[j].y<<" "<<divAuxData.tmp2VecMem[j].z<<", Type = "<<divAuxData.tmp2VecMemNodeType[j]<<std::endl;
		// }
	}
	else if (divAuxData.isMotherCellBehind[i] == true){
		for (uint j = 0; j < ptsBetween2.size(); j++) {
			divAuxData.tmp2VecMem.push_back(ptsBetween2[j]);
			// std::cout<<"ptsBtween2 "<<ptsBetween2[j].x<<" "<<ptsBetween2[j].y<<std::endl;
			divAuxData.tmp2VecMemNodeType.push_back(lateralB);
		}
		std::cout<<"size of tmp2VecMemNodeType = "<<divAuxData.tmp2VecMemNodeType.size()<<std::endl;
		// for (uint j = 0; j < divAuxData.tmp2VecMemNodeType.size();j++){
		// 	std::cout<<"divAuxData.tmp2VecMem["<<j<<"] = "<<divAuxData.tmp2VecMem[j].x<<" "<<divAuxData.tmp2VecMem[j].y<<" "<<divAuxData.tmp2VecMem[j].z<<", Type = "<<divAuxData.tmp2VecMemNodeType[j]<<std::endl;
		// }
	}

	//Here we will try to adjust the positions of newly added nodes (from ptsBetween1 and ptsBetween2) according to the 'center line' of the cell.

	assert(divAuxData.tmp1VecMem.size() <= membThreshold);
	assert(divAuxData.tmp2VecMem.size() <= membThreshold);
}

void SceCells::obtainMembrAndIntnlNodes(uint i, vector<CVector>& membrNodes,
		vector<CVector>& intnlNodes) {
	membrNodes.clear();
	intnlNodes.clear();

	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (divAuxData.tmpIsActive_M[index] != true) {
			continue;
		}
		double posX = divAuxData.tmpNodePosX_M[index];
		double posY = divAuxData.tmpNodePosY_M[index];
		if (j < membThreshold) {
			// means node type is membrane
			CVector memPos(posX, posY, 0);
			membrNodes.push_back(memPos);
		} else {
			CVector intnlPos(posX, posY, 0);
			intnlNodes.push_back(intnlPos);
		}
	}
}
//Ali

void SceCells::obtainMembrAndIntnlNodesPlusNodeType(uint i, vector<CVector>& membrNodes,
		vector<CVector>& intnlNodes, vector<MembraneType1> & nodeTypeIndxDiv) {
	membrNodes.clear();
	intnlNodes.clear();
	nodeTypeIndxDiv.clear() ; 

	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (divAuxData.tmpIsActive_M[index] != true) {
			continue;
		}
		double posX = divAuxData.tmpNodePosX_M[index];
		double posY = divAuxData.tmpNodePosY_M[index];
		MembraneType1    nodeTypeI = divAuxData.tmpNodeType[index];
		if (j < membThreshold) {
			// means node type is membrane
			CVector memPos(posX, posY, 0);
			membrNodes.push_back(memPos);
			nodeTypeIndxDiv.push_back(nodeTypeI) ; 
		} else {
			CVector intnlPos(posX, posY, 0);
			intnlNodes.push_back(intnlPos);
		}
	}
}

void SceCells::obtainMembrAndIntnlNodesPlusNodeType2(uint i, vector<CVector>& membrNodes,
		vector<CVector>& intnlNodes, vector<MembraneType1> & nodeTypeIndxDiv, 
		vector<CVector>& cellCenterLine_Basal2Apical, vector<CVector>& cellCenterLine_Apical2Basal, 
		vector<CVector>& cellCenterLine_Basal2Apical_leftShift, vector<CVector>& cellCenterLine_Basal2Apical_rightShift,
		vector<CVector>& cellCenterLine_Apical2Basal_leftShift, vector<CVector>& cellCenterLine_Apical2Basal_rightShift,
		vector<double>& cellCenterLine_MirrorLength_Basal2Apical, vector<double>& cellCenterLine_MirrorLength_Apical2Basal) {
	membrNodes.clear();
	intnlNodes.clear();
	nodeTypeIndxDiv.clear() ; 
	cellCenterLine_Basal2Apical.clear();
	cellCenterLine_Apical2Basal.clear();
	cellCenterLine_Apical2Basal_leftShift.clear();
	cellCenterLine_Apical2Basal_rightShift.clear();
	cellCenterLine_Basal2Apical_leftShift.clear();
	cellCenterLine_Basal2Apical_rightShift.clear();
	cellCenterLine_MirrorLength_Basal2Apical.clear();
	cellCenterLine_MirrorLength_Apical2Basal.clear();

	vector<CVector> tmpCellCenterLine_Basal2Apical;
	vector<CVector> tmpCellCenterLine_Basal2Apical_leftShift;
	vector<CVector> tmpCellCenterLine_Basal2Apical_rightShift;
	vector<double> tmpCellCenterLine_MirrorLength_Basal2Apical;
	vector<CVector> tmpCellCenterLine_Apical2Basal;
	vector<CVector> tmpCellCenterLine_Apical2Basal_leftShift;
	vector<CVector> tmpCellCenterLine_Apical2Basal_rightShift;
	vector<double> tmpCellCenterLine_MirrorLength_Apical2Basal;

	tmpCellCenterLine_Basal2Apical.clear();
	tmpCellCenterLine_Basal2Apical_leftShift.clear();
	tmpCellCenterLine_Basal2Apical_rightShift.clear();
	tmpCellCenterLine_MirrorLength_Basal2Apical.clear();
	tmpCellCenterLine_Apical2Basal.clear();
	tmpCellCenterLine_Apical2Basal_leftShift.clear();
	tmpCellCenterLine_Apical2Basal_rightShift.clear();
	tmpCellCenterLine_MirrorLength_Apical2Basal.clear();

	bool firstloop_Basal2Apical = true;
	bool firstloop_Apical2Basal = true;
	bool lateralA_earlyShift = false;
	bool lateralB_earlyShift = false;

	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index, nextIndex;
	int lateralACount = 0;
	int lateralBCount = 0;
	int initLateralACount = 0;
	int initLateralBCount = 0;
	std::cout<<"i = "<<i<<std::endl;

	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (divAuxData.tmpNodeType[index]==lateralA){
			initLateralACount+=1;
		}
		if (divAuxData.tmpNodeType[index]==lateralB){
			initLateralBCount+=1;
		}
	}
	std::cout<<"initLateralACount = "<<initLateralACount<<std::endl;
	std::cout<<"initLateralBCount = "<<initLateralBCount<<std::endl;

	for (uint j = 0; j < maxAllNodePerCell; j++) {
		// std::cout<<"j = "<<j<<std::endl;
		index = i * maxAllNodePerCell + j;
		nextIndex = index + 1;
		if (nextIndex >= i*maxAllNodePerCell+maxAllNodePerCell){
			nextIndex = i*maxAllNodePerCell + 0;
		}
		if (divAuxData.tmpIsActive_M[index] != true) {
			continue;
		}
		double posX = divAuxData.tmpNodePosX_M[index];
		// std::cout<<"posX = "<<posX<<std::endl;
		double posY = divAuxData.tmpNodePosY_M[index];
		// std::cout<<"posY = "<<posY<<std::endl;
		MembraneType1    nodeTypeI = divAuxData.tmpNodeType[index];
		// std::cout<<"nodetype = "<<nodeTypeI<<std::endl;
		// std::cout<<"nodeTypeI = "<<nodeTypeI<<std::endl;
		if (j < membThreshold) {
			// means node type is membrane
			CVector memPos(posX, posY, 0);
			membrNodes.push_back(memPos);
			nodeTypeIndxDiv.push_back(nodeTypeI) ; 
			CVector cellCenterLinePos;
			CVector cellCenterLinePos_leftUnitDir;
			CVector cellCenterLinePos_rightUnitDir;
			//Since we start by using nodes with label 'lateralA', the resulting cellCenterLinePos shoould orient from basal to apical.
			//But since mother and daughter cell will go through this center line in different orientation (even though both cells use
			//counterclockwise orientation). //Kevin
			if (nodeTypeI == lateralA){
				// int mirrorNode = nodes->getInfoVecs().nodeAdhereIndex[j];
				// std::cout<<"divAuxData.tmpNodeMemMirrorIndex_M size = "<<divAuxData.tmpNodeMemMirrorIndex_M.size()<<std::endl;
				cellCenterLinePos.x = (divAuxData.tmpNodePosX_M[index] + nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				cellCenterLinePos.y = (divAuxData.tmpNodePosY_M[index] + nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				// std::cout<<"cellCenterLinePos = "<<cellCenterLinePos.x<<" "<<cellCenterLinePos.y<<std::endl;
				double length = sqrt((divAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(divAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]]) +
								(divAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(divAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]]));

				cellCenterLinePos_rightUnitDir.x = divAuxData.tmpNodePosX_M[index] - cellCenterLinePos.x;
				// std::cout<<"tmpNodePosX = "<<divAuxData.tmpNodePosX_M[index]<<" , cellCenterLinePos.x = "<<cellCenterLinePos.x<<std::endl;
				cellCenterLinePos_rightUnitDir.y = divAuxData.tmpNodePosY_M[index] - cellCenterLinePos.y;
				// std::cout<<"tmpNodePosY = "<<divAuxData.tmpNodePosY_M[index]<<" , cellCenterLinePos.y = "<<cellCenterLinePos.y<<std::endl;
				double rightLength = sqrt(cellCenterLinePos_rightUnitDir.x*cellCenterLinePos_rightUnitDir.x + cellCenterLinePos_rightUnitDir.y*cellCenterLinePos_rightUnitDir.y);
				// std::cout<<"rightLength = "<<rightLength<<std::endl;
				cellCenterLinePos_rightUnitDir.x = cellCenterLinePos_rightUnitDir.x/rightLength;
				cellCenterLinePos_rightUnitDir.y = cellCenterLinePos_rightUnitDir.y/rightLength;

				cellCenterLinePos_leftUnitDir.x = nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]] - cellCenterLinePos.x;
				// std::cout<<"nodeLocX[tmpNodeMirrorIndex] = "<<nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]]<<" , cellCenterLinePos.x = "<<cellCenterLinePos.x<<std::endl;
				cellCenterLinePos_leftUnitDir.y = nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]] - cellCenterLinePos.y;
				// std::cout<<"nodeLocY[tmpNodeMirrorIndex] = "<<nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]]<<" , cellCenterLinePos.y = "<<cellCenterLinePos.y<<std::endl;
				double leftLength = sqrt(cellCenterLinePos_leftUnitDir.x*cellCenterLinePos_leftUnitDir.x + cellCenterLinePos_leftUnitDir.y*cellCenterLinePos_leftUnitDir.y);
				// std::cout<<"leftLength = "<<leftLength<<std::endl;
				cellCenterLinePos_leftUnitDir.x = cellCenterLinePos_leftUnitDir.x/leftLength;
				cellCenterLinePos_leftUnitDir.y = cellCenterLinePos_leftUnitDir.y/leftLength;

				cellCenterLine_Basal2Apical.push_back(cellCenterLinePos);
				lateralACount+=1;
				cellCenterLine_Basal2Apical_leftShift.push_back(cellCenterLinePos_leftUnitDir);
				// std::cout<<"cellCenterLine Basal2Apical left shift: "<<cellCenterLinePos_leftUnitDir.x<<" "<<cellCenterLinePos_leftUnitDir.y<<std::endl;
				cellCenterLine_Basal2Apical_rightShift.push_back(cellCenterLinePos_rightUnitDir);
				// std::cout<<"cellCenterLine Basal2Apical right shift: "<<cellCenterLinePos_rightUnitDir.x<<" "<<cellCenterLinePos_rightUnitDir.y<<std::endl;
				cellCenterLine_MirrorLength_Basal2Apical.push_back(length);

				if (firstloop_Basal2Apical == true){
					tmpCellCenterLine_Basal2Apical.push_back(cellCenterLinePos);
					tmpCellCenterLine_Basal2Apical_leftShift.push_back(cellCenterLinePos_leftUnitDir);
					tmpCellCenterLine_Basal2Apical_rightShift.push_back(cellCenterLinePos_rightUnitDir);
					tmpCellCenterLine_MirrorLength_Basal2Apical.push_back(length);
				}

				if (nodeTypeI == lateralA && divAuxData.tmpNodeType[nextIndex] != lateralA){
					// std::cout<<"Node type shift away from lateralA."<<std::endl;
					//This means that we hit a transition point from lateralA to apical1. Need to make sure if we have covered all lateralA nodes or not.
					//If not, this means that we started counting from the middle of lateralA nodes, which will causes serious problem later.
					//Hence we will try to remedy this issue here //Kevin
					if (lateralACount != initLateralACount){
						lateralA_earlyShift = true;
						// std::cout<<"Node type shift away from lateralA before covering all lateralA type node."<<std::endl;
						cellCenterLine_Basal2Apical.clear();
						cellCenterLine_Basal2Apical_leftShift.clear();
						cellCenterLine_Basal2Apical_rightShift.clear();
						cellCenterLine_MirrorLength_Basal2Apical.clear();
					}
					firstloop_Basal2Apical = false;
					// std::cout<<"firstloop_Basal2Apical = "<<firstloop_Basal2Apical<<std::endl;
				}
			}
			if (nodeTypeI == lateralB){
				// int mirrorNode = nodes->getInfoVecs().nodeAdhereIndex[j];
				cellCenterLinePos.x = (divAuxData.tmpNodePosX_M[index] + nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				cellCenterLinePos.y = (divAuxData.tmpNodePosY_M[index] + nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				double length = sqrt((divAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(divAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]]) +
								(divAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(divAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]]));

				cellCenterLinePos_leftUnitDir.x = divAuxData.tmpNodePosX_M[index] - cellCenterLinePos.x;
				cellCenterLinePos_leftUnitDir.y = divAuxData.tmpNodePosY_M[index] - cellCenterLinePos.y;
				double leftLength = sqrt(cellCenterLinePos_leftUnitDir.x*cellCenterLinePos_leftUnitDir.x + cellCenterLinePos_leftUnitDir.y*cellCenterLinePos_leftUnitDir.y);
				cellCenterLinePos_leftUnitDir.x = cellCenterLinePos_leftUnitDir.x/leftLength;
				cellCenterLinePos_leftUnitDir.y = cellCenterLinePos_leftUnitDir.y/leftLength;

				cellCenterLinePos_rightUnitDir.x = nodes->getInfoVecs().nodeLocX[divAuxData.tmpNodeMemMirrorIndex_M[index]] - cellCenterLinePos.x;
				cellCenterLinePos_rightUnitDir.y = nodes->getInfoVecs().nodeLocY[divAuxData.tmpNodeMemMirrorIndex_M[index]] - cellCenterLinePos.y;
				double rightLength = sqrt(cellCenterLinePos_rightUnitDir.x*cellCenterLinePos_rightUnitDir.x + cellCenterLinePos_rightUnitDir.y*cellCenterLinePos_rightUnitDir.y);
				cellCenterLinePos_rightUnitDir.x = cellCenterLinePos_rightUnitDir.x/rightLength;
				cellCenterLinePos_rightUnitDir.y = cellCenterLinePos_rightUnitDir.y/rightLength;

				cellCenterLine_Apical2Basal.push_back(cellCenterLinePos);
				lateralBCount+=1;
				cellCenterLine_Apical2Basal_leftShift.push_back(cellCenterLinePos_leftUnitDir);
				// std::cout<<"cellCenterLine Apical2Basal left shift: "<<cellCenterLinePos_leftUnitDir.x<<" "<<cellCenterLinePos_leftUnitDir.y<<std::endl;
				cellCenterLine_Apical2Basal_rightShift.push_back(cellCenterLinePos_rightUnitDir);
				// std::cout<<"cellCenterLine Apical2Basal right shift: "<<cellCenterLinePos_rightUnitDir.x<<" "<<cellCenterLinePos_rightUnitDir.y<<std::endl;
				cellCenterLine_MirrorLength_Apical2Basal.push_back(length);

				if (firstloop_Apical2Basal == true){
					tmpCellCenterLine_Apical2Basal.push_back(cellCenterLinePos);
					tmpCellCenterLine_Apical2Basal_leftShift.push_back(cellCenterLinePos_leftUnitDir);
					tmpCellCenterLine_Apical2Basal_rightShift.push_back(cellCenterLinePos_rightUnitDir);
					tmpCellCenterLine_MirrorLength_Apical2Basal.push_back(length);
				}

				if (nodeTypeI == lateralB && divAuxData.tmpNodeType[nextIndex] != lateralB){
					// std::cout<<"Node type shift away from lateralB."<<std::endl;
					//This means that we hit a transition point from lateralA to apical1. Need to make sure if we have covered all lateralA nodes or not.
					//If not, this means that we started counting from the middle of lateralA nodes, which will causes serious problem later.
					//Hence we will try to remedy this issue here //Kevin
					if (lateralBCount != initLateralBCount){
						lateralB_earlyShift = true;
						// std::cout<<"Node type shift away from lateralB before covering all lateralB type node."<<std::endl;
						cellCenterLine_Apical2Basal.clear();
						cellCenterLine_Apical2Basal_leftShift.clear();
						cellCenterLine_Apical2Basal_rightShift.clear();
						cellCenterLine_MirrorLength_Apical2Basal.clear();
					}
					firstloop_Apical2Basal = false;
					// std::cout<<"firstloop_Apical2Basal = "<<firstloop_Apical2Basal<<std::endl;
				}
			}
		} else {
			CVector intnlPos(posX, posY, 0);
			intnlNodes.push_back(intnlPos);
		}
	}

	if (lateralA_earlyShift == true){
		for (int p = 0; p < tmpCellCenterLine_Basal2Apical.size(); p++){
			cellCenterLine_Basal2Apical.push_back(tmpCellCenterLine_Basal2Apical[p]);
			cellCenterLine_Basal2Apical_leftShift.push_back(tmpCellCenterLine_Basal2Apical_leftShift[p]);
			cellCenterLine_Basal2Apical_rightShift.push_back(tmpCellCenterLine_Basal2Apical_rightShift[p]);
			cellCenterLine_MirrorLength_Basal2Apical.push_back(tmpCellCenterLine_MirrorLength_Basal2Apical[p]);
		}
	}
	
	// std::cout<<"cellCenterLine_Basal2Apical size = "<<cellCenterLine_Basal2Apical.size()<<std::endl;
	// std::cout<<"cellCenterLine_Basal2Apical_leftShift size = "<<cellCenterLine_Basal2Apical_leftShift.size()<<std::endl;
	// std::cout<<"cellCenterLine_Basal2Apical_rightShift size = "<<cellCenterLine_Basal2Apical_rightShift.size()<<std::endl;
	// std::cout<<"cellCenterLine_MirrorLength_Basal2Apical size = "<<cellCenterLine_MirrorLength_Basal2Apical.size()<<std::endl;

	if (lateralB_earlyShift == true){
		for (int p = 0; p < tmpCellCenterLine_Apical2Basal.size(); p++){
			cellCenterLine_Apical2Basal.push_back(tmpCellCenterLine_Apical2Basal[p]);
			cellCenterLine_Apical2Basal_leftShift.push_back(tmpCellCenterLine_Apical2Basal_leftShift[p]);
			cellCenterLine_Apical2Basal_rightShift.push_back(tmpCellCenterLine_Apical2Basal_rightShift[p]);
			cellCenterLine_MirrorLength_Apical2Basal.push_back(tmpCellCenterLine_MirrorLength_Apical2Basal[p]);
			}
	}

	// std::cout<<"cellCenterLine_Apical2Basal size = "<<cellCenterLine_Apical2Basal.size()<<std::endl;
	// std::cout<<"cellCenterLine_Apical2Basal_leftShift size = "<<cellCenterLine_Apical2Basal_leftShift.size()<<std::endl;
	// std::cout<<"cellCenterLine_Apical2Basal_rightShift size = "<<cellCenterLine_Apical2Basal_rightShift.size()<<std::endl;
	// std::cout<<"cellCenterLine_MirrorLength_Apical2Basal size = "<<cellCenterLine_MirrorLength_Apical2Basal.size()<<std::endl;
	//Now we check to see if the cellCenterLine data structure is built correctly. The two should be identical except with opposite ordering.
}

void SceCells::obtainMembrAndIntnlNodesPlusNodeType2_printingOnly(uint i, vector<CVector>& membrNodes,
		vector<CVector>& intnlNodes, vector<MembraneType1> & nodeTypeIndxDiv, 
		vector<CVector>& cellCenterLine_Basal2Apical, vector<CVector>& cellCenterLine_Apical2Basal) {
	membrNodes.clear();
	intnlNodes.clear();
	nodeTypeIndxDiv.clear() ; 
	cellCenterLine_Basal2Apical.clear();
	cellCenterLine_Apical2Basal.clear();

	vector<CVector> tmpCellCenterLine_Basal2Apical;
	
	vector<CVector> tmpCellCenterLine_Apical2Basal;
	
	tmpCellCenterLine_Basal2Apical.clear();
	tmpCellCenterLine_Apical2Basal.clear();
	
	bool firstloop_Basal2Apical = true;
	bool firstloop_Apical2Basal = true;
	bool lateralA_earlyShift = false;
	bool lateralB_earlyShift = false;

	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index, nextIndex;
	int lateralACount = 0;
	int lateralBCount = 0;
	int initLateralACount = 0;
	int initLateralBCount = 0;
	std::cout<<"i = "<<i<<std::endl;

	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (printAuxData.tmpNodeType[index]==lateralA){
			initLateralACount+=1;
		}
		if (printAuxData.tmpNodeType[index]==lateralB){
			initLateralBCount+=1;
		}
	}
	std::cout<<"initLateralACount_printOnly = "<<initLateralACount<<std::endl;
	std::cout<<"initLateralBCount_printOnly = "<<initLateralBCount<<std::endl;

	for (uint j = 0; j < maxAllNodePerCell; j++) {
		// std::cout<<"j = "<<j<<std::endl;
		index = i * maxAllNodePerCell + j;
		nextIndex = index + 1;
		if (nextIndex >= i*maxAllNodePerCell+maxAllNodePerCell){
			nextIndex = i*maxAllNodePerCell + 0;
		}
		if (printAuxData.tmpIsActive_M[index] != true) {
			continue;
		}
		// std::cout<<"debug 1"<<std::endl;
		double posX = printAuxData.tmpNodePosX_M[index];
		double posY = printAuxData.tmpNodePosY_M[index];
		MembraneType1    nodeTypeI = printAuxData.tmpNodeType[index];
		// std::cout<<"debug 2"<<std::endl;
		if (j < membThreshold) {
			CVector memPos(posX, posY, 0);
			membrNodes.push_back(memPos);
			nodeTypeIndxDiv.push_back(nodeTypeI) ; 
			CVector cellCenterLinePos;
			CVector cellCenterLinePos_leftUnitDir;
			CVector cellCenterLinePos_rightUnitDir;
			if (nodeTypeI == lateralA){
				cellCenterLinePos.x = (printAuxData.tmpNodePosX_M[index] + nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				cellCenterLinePos.y = (printAuxData.tmpNodePosY_M[index] + nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				double length = sqrt((printAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(printAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]]) +
								(printAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(printAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]]));
				cellCenterLine_Basal2Apical.push_back(cellCenterLinePos);
				lateralACount+=1;

				if (firstloop_Basal2Apical == true){
					tmpCellCenterLine_Basal2Apical.push_back(cellCenterLinePos);
				}

				if (nodeTypeI == lateralA && printAuxData.tmpNodeType[nextIndex] != lateralA){
					if (lateralACount != initLateralACount){
						lateralA_earlyShift = true;
						cellCenterLine_Basal2Apical.clear();
					}
					firstloop_Basal2Apical = false;
				}
			}
			if (nodeTypeI == lateralB){
				cellCenterLinePos.x = (printAuxData.tmpNodePosX_M[index] + nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				cellCenterLinePos.y = (printAuxData.tmpNodePosY_M[index] + nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]])/2.0;
				double length = sqrt((printAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(printAuxData.tmpNodePosX_M[index] - nodes->getInfoVecs().nodeLocX[printAuxData.tmpNodeMemMirrorIndex_M[index]]) +
								(printAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]])*
								(printAuxData.tmpNodePosY_M[index] - nodes->getInfoVecs().nodeLocY[printAuxData.tmpNodeMemMirrorIndex_M[index]]));
				cellCenterLine_Apical2Basal.push_back(cellCenterLinePos);
				lateralBCount+=1;

				if (firstloop_Apical2Basal == true){
					tmpCellCenterLine_Apical2Basal.push_back(cellCenterLinePos);
				}

				if (nodeTypeI == lateralB && printAuxData.tmpNodeType[nextIndex] != lateralB){
					if (lateralBCount != initLateralBCount){
						lateralB_earlyShift = true;
						cellCenterLine_Apical2Basal.clear();
					}
					firstloop_Apical2Basal = false;
				}
			}
		} else {
			CVector intnlPos(posX, posY, 0);
			intnlNodes.push_back(intnlPos);
		}
		// std::cout<<"debug 3"<<std::endl;
	}

	if (lateralA_earlyShift == true){
		for (int p = 0; p < tmpCellCenterLine_Basal2Apical.size(); p++){
			cellCenterLine_Basal2Apical.push_back(tmpCellCenterLine_Basal2Apical[p]);
		}
	}

	if (lateralB_earlyShift == true){
		for (int p = 0; p < tmpCellCenterLine_Apical2Basal.size(); p++){
			cellCenterLine_Apical2Basal.push_back(tmpCellCenterLine_Apical2Basal[p]);
			}
	}

}


/* Ali 
CVector SceCells::obtainCenter(uint i) {
	double oldCenterX = divAuxData.tmpCenterPosX_M[i];
	double oldCenterY = divAuxData.tmpCenterPosY_M[i];
	CVector centerPos(oldCenterX, oldCenterY, 0);
	return centerPos;
}
*/
CVector SceCells::obtainCellCenter(uint i) {
	double oldCenterX = divAuxData.tmpCenterPosX_M[i];
	double oldCenterY = divAuxData.tmpCenterPosY_M[i];
	CVector centerPos(oldCenterX, oldCenterY, 0);
	return centerPos;
}

CVector SceCells::obtainNucleusCenter(uint i, vector<CVector> IntnlNodes){
	for (int j = 0; j < IntnlNodes.size(); j++){
		divAuxData.tmpNucleusCenterPosX_M[i] += IntnlNodes[j].x;
		divAuxData.tmpNucleusCenterPosY_M[i] += IntnlNodes[j].y;
	}
	double oldCenterX = divAuxData.tmpNucleusCenterPosX_M[i]/IntnlNodes.size();
	double oldCenterY = divAuxData.tmpNucleusCenterPosY_M[i]/IntnlNodes.size();
	CVector centerPos(oldCenterX, oldCenterY, 0);
	return centerPos;
} //Kevin 

CVector SceCells::obtainIntCenter(uint i) {
	double oldCenterX = divAuxData.tmpCenterPosX_M[i];
	double oldCenterY = divAuxData.tmpCenterPosY_M[i];
	CVector centerPos(oldCenterX, oldCenterY, 0);
	return centerPos;
}

/*
CVector SceCells::calDivDir_MajorAxis(CVector center,
		vector<CVector>& membrNodes, double& lenAlongMajorAxis) {
// not the optimal algorithm but easy to code
	double maxDiff = 0;
	CVector majorAxisDir;
	for (uint i = 0; i < membrNodes.size(); i++) {
		CVector tmpDir = membrNodes[i] - center;
		CVector tmpUnitDir = tmpDir.getUnitVector();
		double min = 0, max = 0;
		for (uint j = 0; j < membrNodes.size(); j++) {
			CVector tmpDir2 = membrNodes[j] - center;
			double tmpVecProduct = tmpDir2 * tmpUnitDir;
			if (tmpVecProduct < min) {
				min = tmpVecProduct;
			}
			if (tmpVecProduct > max) {
				max = tmpVecProduct;
			}
		}
		double diff = max - min;
		if (diff > maxDiff) {
			maxDiff = diff;
			majorAxisDir = tmpUnitDir;
		}
	}
	lenAlongMajorAxis = maxDiff;
	return majorAxisDir;
}
*/

CVector SceCells::calDivDir_MajorAxis(CVector center,
		vector<CVector>& membrNodes, double& lenAlongMajorAxis) {
// not the optimal algorithm but easy to code
	double minDiff = 10000;
	CVector minorAxisDir;
	for (uint i = 0; i < membrNodes.size(); i++) {
		CVector tmpDir = membrNodes[i] - center;
		CVector tmpUnitDir = tmpDir.getUnitVector();
		double min = 0, max = 0;
		for (uint j = 0; j < membrNodes.size(); j++) {
			CVector tmpDir2 = membrNodes[j] - center;
			double tmpVecProduct = tmpDir2 * tmpUnitDir;
			if (tmpVecProduct < min) {
				min = tmpVecProduct;
			}
			if (tmpVecProduct > max) {
				max = tmpVecProduct;
			}
		}
		double diff = max - min;
		if (diff < minDiff) {
			minDiff = diff;
			minorAxisDir = tmpUnitDir;
		}
	}
	lenAlongMajorAxis = minDiff;
	return minorAxisDir;
}

CVector SceCells::calDivDir_ApicalBasal(CVector center,
		vector<CVector>& membrNodes, double& lenAlongMajorAxis, vector<MembraneType1> & nodeTypeIndxDiv) {
// not the optimal algorithm but easy to code
	double minDiff = 10000;
	CVector minorAxisDir;
	int minPointAdhIndex ; 
	int maxPointAdhIndex; 

	//for (uint i = 0; i < membrNodes.size(); i++) {
	//	cout <<"adhesion index for dividing cell node"<<i<<"is" << adhIndxDiv[i] <<endl; 
//	}	
	//return 0 ; 
	for (uint i = 0; i < membrNodes.size(); i++) {
		if ( (nodeTypeIndxDiv[i]!=lateralA) &&  (nodeTypeIndxDiv[i]!=lateralB) ) {
			continue ; 
		} 
		CVector tmpDir = membrNodes[i] - center;
		CVector tmpUnitDir = tmpDir.getUnitVector();
		double min = 0, max = 0;
		//distance finder for node i to the opposite nodes //Ali 
		for (uint j = 0; j < membrNodes.size(); j++) {
			CVector tmpDir2 = membrNodes[j] - center;
			double tmpVecProduct = tmpDir2 * tmpUnitDir;
			if (tmpVecProduct < min) {
				min = tmpVecProduct;
			}
			if (tmpVecProduct > max) {
				max = tmpVecProduct;
			}
		}
		double diff = max - min;
		// minimum distance finder for each cells to be used for cell center shifting. It should also need to be a node that have neighbor
		if (diff < minDiff ) {
			minDiff = diff;
			minorAxisDir = tmpUnitDir;
//			adhesionIndexFinal=adhIndxDiv[i]; 
		}
	}

	lenAlongMajorAxis = minDiff;
	return minorAxisDir;
}

std::pair <int ,int> SceCells::calApicalBasalRingIds(CVector divDir, CVector center,vector<CVector>& membrNodes, vector<MembraneType1> & nodeTypeIndxDiv) {

	int idMin, idMax ; 
	CVector splitDir = divDir.rotateNintyDeg_XY_CC();

	double min = 0, max = 0;
	for (uint j = 0; j < membrNodes.size(); j++) {
		CVector tmpDir2 = membrNodes[j] - center;
		CVector tmpUnitDir2 = tmpDir2.getUnitVector();
		double tmpVecProduct = splitDir * tmpUnitDir2;
		if (tmpVecProduct < min) {
			min = tmpVecProduct;
			idMin=j ; 
		}
		if (tmpVecProduct > max) {
			max = tmpVecProduct;
			idMax=j ; 
		}
	}
	cout << " contractile node location is " << membrNodes[idMin].x << " ," << membrNodes[idMin].y << endl ; 
	cout << " contractile node location is " << membrNodes[idMax].x << " ," << membrNodes[idMax].y << endl ; 

	if (nodeTypeIndxDiv[idMin]==apical1) {
		return make_pair(idMin,idMax) ;
	} else {

		return make_pair(idMax,idMin) ;
	}
}


//A&A
double SceCells::calLengthAlongHertwigAxis(CVector divDir, CVector cellCenter,
		vector<CVector>& membrNodes) {

	CVector divDirUnit = divDir.getUnitVector();


	double minUnit = 0, maxUnit = 0;
	double minOveral = 0, maxOveral = 0;
	for (uint i = 0; i < membrNodes.size(); i++) {
		CVector tmpDir = membrNodes[i] - cellCenter; //it is cell center
		CVector tmpUnitDir = tmpDir.getUnitVector();
			double tmpVecProductUnit = divDirUnit * tmpUnitDir;
			double tmpVecProductOveral = divDirUnit * tmpDir;
			if (tmpVecProductUnit < minUnit) {
				minUnit = tmpVecProductUnit;
				minOveral = tmpVecProductOveral;
			}
			if (tmpVecProductUnit > maxUnit) {
				maxUnit = tmpVecProductUnit;
				maxOveral = tmpVecProductOveral;
			}
	}
	
		double lenAlongHertwigAxis = maxOveral - minOveral;
	return lenAlongHertwigAxis; // it is minor axis
}


void SceCells::obtainTwoNewIntCenters(CVector& oldIntCenter, CVector& divDir,
		double len_MajorAxis, CVector& intCenterNew1, CVector& intCenterNew2) {

	CVector divDirUnit = divDir.getUnitVector();
	std::cout<<"divDirUnit = "<<divDirUnit.x<<" "<<divDirUnit.y<<std::endl;
	double lenChange = len_MajorAxis / 2.0 * centerShiftRatio; // this means small axis
	intCenterNew1 = oldIntCenter + lenChange * divDirUnit;   // it should be nucleus center
	intCenterNew2 = oldIntCenter - lenChange * divDirUnit;   // it should be nulceus center
	CVector centerTissue ;  //Ali
	centerTissue=CVector (40.0, 40.0, 0.0) ; //Ali should be imported
	CVector tmpVec1=intCenterNew1-centerTissue ;  //Ali // assuming New1 is mother cell
	CVector tmpVec2=intCenterNew2-centerTissue ;  //Ali 
	CVector tmpDiff=tmpVec2-tmpVec1 ; 
	CVector tmpCross=Cross(tmpVec1,tmpVec2) ; //Ali 


	bool isMotherCellBehindInt=false ;  //Ali 
	// assuming CCW is the initial order of cell ranks
	//if (tmpCross.z>0){
	// if (tmpDiff.x>0){
		isMotherCellBehindInt=true  ; 
	// }
//Ali
	std::cout<<"isMotherCellBehindInt = "<<isMotherCellBehindInt<<std::endl;
   divAuxData.isMotherCellBehind.push_back(isMotherCellBehindInt) ; 
}

void SceCells::prepareTmpVec(uint i, CVector divDir, CVector oldCellCenter,CVector oldIntCenter
		,std::vector<VecValT>& tmp1, std::vector<VecValT>& tmp2) {
	tmp1.clear(); // is for membrane node of first cell
	tmp2.clear(); // is for membrane node of the second cell
	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	VecValT tmpData;
	CVector splitDir = divDir.rotateNintyDeg_XY_CC();
	std::cout<<"splitDir = "<<splitDir.x<<" "<<splitDir.y<<std::endl;
	std::vector<double> oldCellCenter2BasalVec(2,0.0);
	std::vector<double> oldCellCenter2ApicalVec(2,0.0);
	oldCellCenter2BasalVec[0] = divAuxData.tmpBasalLoc[0] - oldCellCenter.x;
	oldCellCenter2BasalVec[1] = divAuxData.tmpBasalLoc[1] - oldCellCenter.y;
	oldCellCenter2ApicalVec[0] = divAuxData.tmpApicalLoc[0] - oldCellCenter.x;
	oldCellCenter2ApicalVec[1] = divAuxData.tmpApicalLoc[1] - oldCellCenter.y;
	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < membThreshold) {
			// means node type is membrane
			if (divAuxData.tmpIsActive_M[index] == true) {
				CVector memPos(divAuxData.tmpNodePosX_M[index],
						divAuxData.tmpNodePosY_M[index], 0);
				CVector centerToPosDir = memPos - oldCellCenter;  
				CVector centerToPosUnit = centerToPosDir.getUnitVector();
				double dotProduct = centerToPosUnit * splitDir;
				tmpData.val = dotProduct; 
				tmpData.vec = memPos;
				tmpData.type=divAuxData.tmpNodeType[index] ; 

				if (divAuxData.tmpNodeType[index] == lateralA){
					tmp2.push_back(tmpData);
				}
				else if (divAuxData.tmpNodeType[index] == lateralB){
					tmp1.push_back(tmpData);
				}
				else if (divAuxData.tmpNodeType[index] == apical1){
					double crossProduct_withApical = centerToPosDir.x*oldCellCenter2ApicalVec[1] - centerToPosDir.y*oldCellCenter2ApicalVec[0];
					// double crossProduct_withBasal = centerToPosDir.x*oldCellCenter2BasalVec[1] - centerToPosDir.y*oldCellCenter2BasalVec[0];
					if (crossProduct_withApical >= 0){
						tmp2.push_back(tmpData);
					}
					else if (crossProduct_withApical < 0){
						tmp1.push_back(tmpData);
					}
					else {std::cout<<"None of the condition to determine if an apical node belongs to daughter cell or mother cell is met! Something is wrong!"<<std::endl;}
				}
				else if (divAuxData.tmpNodeType[index] == basal1){
					// double crossProduct_withApical = centerToPosDir.x*oldCellCenter2ApicalVec[1] - centerToPosDir.y*oldCellCenter2ApicalVec[0] ;
					double crossProduct_withBasal = centerToPosDir.x*oldCellCenter2BasalVec[1] - centerToPosDir.y*oldCellCenter2BasalVec[0];
					if (crossProduct_withBasal < 0){
						tmp2.push_back(tmpData);
					}
					else if (crossProduct_withBasal >= 0){
						tmp1.push_back(tmpData);
					}
					else {std::cout<<"None of the condition to determine if a basal node belongs to daughter cell or mother cell is met! Something is wrong!"<<std::endl;}
				}
				else{
					std::cout<<"Active membrane node not assigned with any node type is present! Something is wrong!"<<std::endl;
				}

				/*CVector memPos(divAuxData.tmpNodePosX_M[index],
						divAuxData.tmpNodePosY_M[index], 0);
				CVector centerToPosDir = memPos - oldCellCenter;  // Ali it should be center of cells
				CVector centerToPosUnit = centerToPosDir.getUnitVector();
				CVector crossProduct = Cross(centerToPosDir, splitDir);
				double dotProduct = centerToPosUnit * splitDir;//This is the original way to determine which node belongs to which node belongs to mother and daughter cell
				// double dotProduct = centerToPosUnit.getUnitVector() * splitDir.getUnitVector();
				// double dotProduct2 = centerToPosDir * divDir;
				tmpData.val = dotProduct; // for sorting the membrane nodes
				tmpData.vec = memPos;
				tmpData.type=divAuxData.tmpNodeType[index] ; 
				if (crossProduct.z >= 0) {
					// counter-cloce wise
					tmp1.push_back(tmpData);
				} else {
					// cloce wise
					tmp2.push_back(tmpData);
				}*/ //This is the original way to determine which node belongs to which node belongs to mother and daughter cell //Kevin

				// if (dotProduct2 >= 0){
				// 	tmp1.push_back(tmpData);
				// }
				// else{
				// 	tmp2.push_back(tmpData);
				// }
			}
		} else {// shrink the internal nodes around the internal node center
			if (divAuxData.tmpIsActive_M[index] == true) {
				CVector internalPos(divAuxData.tmpNodePosX_M[index],
						divAuxData.tmpNodePosY_M[index], 0);
				CVector centerToPosDir = internalPos - oldIntCenter;  // center of nucleus is more biological
				CVector shrinkedPos = centerToPosDir * shrinkRatio + oldIntCenter;

		       // CVector unitDivDir = divDir.getUnitVector(); // Ali 
			//	double  AmpTanget=centerToPosDir*unitDivDir ;  // Ali dot product of two vectors
			//	double  shrinkedAmpTanget=shrinkRatio*AmpTanget; // multiply two doubles //Ali

			//	CVector TangetVShrink=unitDivDir*shrinkedAmpTanget; // shrink the tanget component //Ali
			//	CVector TangetV=      unitDivDir*        AmpTanget; // get the tanget component to compute the normal vector  //Ali
			//	CVector  NormV=centerToPosDir-TangetV ;             // compute the normal vector Ali

			//	CVector  polarShrinkedPos=NormV+TangetVShrink ;  // summation of shrinked tanget and as previous vector in the normal direction to division axis//Ali
			//	CVector  updatedV=polarShrinkedPos+oldCenter ; //Ali 

				double dotProduct = centerToPosDir * divDir; 
				//double dotProduct = polarShrinkedPos   * divDir; //Ali 
				if (dotProduct >= 0) {
					if (divAuxData.tmp1IntnlVec.size() >= 24){
						divAuxData.tmp2IntnlVec.push_back(shrinkedPos);
					}
					else{
						divAuxData.tmp1IntnlVec.push_back(shrinkedPos);
					}
				} else {
					if (divAuxData.tmp2IntnlVec.size() >= 24){
						divAuxData.tmp1IntnlVec.push_back(shrinkedPos);
					}
					else{
						divAuxData.tmp2IntnlVec.push_back(shrinkedPos);
					}
				}
			}
		}
	}
	for (int i = 0; i < tmp1.size(); i++){
		// std::cout<<"tmp1 ["<<i<<"] "<<tmp1[i].vec.x<<" "<<tmp1[i].vec.y<<" "<<tmp1[i].val<<" "<<tmp1[i].type<< std::endl;
	}
	for (int i = 0; i < tmp2.size(); i++){
		// std::cout<<"tmp2 ["<<i<<"] "<<tmp2[i].vec.x<<" "<<tmp2[i].vec.y<<" "<<tmp2[i].val<<" "<<tmp2[i].type<< std::endl;
	}

	int targetId_tmp2;
	for (int i = 0; i < tmp2.size(); i++){
		int iNext = i+1;
		if (i == (tmp2.size()-1)){
			iNext = 0;
		}
		// std::cout<<"tmp2 i = "<<i<<std::endl;
		// std::cout<<"temp2 iNext = "<<iNext<<std::endl;
		double dotProduct = ((tmp2[i].vec.x  - oldCellCenter.x) * (tmp2[iNext].vec.x  - oldCellCenter.x)) + (tmp2[i].vec.y - oldCellCenter.y) * (tmp2[iNext].vec.y - oldCellCenter.y) ;
		if (dotProduct < 0){
			targetId_tmp2 = i;
			std::cout<<"i for dotProduct < 0, tmp2 = "<<i<<std::endl;
		}
	}
	// std::cout<<"tmp2 before rearranging"<<std::endl;
	// for (int i = 0; i < tmp2.size(); i++){
	// 	std::cout<<tmp2[i].vec.x<<" "<<tmp2[i].vec.y<<std::endl;
	// }
	if (targetId_tmp2 != (tmp2.size()-1)){
		vector<VecValT> tmp_tmp2;
		int currentId = targetId_tmp2+1;
		for (int i = 0; i<tmp2.size(); i++){
			tmp_tmp2.push_back(tmp2[currentId]);
			currentId += 1;
			if (currentId >= tmp2.size()){
				currentId = currentId - tmp2.size();
			}
		}
		tmp2.clear();
		for (int i = 0; i < tmp_tmp2.size(); i++){
			tmp2.push_back(tmp_tmp2[i]);
		}
	}
	else{
		std::cout<<"No need to reorganize the order of tmp2Membr before adding new nodes"<<std::endl;
	}

	int targetId_tmp1;
	for (int i = 0; i < tmp1.size(); i++){
		int iNext = i+1;
		if (i == (tmp1.size()-1)){
			iNext = 0;
		}
		double dotProduct = (tmp1[i].vec.x - oldCellCenter.x) * (tmp1[iNext].vec.x - oldCellCenter.x) + (tmp1[i].vec.y - oldCellCenter.y) * (tmp1[iNext].vec.y - oldCellCenter.y) ;
		if (dotProduct < 0){
			targetId_tmp1 = i;
			std::cout<<"i for dotProduct < 0, tmp1 = "<<i<<std::endl;
		}
	}
	if (targetId_tmp1 != (tmp1.size()-1)){
		vector<VecValT> tmp_tmp1;
		int currentId = targetId_tmp1+1;
		for (int i = 0; i<tmp1.size(); i++){
			tmp_tmp1.push_back(tmp1[currentId]);
			currentId += 1;
			if (currentId >= tmp1.size()){
				currentId = currentId - tmp1.size();
			}
		}
		tmp1.clear();
		for (int i = 0; i < tmp_tmp1.size(); i++){
			tmp1.push_back(tmp_tmp1[i]);
		}
	}
	else{
		std::cout<<"No need to reorganize the order of tmp1Membr before adding new nodes"<<std::endl;
	}

	// oldCellCenter2BasalVec[0] = divAuxData.tmpBasalLoc[0] - oldCellCenter.x;
	// oldCellCenter2BasalVec[1] = divAuxData.tmpBasalLoc[1] - oldCellCenter.y;
	// double dotProduct_withBasal = -10000.0;
	// int targetId_tmp2 = -100;
	// for (int i = 0; i < tmp2.size(); i++){
	// 	CVector tmp2_vec = tmp2[i].vec - oldCellCenter;
	// 	double tmp2_dotProduct_withBasal = tmp2_vec.x*oldCellCenter2BasalVec[0] + tmp2_vec.y*oldCellCenter2BasalVec[1];
	// 	if (tmp2_dotProduct_withBasal >= dotProduct_withBasal){
	// 		dotProduct_withBasal = tmp2_dotProduct_withBasal;
	// 		targetId_tmp2 = i;
	// 	}
	// }
	// vector<VecValT> tmp_tmp2;
	// int currentId = targetId_tmp2;
	// for (int i = 0; i<tmp2.size(); i++){
	// 	tmp_tmp2.push_back(tmp2[currentId]);
	// 	currentId += 1;
	// 	if (currentId >= tmp2.size()){
	// 		currentId = currentId - tmp2.size();
	// 	}
	// }
	// tmp2.clear();
	// for (int i = 0; i < tmp_tmp2.size(); i++){
	// 	tmp2.push_back(tmp_tmp2[i]);
	// }
}

void SceCells::calCellArea() {
	thrust::counting_iterator<uint> iBegin(0), iBegin2(0);
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));

	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(
									thrust::make_permutation_iterator(
											cellInfoVecs.activeMembrNodeCounts.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell)),
									make_transform_iterator(iBegin,
											ModuloFunctor(maxAllNodePerCell)),
									make_permutation_iterator(
											cellInfoVecs.centerCoordX.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))),
									make_permutation_iterator(
											cellInfoVecs.centerCoordY.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))))),
					CalTriArea(maxAllNodePerCell, nodeIsActiveAddr,
							nodeLocXAddr, nodeLocYAddr)),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.cellAreaVec.begin(), thrust::equal_to<uint>(),
			thrust::plus<double>());
}


   //AAMIRI added to calculate Perimeter of each cell
 void SceCells::calCellPerim() {
 	thrust::counting_iterator<uint> iBegin(0), iBegin2(0);
 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
 			* allocPara_m.maxAllNodePerCell;
 	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
 
 	double* nodeLocXAddr = thrust::raw_pointer_cast(
 			&(nodes->getInfoVecs().nodeLocX[0]));
 	double* nodeLocYAddr = thrust::raw_pointer_cast(
 			&(nodes->getInfoVecs().nodeLocY[0]));
 	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
 			&(nodes->getInfoVecs().nodeIsActive[0]));
 
 	thrust::reduce_by_key(
 			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell)),
 			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell))
 					+ totalNodeCountForActiveCells,
 			thrust::make_transform_iterator(
 					thrust::make_zip_iterator(
 							thrust::make_tuple(
 									thrust::make_permutation_iterator(
 											cellInfoVecs.activeMembrNodeCounts.begin(),
 											make_transform_iterator(iBegin,
 													DivideFunctor(
 															maxAllNodePerCell))),
 									make_transform_iterator(iBegin,
 											DivideFunctor(maxAllNodePerCell)),
 									make_transform_iterator(iBegin,
 											ModuloFunctor(maxAllNodePerCell)),
 									make_permutation_iterator(
 											cellInfoVecs.centerCoordX.begin(),
 											make_transform_iterator(iBegin,
 													DivideFunctor(
 															maxAllNodePerCell))),
 									make_permutation_iterator(
 											cellInfoVecs.centerCoordY.begin(),
 											make_transform_iterator(iBegin,
 													DivideFunctor(
 															maxAllNodePerCell))))),
 					CalPerim(maxAllNodePerCell, nodeIsActiveAddr,
 							nodeLocXAddr, nodeLocYAddr)),
 			cellInfoVecs.cellRanksTmpStorage.begin(),
 			cellInfoVecs.cellPerimVec.begin(), thrust::equal_to<uint>(),
 			thrust::plus<double>());
 }
 
   //Ali added to calculate pressure of each cell
 void SceCells::calCellPressure() {
 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
 			* allocPara_m.maxAllNodePerCell;
 
 	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
 	thrust::counting_iterator<uint> iBegin(0) ;

			// if the membrane node is not active or if this is an internal node it will automatically add with zero.
			thrust::reduce_by_key(
									make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
									make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeF_MI_M_N.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.sumF_MI_M_N.begin());  

			thrust::reduce_by_key(
									make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
									make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().lagrangeFN.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.sumLagrangeFN.begin());  

			thrust:: transform(cellInfoVecs.sumF_MI_M_N.begin(), cellInfoVecs.sumF_MI_M_N.begin()+allocPara_m.currentActiveCellCount, 
			                   cellInfoVecs.sumLagrangeFN.begin(), cellInfoVecs.cellPressure.begin(),thrust::plus<float>()) ;
			thrust:: transform(cellInfoVecs.cellPressure.begin(), cellInfoVecs.cellPressure.begin()+allocPara_m.currentActiveCellCount, 
			                   cellInfoVecs.cellPerimVec.begin(), cellInfoVecs.cellPressure.begin(),thrust::divides<float>()) ;



 }
 












CellsStatsData SceCells::outputPolyCountData() {
       
        cout << " I am at begining of outpolycount"<< std::flush  ; 
	std::cout.flush();
       double sumX,sumY,cntr_X_Domain,cntr_Y_Domain ; 
       int BdryApproach ; 
       BdryApproach=1 ; 
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
        cout << " I am before cells area"<< endl ; 
	calCellArea();
	for (int k = 0; k < cellInfoVecs.cellAreaVec.size(); k++){
		if (cellInfoVecs.cellAreaVec[k] < 0){
			cellInfoVecs.cellAreaVec[k] = -1.0*cellInfoVecs.cellAreaVec[k];
		}
	} // !!!!!NOTE!!!!!! this is currently an ad hoc way to make sure area of each triangle is positive.
        cout << " I am after cells area" << endl ;

        calCellPerim();//AAMIRI
		calCellPressure() ; // Ali 
    	//computeBasalLoc();  //Ali  we call it here to compute the length of the cells
	CellsStatsData result;

        cout << " I am after result" << endl ; 
	uint bdryCriteria =
			globalConfigVars.getConfigValue("BdryCellCriteria").toInt();
	// already on host; no need to call thrust::copy
	// thrust::host_vector<int> adhIndxHost =
	// 		nodes->getInfoVecs().nodeAdhIndxHostCopy;
	thrust::host_vector<int> adhIndxHost =
			nodes->getInfoVecs().nodeAdhereIndex;

	thrust::host_vector<double> growthProVecHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			growthProVecHost.begin());
	thrust::host_vector<double> growthProMembrVecHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.membrGrowProgress.begin(),
			cellInfoVecs.membrGrowProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			growthProMembrVecHost.begin());
	thrust::host_vector<uint> activeMembrNodeCountHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.activeMembrNodeCounts.begin(),
			cellInfoVecs.activeMembrNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount,
			activeMembrNodeCountHost.begin());
	thrust::host_vector<uint> activeIntnlNodeCountHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.activeIntnlNodeCounts.begin(),
			cellInfoVecs.activeIntnlNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount,
			activeIntnlNodeCountHost.begin());
//////////////
	thrust::host_vector<double> centerCoordXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> centerCoordYHost(
			allocPara_m.currentActiveCellCount);
	
	thrust::copy(cellInfoVecs.centerCoordX.begin(),
			cellInfoVecs.centerCoordX.begin()
					+ allocPara_m.currentActiveCellCount,
			centerCoordXHost.begin());
	thrust::copy(cellInfoVecs.centerCoordY.begin(),
			cellInfoVecs.centerCoordY.begin()
					+ allocPara_m.currentActiveCellCount,
			centerCoordYHost.begin());
/////////////

// std::cout<<"ISOLATE 1"<<std::endl;

///////////////
	thrust::host_vector<double> apicalLocXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> apicalLocYHost(
			allocPara_m.currentActiveCellCount);

	thrust::copy(cellInfoVecs.apicalLocX.begin(),
			cellInfoVecs.apicalLocX.begin()
					+ allocPara_m.currentActiveCellCount,
			apicalLocXHost.begin());
	thrust::copy(cellInfoVecs.apicalLocY.begin(),
			cellInfoVecs.apicalLocY.begin()
					+ allocPara_m.currentActiveCellCount,
			apicalLocYHost.begin());

///////////

// std::cout<<"ISOLATE 2"<<std::endl;

///////////////
	thrust::host_vector<double> basalLocXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> basalLocYHost(
			allocPara_m.currentActiveCellCount);

	thrust::copy(cellInfoVecs.basalLocX.begin(),
			cellInfoVecs.basalLocX.begin()
					+ allocPara_m.currentActiveCellCount,
			basalLocXHost.begin());
	thrust::copy(cellInfoVecs.basalLocY.begin(),
			cellInfoVecs.basalLocY.begin()
					+ allocPara_m.currentActiveCellCount,
			basalLocYHost.begin());

///////////


// std::cout<<"ISOLATE 3"<<std::endl;

//////
	thrust::host_vector<double> InternalAvgXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> InternalAvgYHost(
			allocPara_m.currentActiveCellCount);

	thrust::host_vector<double> InternalMaxXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> InternalMaxYHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> InternalMinXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> InternalMinYHost(
			allocPara_m.currentActiveCellCount);

	thrust::copy(cellInfoVecs.InternalAvgX.begin(),
			cellInfoVecs.InternalAvgX.begin()
					+ allocPara_m.currentActiveCellCount,
			InternalAvgXHost.begin());
	thrust::copy(cellInfoVecs.InternalAvgY.begin(),
			cellInfoVecs.InternalAvgY.begin()
					+ allocPara_m.currentActiveCellCount,
			InternalAvgYHost.begin());
/////

// std::cout<<"ISOLATE 4"<<std::endl;

	thrust::host_vector<double> cellAreaHost(
			allocPara_m.currentActiveCellCount);

        thrust::host_vector<double> cellPerimHost(
 			allocPara_m.currentActiveCellCount);//AAMIRI

        thrust::host_vector<double> cellPressureHost(
 			allocPara_m.currentActiveCellCount);//Ali

	thrust::copy(cellInfoVecs.cellAreaVec.begin(),
			cellInfoVecs.cellAreaVec.begin()
					+ allocPara_m.currentActiveCellCount, cellAreaHost.begin());

        thrust::copy(cellInfoVecs.cellPerimVec.begin(),
 			cellInfoVecs.cellPerimVec.begin()
					+ allocPara_m.currentActiveCellCount, cellPerimHost.begin());//AAMIRI

        thrust::copy(cellInfoVecs.cellPressure.begin(),
 			cellInfoVecs.cellPressure.begin()
					+ allocPara_m.currentActiveCellCount, cellPressureHost.begin());//Ali

// std::cout<<"ISOLATE 4"<<std::endl;
//LOOKS LIKE AN ERROR POP UP IN THIS SECTION..............
        sumX=0 ; 
        sumY=0 ; 
	for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
		CellStats cellStatsData;
		cellStatsData.cellGrowthProgress = growthProVecHost[i];
		cellStatsData.cellRank = i;
		bool isBdry = false;
		std::set<int> neighbors;
		std::vector<int> neighborsV; //Ali
                int neighborStrength[10]; //Ali
		int continousNoAdh = 0;
                map <int, int> cellAndNeighborRank ;  //Ali
		// std::cout << "printing adhesion indicies ";
                //for (int ii=0 ; ii<neighborStrength.size() ; ii++){
                for (int ii=0 ; ii< 10; ii++){ //Ali
                      
                  neighborStrength[ii]=0  ;
                }
                          
                cellAndNeighborRank.clear();  //Ali
		// std::cout << "printing adhesion indicies 2";
		for (uint j = 0; j < activeMembrNodeCountHost[i]; j++) {
			uint index = i * allocPara_m.maxAllNodePerCell + j;
			//std::cout << adhIndxHost[index] << ",";
                        
			if (adhIndxHost[index] != -1) {
				uint adhCellRank = adhIndxHost[index]
						/ allocPara_m.maxAllNodePerCell;
				//std::cout << adhCellRank << " ";
				neighbors.insert(adhCellRank);
                                 map <int, int>:: iterator iteratorMap=cellAndNeighborRank.find(adhCellRank); //Ali
                                 if (iteratorMap==cellAndNeighborRank.end()) {  //Ali
                                   int NewneighborRank= neighbors.size()-1; //Ali
                                   cellAndNeighborRank[adhCellRank]=NewneighborRank; //Ali
                                   neighborStrength[NewneighborRank]=neighborStrength[NewneighborRank]+1 ; //Ali
				   neighborsV.push_back(adhCellRank); //Ali
                                   }
                                 else {   //Ali
                                   int oldNeighborRank=(*iteratorMap).second ; 
                                   neighborStrength[oldNeighborRank]=neighborStrength[oldNeighborRank]+1 ; //Ali
                                 }      
				continousNoAdh = 0;
				// std::cout << "printing adhesion indicies 3";
			} else {
				continousNoAdh = continousNoAdh + 1;
				if (continousNoAdh > bdryCriteria) {
					isBdry = true;
				}
				// std::cout << "printing adhesion indicies 4";
			}
			if (j == activeMembrNodeCountHost[i] - 1
					&& adhIndxHost[index] == -1) {
				int k = 0;
				uint indexNew;
				while (k < activeMembrNodeCountHost[i] - 1) {
					indexNew = i * allocPara_m.maxAllNodePerCell + k;
					if (adhIndxHost[indexNew] == -1) {
						continousNoAdh = continousNoAdh + 1;
						if (continousNoAdh > bdryCriteria) {
							isBdry = true;
						}
						k++;
					} else {
						break;
					}
				}
				// std::cout << "printing adhesion indicies 5";
			}

		}

		// std::cout<<"ISOLATE 4.5"<<std::endl; 
                
		cellStatsData.isBdryCell = isBdry;
		cellStatsData.numNeighbors = neighbors.size();
		cellStatsData.currentActiveMembrNodes = activeMembrNodeCountHost[i];
		cellStatsData.currentActiveIntnlNodes = activeIntnlNodeCountHost[i];
		cellStatsData.neighborVec = neighbors;
		cellStatsData.neighborVecV = neighborsV; //Ali
                for (int iiii=0; iiii<10 ; iiii++){
		cellStatsData.cellNeighborStrength[iiii] = neighborStrength[iiii];
                 }    //Ali
		cellStatsData.membrGrowthProgress = growthProMembrVecHost[i];
		cellStatsData.cellCenter = CVector(centerCoordXHost[i],
				centerCoordYHost[i], 0);

		cellStatsData.cellApicalLoc= CVector(apicalLocXHost[i],
				apicalLocYHost[i], 0);  //Ali

		cellStatsData.cellBasalLoc= CVector(basalLocXHost[i],
				basalLocYHost[i], 0);  //Ali
		cellStatsData.cellNucleusLoc = CVector(InternalAvgXHost[i],
				InternalAvgYHost[i], 0);  // Ali
		cellStatsData.cellNucleusLocMax = CVector(InternalMaxXHost[i],
				InternalMaxYHost[i], 0);
		cellStatsData.cellNucleusLocMin = CVector(InternalMinXHost[i],
				InternalMinYHost[i], 0);
		cellStatsData.cellArea = cellAreaHost[i];
        cellStatsData.cellPerim = cellPerimHost[i];//AAMIRI
        cellStatsData.cellPressure = cellPressureHost[i];//Ali
		result.cellsStats.push_back(cellStatsData);
                sumX=sumX+cellStatsData.cellCenter.x ; 
                sumY=sumY+cellStatsData.cellCenter.y ;
                
	}

	// std::cout<<"ISOLATE 5"<<std::endl;
//Ali
        if (BdryApproach==2) {  
          cout << "sumX=" << sumX << endl ; 
          cout << "sumY=" << sumY << endl ; 
          cntr_X_Domain=sumX/result.cellsStats.size() ; 
          cntr_Y_Domain=sumY/result.cellsStats.size() ;  
          cout << "cntr_X=" << cntr_X_Domain << endl ; 
          cout << "cntr_Y=" << cntr_Y_Domain << endl ;

          double R_Max ;
          double Distance ;
          R_Max=0 ;  
	  for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
            Distance=sqrt( pow(centerCoordXHost[i]-cntr_X_Domain,2) +pow(centerCoordYHost[i]-cntr_Y_Domain,2) ) ; 
            if (Distance > R_Max) {
              R_Max=Distance ; 
            }
          }
        
          cout << "R_Max=" << R_Max << endl ;

	  for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
            Distance=sqrt( pow(centerCoordXHost[i]-cntr_X_Domain,2) +pow(centerCoordYHost[i]-cntr_Y_Domain,2) ) ; 
            if (Distance > 0.9* R_Max) {
	      result.cellsStats[i].isBdryCell = true;
              cout << "isBdryCell"<< i<< endl ; 
            }
            else {
	      result.cellsStats[i].isBdryCell = false;
              cout << "isNormalCell"<< i << endl ; 
            }
          }
        }
		// std::cout<<"ISOLATE 6"<<std::endl;
        	return result;
}

SingleCellData SceCells::OutputStressStrain() {

	SingleCellData result ; 
    vector <double> nodeExtForceXHost; 
    vector <double> nodeExtForceYHost; 

	nodeExtForceXHost.resize(totalNodeCountForActiveCells); 
	nodeExtForceYHost.resize(totalNodeCountForActiveCells); 
	thrust::copy ( nodes->getInfoVecs().nodeExtForceX.begin(),
	               nodes->getInfoVecs().nodeExtForceX.begin()+ totalNodeCountForActiveCells,
				   nodeExtForceXHost.begin()); 
	thrust::copy ( nodes->getInfoVecs().nodeExtForceY.begin(),
	               nodes->getInfoVecs().nodeExtForceY.begin()+ totalNodeCountForActiveCells,
				   nodeExtForceYHost.begin()); 
         // There is a compiling issue with using count_if on GPU. 
    int numPositiveForces     =     count_if(nodeExtForceXHost.begin(),nodeExtForceXHost.end(),isGreaterZero() ) ;
    double totalExtPositiveForce =accumulate(nodeExtForceXHost.begin(),nodeExtForceXHost.end(),0.0, SumGreaterZero() ) ;
	cout << "number of positive external forces are=" <<numPositiveForces<<endl ; 
	cout << "Total external forces are=" <<totalExtPositiveForce<<endl ; 


	//thrust::device_vector<double>::iterator  
	double MinX=*thrust::min_element(nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount,
                                     nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
    //thrust::device_vector<double>::iterator 
	double MaxX=*thrust::max_element(nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount,
                                     nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;

    result.Cells_Extrem_Loc[0]=MinX; 
    result.Cells_Extrem_Loc[1]=MaxX; 
    result.F_Ext_Out=totalExtPositiveForce ; 
    return result ; 

}

__device__ bool bigEnough(double& num) {
	if (num > minDivisor) {
		return true;
	} else {
		return false;
	}
}

__device__ double cross_Z(double vecA_X, double vecA_Y, double vecB_X,
		double vecB_Y) {
	return vecA_X * vecB_Y - vecA_Y * vecB_X;
}
/*
__device__ double calBendMulti(double& angle, uint activeMembrCt) {
	double equAngle = PI - PI / activeMembrCt;
	return bendCoeff * (angle - equAngle);
}
*/

//AAMIRI
__device__ double calBendMulti_Mitotic(double& angle, uint activeMembrCt, double& progress, double mitoticCri) {

	//double equAngle = PI - PI / activeMembrCt;
	double equAngle = PI ; // - PI / activeMembrCt;
	
	if (progress <= mitoticCri){
		return bendCoeff * (angle - equAngle);}
	else{
		return (angle - equAngle)*(bendCoeff + (bendCoeff_Mitotic - bendCoeff) * (progress - mitoticCri)/(1.0 - mitoticCri));
	}
}

__device__ double CalMembrBendSpringEnergy(double& angle, uint activeMembrCt, double& progress, double mitoticCri) {

	double equAngle = PI - PI / activeMembrCt;
	
	if (progress <= mitoticCri){
		return ( 0.5*bendCoeff * (angle - equAngle)*(angle - equAngle) );
		}
	else{
		return ( 0.5*bendCoeff *  (angle - equAngle)*(angle - equAngle)*
		                         (bendCoeff + (bendCoeff_Mitotic - bendCoeff) * (progress - mitoticCri)/(1.0 - mitoticCri)) );
	}
}




void SceCells::applySceCellDisc_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));
	bool* isEnteringMitotic = thrust::raw_pointer_cast(
			&(cellInfoVecs.isEnteringMitotic[0]));

	//double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
	//		- growthAuxData.prolifDecay
	//				* (growthAuxData.grthProgrEndCPU
	//						- growthAuxData.grthPrgrCriVal_M_Ori);

	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeIntnlNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeIntnlNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							   nodes->getInfoVecs().nodeVelY.begin(),
							   nodes->getInfoVecs().nodeF_MI_M_x.begin(),  //Ali added for cell pressure calculation 
							   nodes->getInfoVecs().nodeF_MI_M_y.begin(),// ALi added for cell pressure calculation
							   nodes->getInfoVecs().nodeIIEnergy.begin(),
							   nodes->getInfoVecs().nodeIMEnergy.begin())),
			AddSceCellForce(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,
					nodeLocYAddr, nodeIsActiveAddr, grthPrgrCriVal_M, isEnteringMitotic));
/*
		for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++) {
				cout << "for cell rank "<<i<< " cell apical location is " << cellInfoVecs.apicalLocX[i] <<" , " <<cellInfoVecs.apicalLocY[i] << endl ; 
				cout << "for cell rank "<<i<< " cell nucleus distance from apical is " << cellInfoVecs.nucDesireDistApical[i]  << endl ; 
		}

*/
}


// void SceCells::applyMembContraction() {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;
// 	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
// 	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
// 	thrust::counting_iterator<uint> iBegin2(0);

// 	double* nodeLocXAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeLocX[0]));
// 	double* nodeLocYAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeLocY[0]));
// 	MembraneType1* nodeTypeAddr=thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().memNodeType1[0]));
// 	int* nodeMemMirrorIndexAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeMemMirrorIndex[0]));

// 			thrust::transform(
// 				thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.nucDesireDistApical.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.activeMembrNodeCounts.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocX.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocY.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell)),
// 									make_transform_iterator(iBegin2,
// 											ModuloFunctor(maxAllNodePerCell)),
// 									nodes->getInfoVecs().nodeIsActive.begin(),
// 									nodes->getInfoVecs().nodeVelX.begin(),
// 									nodes->getInfoVecs().nodeVelY.begin(),
// 									nodes->getInfoVecs().memNodeType1.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.nucDesireDistApical.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.activeMembrNodeCounts.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocX.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocY.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 									make_transform_iterator(iBegin2,
// 										DivideFunctor(maxAllNodePerCell)),
// 									make_transform_iterator(iBegin2,
// 										ModuloFunctor(maxAllNodePerCell)),
// 									nodes->getInfoVecs().nodeIsActive.begin(),
// 									nodes->getInfoVecs().nodeVelX.begin(),
// 									nodes->getInfoVecs().nodeVelY.begin(),
// 									nodes->getInfoVecs().memNodeType1.begin()))
// 									+ totalNodeCountForActiveCells,
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 								nodes->getInfoVecs().nodeVelX.begin(),
// 							   	nodes->getInfoVecs().nodeVelY.begin(),
// 							    nodes->getInfoVecs().nodeF_MM_C_X.begin(),   
// 							    nodes->getInfoVecs().nodeF_MM_C_Y.begin(),
// 							    nodes->getInfoVecs().nodeContractEnergyT.begin(),
// 							    nodes->getInfoVecs().basalContractPair.begin())),
// 			AddMemContractForce(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,nodeLocYAddr, nodeTypeAddr,nodeMemMirrorIndexAddr));

	
// }
// void SceCells::applyMembContraction(double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_Division, double distFromNucleus_max, double distFromNucleus_min, double percentage_before_timeRatio_Crit_Division_scaling) {
// 	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
// 			* allocPara_m.maxAllNodePerCell;
// 	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
// 	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
// 	thrust::counting_iterator<uint> iBegin2(0);

// 	double* nodeLocXAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeLocX[0]));
// 	double* nodeLocYAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeLocY[0]));
// 	// double* nodeLocZAddr = thrust::raw_pointer_cast(
// 	// 		&(nodes->getInfoVecs().nodeLocZ[0]));
// 	/*double* ActomyosinMultipBasal = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeActomyosinMultip_basal[0]));
// 	double* ActomyosinMultipApical = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeActomyosinMultip_apical[0]));*/
// 	MembraneType1* nodeTypeAddr=thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().memNodeType1[0]));
// 	int* nodeMemMirrorIndexAddr = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().nodeMemMirrorIndex[0]));
// 	double* contractActomyo_multip = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().contractActomyo_multip[0]));
// 	double* contractActomyo_multip_apical = thrust::raw_pointer_cast(
// 			&(nodes->getInfoVecs().contractActomyo_multip_apical[0]));
// 	bool* isEnteringMitotic = thrust::raw_pointer_cast(
// 			&(cellInfoVecs.isEnteringMitotic[0]));
// 	double* contractileSpringGrowthProgress = thrust::raw_pointer_cast(
// 			&(cellInfoVecs.contractileSpringGrowthProgress[0]));
// 	double* distFromNucleus_normal = thrust::raw_pointer_cast(
// 			&(cellInfoVecs.distFromNucleus_normal[0]));
// 	double* distFromNucleus_normal_apical = thrust::raw_pointer_cast(
// 			&(cellInfoVecs.distFromNucleus_normal_apical[0]));

// 			thrust::transform(
// 				thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.nucDesireDistApical.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.activeMembrNodeCounts.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocX.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocY.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell)),
// 									make_transform_iterator(iBegin2,
// 											ModuloFunctor(maxAllNodePerCell)),
// 									nodes->getInfoVecs().nodeIsActive.begin(),
// 									nodes->getInfoVecs().nodeVelX.begin(),
// 									nodes->getInfoVecs().nodeVelY.begin(),
// 									nodes->getInfoVecs().memNodeType1.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.nucDesireDistApical.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.activeMembrNodeCounts.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocX.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 							thrust::make_permutation_iterator(
// 									cellInfoVecs.apicalLocY.begin(),
// 									make_transform_iterator(iBegin2,
// 											DivideFunctor(maxAllNodePerCell))),
// 									make_transform_iterator(iBegin2,
// 										DivideFunctor(maxAllNodePerCell)),
// 									make_transform_iterator(iBegin2,
// 										ModuloFunctor(maxAllNodePerCell)),
// 									nodes->getInfoVecs().nodeIsActive.begin(),
// 									nodes->getInfoVecs().nodeVelX.begin(),
// 									nodes->getInfoVecs().nodeVelY.begin(),
// 									nodes->getInfoVecs().memNodeType1.begin()))
// 									+ totalNodeCountForActiveCells,
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 								nodes->getInfoVecs().nodeVelX.begin(),
// 							   	nodes->getInfoVecs().nodeVelY.begin(),
// 							    nodes->getInfoVecs().nodeF_MM_C_X.begin(),   
// 							    nodes->getInfoVecs().nodeF_MM_C_Y.begin(),
// 							    nodes->getInfoVecs().nodeContractEnergyT.begin(),
// 							    nodes->getInfoVecs().basalContractPair.begin())),
// 			// AddMemContractForce(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,nodeLocYAddr, ActomyosinMultipBasal,ActomyosinMultipApical, nodeTypeAddr,nodeMemMirrorIndexAddr, 
// 			// 						timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_Division));
// 			AddMemContractForce_tmp(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,nodeLocYAddr, nodeTypeAddr,nodeMemMirrorIndexAddr, 
// 									timeRatio,  contractActomyo_multip, contractActomyo_multip_apical, distFromNucleus_max, distFromNucleus_min, distFromNucleus_normal, distFromNucleus_normal_apical, percentage_before_timeRatio_Crit_Division_scaling, isEnteringMitotic, contractileSpringGrowthProgress));

	
// }
void SceCells::applyMembContraction2(double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_Division, double distFromNucleus_max, double distFromNucleus_min, double mitoRndActomyoStrengthScaling) {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	thrust::counting_iterator<uint> iBegin2(0);

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	// double* nodeLocZAddr = thrust::raw_pointer_cast(
	// 		&(nodes->getInfoVecs().nodeLocZ[0]));
	/*double* ActomyosinMultipBasal = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeActomyosinMultip_basal[0]));
	double* ActomyosinMultipApical = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeActomyosinMultip_apical[0]));*/
	MembraneType1* nodeTypeAddr=thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().memNodeType1[0]));
	int* nodeMemMirrorIndexAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeMemMirrorIndex[0]));
	// double* contractActomyo_multip = thrust::raw_pointer_cast(
	// 		&(nodes->getInfoVecs().contractActomyo_multip[0]));
	// double* contractActomyo_multip_apical = thrust::raw_pointer_cast(
	// 		&(nodes->getInfoVecs().contractActomyo_multip_apical[0]));
	double* contractActomyo_multip = thrust::raw_pointer_cast(
			&(cellInfoVecs.contractActomyo_multip_perCell[0]));
	double* contractActomyo_multip_apical = thrust::raw_pointer_cast(
			&(cellInfoVecs.contractActomyo_multip_apical_perCell[0]));
	bool* isEnteringMitotic = thrust::raw_pointer_cast(
			&(cellInfoVecs.isEnteringMitotic[0]));
	double* contractileSpringGrowthProgress = thrust::raw_pointer_cast(
			&(cellInfoVecs.contractileSpringGrowthProgress[0]));
	double* distFromNucleus_normal = thrust::raw_pointer_cast(
			&(cellInfoVecs.distFromNucleus_normal[0]));
	double* distFromNucleus_normal_apical = thrust::raw_pointer_cast(
			&(cellInfoVecs.distFromNucleus_normal_apical[0]));
	double* individualCellHeight = thrust::raw_pointer_cast(
			&(cellInfoVecs.individualCellHeight[0]));
	double* distFromBasalLoc = thrust::raw_pointer_cast(
			&(cellInfoVecs.distFromBasalLoc[0]));
	double* distFromApicalLoc = thrust::raw_pointer_cast(
			&(cellInfoVecs.distFromApicalLoc[0]));
	

			thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.nucDesireDistApical.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.apicalLocX.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.apicalLocY.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell)),
									make_transform_iterator(iBegin2,
											ModuloFunctor(maxAllNodePerCell)),
									nodes->getInfoVecs().nodeIsActive.begin(),
									nodes->getInfoVecs().nodeVelX.begin(),
									nodes->getInfoVecs().nodeVelY.begin(),
									nodes->getInfoVecs().memNodeType1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.nucDesireDistApical.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.apicalLocX.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.apicalLocY.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor(maxAllNodePerCell))),
									make_transform_iterator(iBegin2,
										DivideFunctor(maxAllNodePerCell)),
									make_transform_iterator(iBegin2,
										ModuloFunctor(maxAllNodePerCell)),
									nodes->getInfoVecs().nodeIsActive.begin(),
									nodes->getInfoVecs().nodeVelX.begin(),
									nodes->getInfoVecs().nodeVelY.begin(),
									nodes->getInfoVecs().memNodeType1.begin()))
									+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
								nodes->getInfoVecs().nodeVelX.begin(),
							   	nodes->getInfoVecs().nodeVelY.begin(),
							    nodes->getInfoVecs().nodeF_MM_C_X.begin(),   
							    nodes->getInfoVecs().nodeF_MM_C_Y.begin(),
							    nodes->getInfoVecs().nodeContractEnergyT.begin(),
							    nodes->getInfoVecs().basalContractPair.begin())),
			// AddMemContractForce(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,nodeLocYAddr, ActomyosinMultipBasal,ActomyosinMultipApical, nodeTypeAddr,nodeMemMirrorIndexAddr, 
			// 						timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_Division));
			AddMemContractForce_tmp2(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,nodeLocYAddr, nodeTypeAddr,nodeMemMirrorIndexAddr, 
									timeRatio,  contractActomyo_multip, contractActomyo_multip_apical, distFromNucleus_max, distFromNucleus_min, 
									distFromNucleus_normal, distFromNucleus_normal_apical, mitoRndActomyoStrengthScaling,//percentage_before_timeRatio_Crit_Division_scaling, 
									isEnteringMitotic, contractileSpringGrowthProgress, individualCellHeight, distFromBasalLoc, distFromApicalLoc));

	
}


// this function is not currently active. it was useful when nucleus is modeled as one point node and then force interaction between nucleus and other nodes was implemented in this function. The force interaction is one-way so only effect of nucleus on other nodes.
void SceCells::applyForceInteractionNucleusAsPoint() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);

	//double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
	//		- growthAuxData.prolifDecay
	//				* (growthAuxData.grthProgrEndCPU
	//						- growthAuxData.grthPrgrCriVal_M_Ori);

	double grthPrgrCriVal_M = growthAuxData.grthPrgrCriVal_M_Ori; 
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.nucleusLocX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.nucleusLocY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.nucleusLocX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.nucleusLocY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))		
							+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
									   nodes->getInfoVecs().nodeVelY.begin())),
			AddNucleusForce(grthPrgrCriVal_M));
	
}

void SceCells::PlotNucleus (int & lastPrintNucleus, int & outputFrameNucleus) {
	lastPrintNucleus=lastPrintNucleus+1 ; 
    if (lastPrintNucleus>=10000) { 
		outputFrameNucleus++ ; 
		lastPrintNucleus=0 ; 
		std::string vtkFileName = "Nucleus_" + patch::to_string(outputFrameNucleus-1) + ".vtk";
		ofstream NucleusOut;
		NucleusOut.open(vtkFileName.c_str());
		NucleusOut<< "# vtk DataFile Version 3.0" << endl;
		NucleusOut<< "Result for paraview 2d code" << endl;
		NucleusOut << "ASCII" << endl;
		NucleusOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
		NucleusOut << "POINTS " << allocPara_m.currentActiveCellCount << " float" << std::endl;
		for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
			NucleusOut << cellInfoVecs.InternalAvgX[i] << " " << cellInfoVecs.InternalAvgY[i] << " "
			<< 0.0 << std::endl;
		}
		NucleusOut<< std::endl;


		NucleusOut.close(); 
	}

}



// __device__
// void calAndAddIB_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
// 		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M) {
// 	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

// 	double forceValue = 0;
// 	// if (growPro > grthPrgrCriEnd_M) {
// 	if (1 < 0){
// 		if (linkLength < sceIBDiv_M[4]) {
// 			forceValue = -sceIBDiv_M[0] / sceIBDiv_M[2]
// 					* exp(-linkLength / sceIBDiv_M[2])
// 					+ sceIBDiv_M[1] / sceIBDiv_M[3]
// 							* exp(-linkLength / sceIBDiv_M[3]);
// 		}
// 	} else if (2< 0){//(growPro > grthPrgrCriVal_M) {
// 		double percent = (growPro - grthPrgrCriVal_M)
// 				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
// 		double lenLimit = percent * (sceIBDiv_M[4])
// 				+ (1.0 - percent) * sceIB_M[4];
// 		if (linkLength < lenLimit) {
// 			double intnlBPara0 = percent * (sceIBDiv_M[0])
// 					+ (1.0 - percent) * sceIB_M[0];
// 			double intnlBPara1 = percent * (sceIBDiv_M[1])
// 					+ (1.0 - percent) * sceIB_M[1];
// 			double intnlBPara2 = percent * (sceIBDiv_M[2])
// 					+ (1.0 - percent) * sceIB_M[2];
// 			double intnlBPara3 = percent * (sceIBDiv_M[3])
// 					+ (1.0 - percent) * sceIB_M[3];
// 			forceValue = -intnlBPara0 / intnlBPara2
// 					* exp(-linkLength / intnlBPara2)
// 					+ intnlBPara1 / intnlBPara3
// 							* exp(-linkLength / intnlBPara3);
// 		}
// 	} else {
// 		if (linkLength < sceIB_M[4]) {
// 			forceValue = -sceIB_M[0] / sceIB_M[2]
// 					* exp(-linkLength / sceIB_M[2])
// 					+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / sceIB_M[3]);
// 		}
// 	}
// 	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
// 	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
// }
__device__
void calAndAddIB_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M, bool enteringMitotic) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	// if (growPro > grthPrgrCriEnd_M) {
	if (1 < 0){
		if (linkLength < sceIBDiv_M[4]) {
			forceValue = -sceIBDiv_M[0] / sceIBDiv_M[2]
					* exp(-linkLength / sceIBDiv_M[2])
					+ sceIBDiv_M[1] / sceIBDiv_M[3]
							* exp(-linkLength / sceIBDiv_M[3]);
		}
	} else if (2< 0){//(growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIBDiv_M[4])
				+ (1.0 - percent) * sceIB_M[4];
		if (linkLength < lenLimit) {
			double intnlBPara0 = percent * (sceIBDiv_M[0])
					+ (1.0 - percent) * sceIB_M[0];
			double intnlBPara1 = percent * (sceIBDiv_M[1])
					+ (1.0 - percent) * sceIB_M[1];
			double intnlBPara2 = percent * (sceIBDiv_M[2])
					+ (1.0 - percent) * sceIB_M[2];
			double intnlBPara3 = percent * (sceIBDiv_M[3])
					+ (1.0 - percent) * sceIB_M[3];
			forceValue = -intnlBPara0 / intnlBPara2
					* exp(-linkLength / intnlBPara2)
					+ intnlBPara1 / intnlBPara3
							* exp(-linkLength / intnlBPara3);
		}
	} else {
		if (enteringMitotic == true){
			if (linkLength < 1.25*sceIB_M[4]) {
				forceValue = (-sceIB_M[0] / sceIB_M[2]
						* exp(-linkLength / (1.25*sceIB_M[2]))
						+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / (1.25*sceIB_M[3])));
			}
		}
		else{
			if (linkLength < sceIB_M[4]) {
				forceValue = -sceIB_M[0] / sceIB_M[2]
						* exp(-linkLength / sceIB_M[2])
						+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / sceIB_M[3]);
			}
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void CalAndAddIMEnergy(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& IMEnergyT, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double IMEnergy = 0;
	if (1<0){//(growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIBDiv_M[4]) {
			IMEnergy = sceIBDiv_M[0] 
					* exp(-linkLength / sceIBDiv_M[2])
					-sceIBDiv_M[1] 
							* exp(-linkLength / sceIBDiv_M[3]);
		}
	} else if (2<0){//(growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIBDiv_M[4])
				+ (1.0 - percent) * sceIB_M[4];
		if (linkLength < lenLimit) {
			double intnlBPara0 = percent * (sceIBDiv_M[0])
					+ (1.0 - percent) * sceIB_M[0];
			double intnlBPara1 = percent * (sceIBDiv_M[1])
					+ (1.0 - percent) * sceIB_M[1];
			double intnlBPara2 = percent * (sceIBDiv_M[2])
					+ (1.0 - percent) * sceIB_M[2];
			double intnlBPara3 = percent * (sceIBDiv_M[3])
					+ (1.0 - percent) * sceIB_M[3];
			IMEnergy = intnlBPara0 
					* exp(-linkLength / intnlBPara2)
					- intnlBPara1 
							* exp(-linkLength / intnlBPara3);
		}
	} else {
		if (linkLength < sceIB_M[4]) {
			IMEnergy = sceIB_M[0] 
					* exp(-linkLength / sceIB_M[2])
					- sceIB_M[1] * exp(-linkLength / sceIB_M[3]);
		}
	}
	IMEnergyT=IMEnergyT+IMEnergy ; 
}




//Ali function added for eventually computing pressure for each cells
// __device__
// void calAndAddIB_M2(double& xPos, double& yPos, double& xPos2, double& yPos2,
// 		double& growPro, double& xRes, double& yRes, double & F_MI_M_x, double & F_MI_M_y, double grthPrgrCriVal_M) {
// 	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

// 	double forceValue = 0;
// 	if (1<0){//if (growPro > grthPrgrCriEnd_M) {
// 		if (linkLength < sceIBDiv_M[4]) {
// 			forceValue = -sceIBDiv_M[0] / sceIBDiv_M[2]
// 					* exp(-linkLength / sceIBDiv_M[2])
// 					+ sceIBDiv_M[1] / sceIBDiv_M[3]
// 							* exp(-linkLength / sceIBDiv_M[3]);
// 		}
// 	} else if (2<0){//(growPro > grthPrgrCriVal_M) {
// 		double percent = (growPro - grthPrgrCriVal_M)
// 				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
// 		double lenLimit = percent * (sceIBDiv_M[4])
// 				+ (1.0 - percent) * sceIB_M[4];
// 		if (linkLength < lenLimit) {
// 			double intnlBPara0 = percent * (sceIBDiv_M[0])
// 					+ (1.0 - percent) * sceIB_M[0];
// 			double intnlBPara1 = percent * (sceIBDiv_M[1])
// 					+ (1.0 - percent) * sceIB_M[1];
// 			double intnlBPara2 = percent * (sceIBDiv_M[2])
// 					+ (1.0 - percent) * sceIB_M[2];
// 			double intnlBPara3 = percent * (sceIBDiv_M[3])
// 					+ (1.0 - percent) * sceIB_M[3];
// 			forceValue = -intnlBPara0 / intnlBPara2
// 					* exp(-linkLength / intnlBPara2)
// 					+ intnlBPara1 / intnlBPara3
// 							* exp(-linkLength / intnlBPara3);
// 		}
// 	} else {
// 		if (linkLength < sceIB_M[4]) {
// 			forceValue = -sceIB_M[0] / sceIB_M[2]
// 					* exp(-linkLength / sceIB_M[2])
// 					+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / sceIB_M[3]);
// 		}
// 	}

// 	F_MI_M_x=F_MI_M_x+forceValue * (xPos2 - xPos) / linkLength;
// 	F_MI_M_y=F_MI_M_y+forceValue * (yPos2 - yPos) / linkLength;
       
// 	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
// 	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
// }

__device__
void calAndAddIB_M2(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double & F_MI_M_x, double & F_MI_M_y, double grthPrgrCriVal_M, bool enteringMitotic) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	if (1<0){//if (growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIBDiv_M[4]) {
			forceValue = -sceIBDiv_M[0] / sceIBDiv_M[2]
					* exp(-linkLength / sceIBDiv_M[2])
					+ sceIBDiv_M[1] / sceIBDiv_M[3]
							* exp(-linkLength / sceIBDiv_M[3]);
		}
	} else if (2<0){//(growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIBDiv_M[4])
				+ (1.0 - percent) * sceIB_M[4];
		if (linkLength < lenLimit) {
			double intnlBPara0 = percent * (sceIBDiv_M[0])
					+ (1.0 - percent) * sceIB_M[0];
			double intnlBPara1 = percent * (sceIBDiv_M[1])
					+ (1.0 - percent) * sceIB_M[1];
			double intnlBPara2 = percent * (sceIBDiv_M[2])
					+ (1.0 - percent) * sceIB_M[2];
			double intnlBPara3 = percent * (sceIBDiv_M[3])
					+ (1.0 - percent) * sceIB_M[3];
			forceValue = -intnlBPara0 / intnlBPara2
					* exp(-linkLength / intnlBPara2)
					+ intnlBPara1 / intnlBPara3
							* exp(-linkLength / intnlBPara3);
		}
	} else {
		if (enteringMitotic == true){
			if (linkLength < 1.25*sceIB_M[4]) {
				forceValue = (-sceIB_M[0] / sceIB_M[2]
						* exp(-linkLength / (1.25*sceIB_M[2]))
						+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / (1.25*sceIB_M[3])));
			}
		}
		else{
			if (linkLength < sceIB_M[4]) {
				forceValue = -sceIB_M[0] / sceIB_M[2]
						* exp(-linkLength / sceIB_M[2])
						+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / sceIB_M[3]);
			}
		}
	}

	F_MI_M_x=F_MI_M_x+forceValue * (xPos2 - xPos) / linkLength;
	F_MI_M_y=F_MI_M_y+forceValue * (yPos2 - yPos) / linkLength;
       
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddMM_ContractRepl(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes, double & F_MM_C_X, double & F_MM_C_Y) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	double sceMM_C[5] ; 
	for (int i=0 ; i<5 ; i++) {
		sceMM_C[i]=sceIIDiv_M[i] ; 
	}
		
	if (linkLength < sceMM_C[4]) {
		forceValue = -sceMM_C[0] / sceMM_C[2]
					* exp(-linkLength / sceMM_C[2])
					+ sceMM_C[1] / sceMM_C[3] * exp(-linkLength / sceMM_C[3]);
		}
	

	F_MM_C_X=F_MM_C_X+forceValue * (xPos2 - xPos) / linkLength;
	F_MM_C_Y=F_MM_C_Y+forceValue * (yPos2 - yPos) / linkLength;
       
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

// __device__
// void calAndAddMM_ContractAdh(double& xPos, double& yPos, double& xPos2, double& yPos2,
// 		double& xRes, double& yRes, double & F_MM_C_X, double & F_MM_C_Y) {
// 	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

// 	double lZero=0.03125 ;
// 	//double kCAdh=30 ; 
// 	double forceValue = 0;
		
// 	if (linkLength > lZero) {
// 		forceValue =kContractMemb*(linkLength-lZero) ; 
// 		}
	

// 	F_MM_C_X=F_MM_C_X+forceValue * (xPos2 - xPos) / linkLength;
// 	F_MM_C_Y=F_MM_C_Y+forceValue * (yPos2 - yPos) / linkLength;
       
// 	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
// 	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
// }
__device__
void calAndAddMM_ContractAdh(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes, double & F_MM_C_X, double & F_MM_C_Y, double& kContrMemb_multip, double& kContrMemb_multip2) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double lZero=0.03125 ;
	//double kCAdh=30 ; 
	double forceValue = 0;
	double scaling = (kContrMemb_multip + kContrMemb_multip2)/2.0;
		
	if (linkLength > lZero) {
		forceValue =scaling*kContractMemb*(linkLength-lZero) ; 
		}
	

	F_MM_C_X=F_MM_C_X+forceValue * (xPos2 - xPos) / linkLength;
	F_MM_C_Y=F_MM_C_Y+forceValue * (yPos2 - yPos) / linkLength;
       
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}


__device__
void calAndAddII_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	if (1<0){//(growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIIDiv_M[4]) {
			forceValue = -sceIIDiv_M[0] / sceIIDiv_M[2]
					* exp(-linkLength / sceIIDiv_M[2])
					+ sceIIDiv_M[1] / sceIIDiv_M[3]
							* exp(-linkLength / sceIIDiv_M[3]);
		}
	} else if (2<0){//(growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIIDiv_M[4])
				+ (1.0 - percent) * sceII_M[4];
		if (linkLength < lenLimit) {
			double intraPara0 = percent * (sceIIDiv_M[0])
					+ (1.0 - percent) * sceII_M[0];
			double intraPara1 = percent * (sceIIDiv_M[1])
					+ (1.0 - percent) * sceII_M[1];
			double intraPara2 = percent * (sceIIDiv_M[2])
					+ (1.0 - percent) * sceII_M[2];
			double intraPara3 = percent * (sceIIDiv_M[3])
					+ (1.0 - percent) * sceII_M[3];
			forceValue = -intraPara0 / intraPara2
					* exp(-linkLength / intraPara2)
					+ intraPara1 / intraPara3 * exp(-linkLength / intraPara3);
		}
	} 
	else {
		if (linkLength < sceII_M[4]) {
			forceValue = -sceII_M[0] / sceII_M[2]
					* exp(-linkLength / sceII_M[2])
					+ sceII_M[1] / sceII_M[3] * exp(-linkLength / sceII_M[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void CalAndAddIIEnergy(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& IIEnergyT, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double IIEnergy = 0;
	if (1<0){//(growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIIDiv_M[4]) {
			IIEnergy= sceIIDiv_M[0]* exp(-linkLength / sceIIDiv_M[2])
					- sceIIDiv_M[1]* exp(-linkLength / sceIIDiv_M[3]);
		}
	} else if (2<0){//(growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIIDiv_M[4])
				+ (1.0 - percent) * sceII_M[4];
		if (linkLength < lenLimit) {
			double intraPara0 = percent * (sceIIDiv_M[0])
					+ (1.0 - percent) * sceII_M[0];
			double intraPara1 = percent * (sceIIDiv_M[1])
					+ (1.0 - percent) * sceII_M[1];
			double intraPara2 = percent * (sceIIDiv_M[2])
					+ (1.0 - percent) * sceII_M[2];
			double intraPara3 = percent * (sceIIDiv_M[3])
					+ (1.0 - percent) * sceII_M[3];
			IIEnergy = intraPara0  * exp(-linkLength / intraPara2)
					  -intraPara1  * exp(-linkLength / intraPara3);
		}
	} else {
		if (linkLength < sceII_M[4]) {
			IIEnergy = sceII_M[0] * exp(-linkLength / sceII_M[2])
					  -sceII_M[1] * exp(-linkLength / sceII_M[3]);
		}
	}
	IIEnergyT=IIEnergyT+IIEnergy ; 
}



__device__
void calAndAddNucleusEffect(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	if (growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceNDiv_M[4]) {
			forceValue = -sceNDiv_M[0] / sceNDiv_M[2]
					* exp(-linkLength / sceNDiv_M[2])
					+ sceNDiv_M[1] / sceNDiv_M[3]
							* exp(-linkLength / sceNDiv_M[3]);
		}
	} else if (growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceNDiv_M[4])
				+ (1.0 - percent) * sceN_M[4];
		if (linkLength < lenLimit) {
			double intraPara0 = percent * (sceNDiv_M[0])
					+ (1.0 - percent) * sceN_M[0];
			double intraPara1 = percent * (sceNDiv_M[1])
					+ (1.0 - percent) * sceN_M[1];
			double intraPara2 = percent * (sceNDiv_M[2])
					+ (1.0 - percent) * sceN_M[2];
			double intraPara3 = percent * (sceNDiv_M[3])
					+ (1.0 - percent) * sceN_M[3];
			forceValue = -intraPara0 / intraPara2
					* exp(-linkLength / intraPara2)
					+ intraPara1 / intraPara3 * exp(-linkLength / intraPara3);
		}
	} else {
		if (linkLength < sceN_M[4]) {
			forceValue = -sceN_M[0] / sceN_M[2]
					* exp(-linkLength / sceN_M[2])
					+ sceN_M[1] / sceN_M[3] * exp(-linkLength / sceN_M[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	

}



void SceCells::writeNucleusIniLocPercent() {
	
	ofstream output ; 
	thrust::host_vector <double> nucleusLocPercentHost ; 
	
	string uniqueSymbolOutput = globalConfigVars.getConfigValue("UniqueSymbol").toString();
	std::string resumeFileName = "./resources/DataFileInitLocNucleusPercent_" + uniqueSymbolOutput + "Resume.cfg";
	output.open(resumeFileName.c_str() );
	nucleusLocPercentHost=cellInfoVecs.nucleusLocPercent ; 

	for (int i=0 ; i<allocPara_m.currentActiveCellCount  ; i++){
		output << i <<"	"<<nucleusLocPercentHost[i] << endl ; 
	}
	
	output.close() ; 
}


// void SceCells::readNucleusIniLocPercent() {

// 	ifstream input ; 
// 	vector <double> nucleusLocPercentHost ; 
// 	int dummy ; 
// 	double percent ;
// 	string uniqueSymbol = globalConfigVars.getConfigValue("UniqueSymbol").toString();
// 	string resumeFileName = "./resources/DataFileInitLocNucleusPercent_" + uniqueSymbol + "Resume.cfg";
// 	input.open(resumeFileName.c_str() );
	
// 	if (input.is_open()) {
// 		cout << " Suceessfully openend resume input file for initial locations of nucleus" << endl ; 
// 	}
// 	else{
// 		throw std::invalid_argument ("Failed openening the resume input file for initial locations of nucleus")  ; 

// 	}

// 	for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++){
// 		input >> dummy >> percent ;  
// 		nucleusLocPercentHost.push_back(percent) ;  
// 	}

// 	input.close() ;

// 	cellInfoVecs.nucleusLocPercent= nucleusLocPercentHost ; 
// } // Original function dealing with nucleus position percentage if there is no division at all.

void SceCells::readNucleusIniLocPercent() {

	ifstream input ; 
	vector <double> nucleusLocPercentHost ; 
	int dummy ; 
	double percent ;
	string uniqueSymbol = globalConfigVars.getConfigValue("UniqueSymbol").toString();
	string resumeFileName = "./resources/DataFileInitLocNucleusPercent_" + uniqueSymbol + "Resume.cfg";
	input.open(resumeFileName.c_str() );
	
	if (input.is_open()) {
		cout << " Suceessfully openend resume input file for initial locations of nucleus" << endl ; 
	}
	else{
		throw std::invalid_argument ("Failed openening the resume input file for initial locations of nucleus")  ; 

	}

	for (int i=0 ; i<allocPara_m.currentActiveCellCount ; i++){
		input >> dummy >> percent ;  
		nucleusLocPercentHost.push_back(percent) ;  
	}
	if (allocPara_m.currentActiveCellCount < allocPara_m.maxCellCount){
		std::cout<<"The number of currently active cell is less than the maximally allowed number, we will reserve space for additional cells."<<std::endl;
		std::cout<<"The nucleus percentage of inactive cells will be set to 0.0 initially, and must be updated when new cells are introduced."<<std::endl;
		int gap = allocPara_m.maxCellCount - allocPara_m.currentActiveCellCount;
		for (int gap_count = 0; gap_count < gap; gap_count++){
			nucleusLocPercentHost.push_back(0.0);
		}

	}

	input.close() ;

	cellInfoVecs.nucleusLocPercent= nucleusLocPercentHost ; 
}





void SceCells::allComponentsMoveImplicitPart() 
{
  vector <int> indexPrev, indexNext ;
#ifdef debugModeECM 
	cudaEvent_t start1, start2, start3,  stop;
	float elapsedTime1, elapsedTime2, elapsedTime3  ; 
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&start3);
	cudaEventCreate(&stop);
	cudaEventRecord(start1, 0);
#endif


  CalRHS();
#ifdef debugModeECM
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif

  EquMotionCoef(indexPrev, indexNext);

#ifdef debugModeECM
	cudaEventRecord(start3, 0);
	cudaEventSynchronize(start3);
	cudaEventElapsedTime(&elapsedTime2, start2, start3);
#endif


  UpdateLocations(indexPrev,indexNext); 
#ifdef debugModeECM
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime3, start3, stop);
	std::cout << "time 1 spent in cell-solver module for moving the membrane node of cells and ECM nodes are: " << elapsedTime1 << endl ; 
	std::cout << "time 2 spent in cell-solver module for moving the membrane node of cells and ECM nodes are: " << elapsedTime2 << endl ; 
	std::cout << "time 3 spent in cell-solver module for moving the membrane node of cells and ECM nodes are: " << elapsedTime3 << endl ; 
#endif


}

void SceCells::StoreNodeOldPositions() {
	//nodes->getInfoVecs().locXOldHost.clear();
	//nodes->getInfoVecs().locYOldHost.clear(); 
	//nodes->getInfoVecs().locXOldHost.resize(totalNodeCountForActiveCells) ; 
	//nodes->getInfoVecs().locYOldHost.resize(totalNodeCountForActiveCells) ;

	thrust::copy (nodes->getInfoVecs().nodeLocX.begin(), nodes->getInfoVecs().nodeLocX.begin() +
				totalNodeCountForActiveCells, nodes->getInfoVecs().locXOldHost.begin()); 
	thrust::copy (nodes->getInfoVecs().nodeLocY.begin() , nodes->getInfoVecs().nodeLocY.begin() +
				totalNodeCountForActiveCells, nodes->getInfoVecs().locYOldHost.begin()); 

}
void SceCells::CalRHS () {
 //   cout << "total node count for active cells in CalRHS function is="<<totalNodeCountForActiveCells << endl ; 

	//nodes->getInfoVecs().rHSXHost.clear() ; 
	//nodes->getInfoVecs().rHSYHost.clear() ; 
	//nodes->getInfoVecs().rHSXHost.resize(totalNodeCountForActiveCells) ; 
	//nodes->getInfoVecs().rHSYHost.resize(totalNodeCountForActiveCells) ; 

	thrust::copy( 
	      thrust::make_zip_iterator(
		       thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(), 
	                              nodes->getInfoVecs().nodeLocY.begin())),
		  thrust::make_zip_iterator(
		       thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(), 
	                              nodes->getInfoVecs().nodeLocY.begin()))+totalNodeCountForActiveCells,
		  thrust::make_zip_iterator(
		   	   thrust::make_tuple(nodes->getInfoVecs().rHSXHost.begin(),
			                      nodes->getInfoVecs().rHSYHost.begin()))); 
}

void SceCells::EquMotionCoef(vector<int> & indexPrev, vector<int> & indexNext)
{
   vector <uint> activeMemCount(allocPara_m.currentActiveCellCount) ;
   double distWithNext[totalNodeCountForActiveCells]  ; 
   double distWithPrev[totalNodeCountForActiveCells]  ;
   int cellRank ; 
   int nodeRank ;

   indexPrev.clear() ; 
   indexNext.clear() ;
   //nodes->getInfoVecs().hCoefD.clear() ; 
   //nodes->getInfoVecs().hCoefLd.clear() ; 
   //nodes->getInfoVecs().hCoefUd.clear() ;
   //nodes->getInfoVecs().nodeIsActiveH.clear(); 
   
   indexPrev.resize(totalNodeCountForActiveCells) ; 
   indexNext.resize(totalNodeCountForActiveCells) ; 
   //nodes->getInfoVecs().hCoefD.resize(totalNodeCountForActiveCells,0.0) ; 
   //nodes->getInfoVecs().hCoefLd.resize(totalNodeCountForActiveCells,0.0) ; 
   //nodes->getInfoVecs().hCoefUd.resize(totalNodeCountForActiveCells,0.0) ;
   //nodes->getInfoVecs().nodeIsActiveH.resize(totalNodeCountForActiveCells) ; 
   
   thrust::copy (nodes->getInfoVecs().nodeIsActive.begin(),
        	     nodes->getInfoVecs().nodeIsActive.begin()+ totalNodeCountForActiveCells,
		         nodes->getInfoVecs().nodeIsActiveH.begin()); 
	
   thrust::copy(cellInfoVecs.activeMembrNodeCounts.begin() , cellInfoVecs.activeMembrNodeCounts.begin()+ 
         allocPara_m.currentActiveCellCount, activeMemCount.begin()); 

   
	//cout << "Maximum all node per cells is " << allocPara_m.maxAllNodePerCell << endl ;  
   for ( int i=0 ;  i< totalNodeCountForActiveCells ; i++) {
	   cellRank=i/allocPara_m.maxAllNodePerCell ; 
	   nodeRank=i%allocPara_m.maxAllNodePerCell ;
	   if ( nodeRank<activeMemCount [cellRank]) {
	      indexNext.at(i)=i+1 ;
	   	  indexPrev.at(i)=i-1 ;
	   	  if ( nodeRank==activeMemCount [cellRank]-1){
	         indexNext.at(i)=cellRank*allocPara_m.maxAllNodePerCell ;
		  //	cout << "index next for cell rank " << cellRank << " is " << indexNext.at(i) << endl ; 
	      }
	      if (nodeRank==0){
	         indexPrev.at(i)=cellRank*allocPara_m.maxAllNodePerCell  +activeMemCount [cellRank]-1  ; 
          //   cout << "Active membrane nodes for cell rank " << cellRank << " is " <<activeMemCount [cellRank]<<endl ;  
		  //   cout << "index previous for cell rank " << cellRank << " is " << indexPrev.at(i) << endl ; 
	      }
	      distWithNext[i]=sqrt( pow(nodes->getInfoVecs().locXOldHost[indexNext.at(i)] - 
		                            nodes->getInfoVecs().locXOldHost[i],2) + 
	                            pow(nodes->getInfoVecs().locYOldHost[indexNext.at(i)] - 
								    nodes->getInfoVecs().locYOldHost[i],2)) ;
	      distWithPrev[i]=sqrt( pow(nodes->getInfoVecs().locXOldHost[indexPrev.at(i)] - 
		                            nodes->getInfoVecs().locXOldHost[i],2) + 
	                            pow(nodes->getInfoVecs().locYOldHost[indexPrev.at(i)] - 
								    nodes->getInfoVecs().locYOldHost[i],2)); 

   	   }
   }

   double sponLen= globalConfigVars.getConfigValue("MembrEquLen").toDouble();
   double k= globalConfigVars.getConfigValue("MembrStiff").toDouble();
   for ( int i=0 ;  i< totalNodeCountForActiveCells ; i++) {

      if (nodes->getInfoVecs().nodeIsActiveH.at(i)==false) {
		  continue ; 
	  }
	  cellRank=i / allocPara_m.maxAllNodePerCell; 
	  nodeRank=i % allocPara_m.maxAllNodePerCell;

	  if (nodeRank<activeMemCount [cellRank]) {
      	nodes->getInfoVecs().hCoefD[i]= 1 + k*dt/Damp_Coef*( 2 - sponLen/(distWithPrev[i]+0.2*sponLen) - sponLen/(distWithNext[i]+0.2*sponLen)) ; 
	  	nodes->getInfoVecs().hCoefLd[i]=    k*dt/Damp_Coef*(-1 + sponLen/(distWithPrev[i]+0.2*sponLen)) ; 
	  	nodes->getInfoVecs().hCoefUd[i]=    k*dt/Damp_Coef*(-1 + sponLen/(distWithNext[i]+0.2*sponLen)) ; 
   	  }
	  else { // no spring between neighboring points exist 
      	nodes->getInfoVecs().hCoefD[i]=1.0 ; 
	  	nodes->getInfoVecs().hCoefLd[i]=0.0 ; 
	  	nodes->getInfoVecs().hCoefUd[i]=0.0 ; 
	  }

  }

}


void SceCells::UpdateLocations(const vector <int> & indexPrev,const vector <int> & indexNext ) {
   vector <double> locXTmpHost=solverPointer->SOR3DiagPeriodic(nodes->getInfoVecs().nodeIsActiveH,
   												     		     nodes->getInfoVecs().hCoefLd, 
													             nodes->getInfoVecs().hCoefD, 
													             nodes->getInfoVecs().hCoefUd,
													             nodes->getInfoVecs().rHSXHost,
																 indexPrev,indexNext,
													             nodes->getInfoVecs().locXOldHost); 
    
   vector <double> locYTmpHost=solverPointer->SOR3DiagPeriodic(nodes->getInfoVecs().nodeIsActiveH,
												     		     nodes->getInfoVecs().hCoefLd, 
													             nodes->getInfoVecs().hCoefD, 
													             nodes->getInfoVecs().hCoefUd,
													             nodes->getInfoVecs().rHSYHost,
																 indexPrev,indexNext,
													             nodes->getInfoVecs().locYOldHost); 
   
   thrust::copy(
            thrust::make_zip_iterator(
                thrust::make_tuple (locXTmpHost.begin(), 
                                    locYTmpHost.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple (locXTmpHost.begin(), 
                                    locYTmpHost.begin()))+totalNodeCountForActiveCells,
		   thrust::make_zip_iterator(
			    thrust::make_tuple (nodes->getInfoVecs().nodeLocX.begin(),
							        nodes->getInfoVecs().nodeLocY.begin()))); 


}

