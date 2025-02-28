//Notes: 1-  In infoVecs.nodeCellRankBehind and  infoVecs.nodeCellRankFront are given sequential values correspond to their actual values in function: SceNodes::allocSpaceForNodesdh
// 2- the algorithm of adhesion won't work if there is no  apical node.
// 3- maxNumAdh is given inside the code as a parameters in .cu file. It should become an input or I should write a function to detect that automatically
// 4- In SceNodes::NumAdhAfter and SceNodes::NumAdhBefore  number of lateral nodes is given manually inside the code it should be automatically calculate from the input parameters/

#include "SceNodes.h"
#include "SceCells.h"

__constant__ double sceInterPara[5];
__constant__ double sceIntraPara[5];
// parameter set for cells that are going to divide
__constant__ double sceIntraParaDiv[5];
__constant__ double sceDivProPara;
__constant__ double sceCartPara[5];
__constant__ double sceInterDiffPara[5];
__constant__ double sceProfilePara[7];
__constant__ double sceECMPara[5];
__constant__ double sceDiffPara[5];

__constant__ double cartGrowDirVec[3];
__constant__ uint ProfilebeginPos;
__constant__ uint ECMbeginPos;
__constant__ uint cellNodeBeginPos;
__constant__ uint nodeCountPerECM;
__constant__ uint nodeCountPerCell;
//
__constant__ uint cellNodeBeginPos_M;
__constant__ uint allNodeCountPerCell_M;
__constant__ uint membrThreshold_M;
__constant__ double sceInterBPara_M[5];
__constant__ int    sceInterBPara_Jones_On_M ;  //Ali
__constant__ double sceInterBPara_Jones_M[3] ;  //Ali
__constant__ double sceIntnlBPara_M[5];
__constant__ double sceIntraPara_M[5];
__constant__ double sceIntraParaDiv_M[5];
__constant__ double growthPrgrCriVal_M;
__constant__ double maxAdhBondLen_M;
__constant__ double minAdhBondLen_M;
__constant__ double bondStiff_M;
__constant__ double bondStiff_Mitotic;
__constant__ double bondAdhCriLen_M;

// #define DebugMode

// This template method expands an input sequence by
// replicating each element a variable number of times. For example,
//
//   expand([2,2,2],[A,B,C]) -> [A,A,B,B,C,C]
//   expand([3,0,1],[A,B,C]) -> [A,A,A,C]
//   expand([1,3,2],[A,B,C]) -> [A,B,B,B,C,C]
//
// The element counts are assumed to be non-negative integers
template<typename InputIterator1, typename InputIterator2,
		typename OutputIterator>
OutputIterator expand(InputIterator1 first1, InputIterator1 last1,
		InputIterator2 first2, OutputIterator output) {
	typedef typename thrust::iterator_difference<InputIterator1>::type difference_type;

	difference_type input_size = thrust::distance(first1, last1);
	difference_type output_size = thrust::reduce(first1, last1);

	// scan the counts to obtain output offsets for each input element
	thrust::device_vector<difference_type> output_offsets(input_size, 0);
	thrust::exclusive_scan(first1, last1, output_offsets.begin());

	// scatter the nonzero counts into their corresponding output positions
	thrust::device_vector<difference_type> output_indices(output_size, 0);
	thrust::scatter_if(thrust::counting_iterator<difference_type>(0),
			thrust::counting_iterator<difference_type>(input_size),
			output_offsets.begin(), first1, output_indices.begin());

	// compute max-scan over the output indices, filling in the holes
	thrust::inclusive_scan(output_indices.begin(), output_indices.end(),
			output_indices.begin(), thrust::maximum<difference_type>());

	// gather input values according to index array (output = first2[output_indices])
	OutputIterator output_end = output;
	thrust::advance(output_end, output_size);
	thrust::gather(output_indices.begin(), output_indices.end(), first2,
			output);

	// return output + output_size
	thrust::advance(output, output_size);
	return output;
}

// Ali. This constructor will be called in Disc_M simulation, but it will be overwritten by another constructor
SceNodes::SceNodes() {
	readDomainPara();
}

int SceNodes::NumAdhBefore(int cellRank,ECellType eCellType) {
	if (eCellType==peri) {
		return 28 ; 
	}
	if (eCellType==bc) {
		return 30 ; 
	}
	if (eCellType==pouch) {
		if (cellRank==0){
			return 36 ; 
		}
			else if ( cellRank==1) {
				return 108 ; 
			}
				else if ( cellRank==64) {
					return 108 ; 					
				}
					else {
						return 180 ; 
					}
	}

}

int SceNodes::NumAdhAfter(int cellRank,ECellType eCellType) {
	if ( eCellType==peri) {
		return 28 ; 
	}
	if (eCellType==bc) {
		return 30 ; 
	}
	if (eCellType==pouch) {
		if (cellRank==64){
			return 36 ; 
		}
			else if ( cellRank==63) {
				return 108 ; 
			}
				else if ( cellRank==0) {
					return 108 ; 
				}
					else {
						return 180 ; 
					}
	}
}




void SceNodes::readDomainPara() {
	domainPara.minX = globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble();
	domainPara.maxX = globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble();
	domainPara.minY = globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble();
	domainPara.maxY = globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble();
	//domainPara.minZ = globalConfigVars.getConfigValue("DOMAIN_ZMIN").toDouble();
	//domainPara.maxZ = globalConfigVars.getConfigValue("DOMAIN_ZMAX").toDouble();
	domainPara.gridSpacing = getMaxEffectiveRange();
	domainPara.XBucketSize = (domainPara.maxX - domainPara.minX)
			/ domainPara.gridSpacing + 1;
	domainPara.YBucketSize = (domainPara.maxY - domainPara.minY)
			/ domainPara.gridSpacing + 1;
	//domainPara.ZBucketSize = (domainPara.maxZ - domainPara.minZ)
	//		/ domainPara.gridSpacing + 1;
}

void SceNodes::readMechPara() {
	double U0 =
			globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_U0_DivFactor").toDouble();
	double V0 =
			globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_V0_DivFactor").toDouble();
	double k1 =
			globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_k1_DivFactor").toDouble();
	double k2 =
			globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_k2_DivFactor").toDouble();

	mechPara.sceInterParaCPU[0] = U0;
	mechPara.sceInterParaCPU[1] = V0;
	mechPara.sceInterParaCPU[2] = k1;
	mechPara.sceInterParaCPU[3] = k2;

	double interLinkEffectiveRange;
	if (controlPara.simuType != Disc_M) {
		interLinkEffectiveRange = globalConfigVars.getConfigValue(
				"InterCellLinkEffectRange").toDouble();
		mechPara.sceInterParaCPU[4] = interLinkEffectiveRange;
	}

	double U0_Intra =
			globalConfigVars.getConfigValue("IntraCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_U0_DivFactor").toDouble();
	double V0_Intra =
			globalConfigVars.getConfigValue("IntraCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_V0_DivFactor").toDouble();
	double k1_Intra =
			globalConfigVars.getConfigValue("IntraCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_k1_DivFactor").toDouble();
	double k2_Intra =
			globalConfigVars.getConfigValue("IntraCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_k2_DivFactor").toDouble();

	mechPara.sceIntraParaCPU[0] = U0_Intra;
	mechPara.sceIntraParaCPU[1] = V0_Intra;
	mechPara.sceIntraParaCPU[2] = k1_Intra;
	mechPara.sceIntraParaCPU[3] = k2_Intra;

	double intraLinkEffectiveRange;
	if (controlPara.simuType != Disc_M) {
		intraLinkEffectiveRange = globalConfigVars.getConfigValue(
				"IntraCellLinkEffectRange").toDouble();
		mechPara.sceIntraParaCPU[4] = intraLinkEffectiveRange;
	}
	if (controlPara.simuType == Disc) {
		double U0_Intra_Div =
				globalConfigVars.getConfigValue("IntraCell_U0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_U0_Div_DivFactor").toDouble();
		double V0_Intra_Div =
				globalConfigVars.getConfigValue("IntraCell_V0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_V0_Div_DivFactor").toDouble();
		double k1_Intra_Div =
				globalConfigVars.getConfigValue("IntraCell_k1_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_k1_Div_DivFactor").toDouble();
		double k2_Intra_Div =
				globalConfigVars.getConfigValue("IntraCell_k2_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_k2_Div_DivFactor").toDouble();
		double growthProgressThreshold = globalConfigVars.getConfigValue(
				"GrowthProgressThreshold").toDouble();

		mechPara.sceIntraParaDivCPU[0] = U0_Intra_Div;
		mechPara.sceIntraParaDivCPU[1] = V0_Intra_Div;
		mechPara.sceIntraParaDivCPU[2] = k1_Intra_Div;
		mechPara.sceIntraParaDivCPU[3] = k2_Intra_Div;
		mechPara.sceIntraParaDivCPU[4] = growthProgressThreshold;
	}
}

// This constructor is not active Ali 
SceNodes::SceNodes(uint totalBdryNodeCount, uint maxProfileNodeCount,
		uint maxCartNodeCount, uint maxTotalECMCount, uint maxNodeInECM,
		uint maxTotalCellCount, uint maxNodeInCell, bool isStab) {

	cout<< " I am inside SceNodes constructor which I beleive is not active" << endl ; 
	initControlPara(isStab);
	readDomainPara();
	uint maxTotalNodeCount;
	if (controlPara.simuType != Disc_M) {
		initNodeAllocPara(totalBdryNodeCount, maxProfileNodeCount,
				maxCartNodeCount, maxTotalECMCount, maxNodeInECM,
				maxTotalCellCount, maxNodeInCell);
		maxTotalNodeCount = totalBdryNodeCount + maxProfileNodeCount
				+ maxCartNodeCount + allocPara.maxTotalECMNodeCount
				+ allocPara.maxTotalCellNodeCount;
	} else {
		uint maxEpiNodeCount = globalConfigVars.getConfigValue(
				"MaxEpiNodeCountPerCell").toInt();
		uint maxInternalNodeCount = globalConfigVars.getConfigValue(
				"MaxAllNodeCountPerCell").toInt() - maxEpiNodeCount;

		initNodeAllocPara_M(totalBdryNodeCount, maxTotalCellCount,
				maxEpiNodeCount, maxInternalNodeCount);
		maxTotalNodeCount = allocPara_M.maxTotalNodeCount;
	}
	//allocSpaceForNodes(maxTotalNodeCount); Ali comment this becasue it is not active in this simulation and I updated the function in the active constructor
	thrust::host_vector<SceNodeType> hostTmpVector(maxTotalNodeCount);
	thrust::host_vector<bool> hostTmpVector2(maxTotalNodeCount);
	thrust::host_vector<int> hostTmpVector3(maxTotalNodeCount);

	if (controlPara.simuType != Disc_M) {

		for (int i = 0; i < maxTotalNodeCount; i++) {
			if (i < allocPara.startPosProfile) {
				hostTmpVector[i] = Boundary;
				hostTmpVector3[i] = 0;
			} else if (i < allocPara.startPosCart) {
				hostTmpVector[i] = Profile;
				hostTmpVector3[i] = 0;
			} else if (i < allocPara.startPosECM) {
				hostTmpVector[i] = Cart;
				hostTmpVector3[i] = 0;
			} else if (i < allocPara.startPosCells) {
				hostTmpVector[i] = ECM;
				hostTmpVector3[i] = (i - allocPara.startPosECM)
						/ allocPara.maxNodePerECM;
			} else {
				// all initialized as FNM
				hostTmpVector[i] = FNM;
				hostTmpVector3[i] = (i - allocPara.startPosCells)
						/ allocPara.maxNodeOfOneCell;
			}
			hostTmpVector2[i] = false;
		}

	} else {
		for (uint i = 0; i < maxTotalNodeCount; i++) {
			if (i < allocPara_M.bdryNodeCount) {
				hostTmpVector[i] = Boundary;
				hostTmpVector3[i] = 0;
			} else {
				uint tmp = i - allocPara_M.bdryNodeCount;
				uint cellRank = tmp / allocPara_M.bdryNodeCount;
				uint nodeRank = tmp % allocPara_M.bdryNodeCount;
				if (nodeRank < allocPara_M.maxMembrNodePerCell) {
					hostTmpVector[i] = CellMembr;
				} else {
					hostTmpVector[i] = CellIntnl;
				}
				hostTmpVector3[i] = cellRank;
			}
			hostTmpVector2[i] = false;
		}
	}
	infoVecs.nodeCellType = hostTmpVector;
	infoVecs.nodeIsActive = hostTmpVector2;
	infoVecs.nodeCellRank = hostTmpVector3;
       std::cout << " I am in SceNodes constructor with long input which includes copyParaToGPUConstMem function " << endl ; 
	copyParaToGPUConstMem();
}

SceNodes::SceNodes(uint maxTotalCellCount, uint maxAllNodePerCell, uint currentActiveCellCount) {
	//initControlPara (isStab);


	cout<< " I am inside SceNodes constructor which I beleive is active" << endl ;

	isMemNodeTypeAssigned=false ; 
	isMemNodeTypeAssigned_postCellDivision=false;
	isMemNodeTypeAssigned_postAddNode=true;
	isActinLevelDisplayed=false;
	isCellAreaDisplayed= false;
	isECMPairPrinted=false;
	isBucketValuesIncludingNeighborPrinted=false;
	isApicalAdhPresent=true ;
	int simuTypeConfigValue =
			globalConfigVars.getConfigValue("SimulationType").toInt();
	controlPara.simuType = parseTypeFromConfig(simuTypeConfigValue);
	controlPara.resumeSimulation = globalConfigVars.getConfigValue(
			"ResumeSimulation").toInt();

	readDomainPara();
	uint maxTotalNodeCount = maxTotalCellCount * maxAllNodePerCell;

	uint maxMembrNodeCountPerCell = globalConfigVars.getConfigValue(
			"MaxMembrNodeCountPerCell").toInt();
	uint maxIntnlNodeCountPerCell = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();

	initNodeAllocPara_M(0, maxTotalCellCount, maxMembrNodeCountPerCell,
			maxIntnlNodeCountPerCell);

	std::cout << "    Number of boundary nodes = " << allocPara_M.bdryNodeCount
			<< std::endl;
	std::cout << "    Max number of cells in domain = "
			<< allocPara_M.maxCellCount << std::endl;
	std::cout << "    Max all nodes per cell = "
			<< allocPara_M.maxAllNodePerCell << std::endl;
	std::cout << "    Max membrane node per cell= "
			<< allocPara_M.maxMembrNodePerCell << std::endl;
	std::cout << "    Max internal node per cell= "
			<< allocPara_M.maxIntnlNodePerCell << std::endl;
	std::cout << "    Max total number of nodes in domain = "
			<< allocPara_M.maxTotalNodeCount << std::endl;

	allocSpaceForNodes(maxTotalNodeCount, allocPara_M.maxCellCount, currentActiveCellCount);
	std::cout<<"Finished allocSpaceForNodes"<<std::endl;
	thrust::host_vector<SceNodeType> hostTmpVector(maxTotalNodeCount);
	thrust::host_vector<bool> hostTmpVector2(maxTotalNodeCount);

	uint nodeRank;
	for (uint i = 0; i < maxTotalNodeCount; i++) {
		if (i < allocPara_M.bdryNodeCount) {
			hostTmpVector[i] = Boundary;
		} else {
			uint tmp = i - allocPara_M.bdryNodeCount;
			nodeRank = tmp % allocPara_M.maxAllNodePerCell;
			if (nodeRank < allocPara_M.maxMembrNodePerCell) {
				hostTmpVector[i] = CellMembr;
				//std::cout << "0";
			} else {
				hostTmpVector[i] = CellIntnl;
				//std::cout << "1";
			}

		}
		hostTmpVector2[i] = false;
		if (nodeRank == 0) {
			//std::cout << std::endl;
		}
	}
	std::cout << "finished" << std::endl;
	//std::cout.flush();
	infoVecs.nodeCellType = hostTmpVector;
	infoVecs.nodeIsActive = hostTmpVector2;

	thrust::host_vector<int> bondVec(maxTotalNodeCount, -1);
	infoVecs.nodeAdhereIndex = bondVec;
	infoVecs.membrIntnlIndex = bondVec;
	infoVecs.nodeAdhIndxHostCopy = bondVec;
	std::cout << "copy finished!" << std::endl;
	//std::cout.flush();
	copyParaToGPUConstMem_M();
        std::cout << " I am in SceNodes constructor with short input which includes copyParaToGPUConstMem_M  function " << endl ; 
	//std::cout << "at the end" << std::endl;
	//std::cout.flush();
	adhUpdate=true ; //Ali 
	cout << "adhesion update is initialized as " << adhUpdate  << endl ; 
}

void SceNodes::copyParaToGPUConstMem() {
	readMechPara();

	cudaMemcpyToSymbol(sceInterPara, mechPara.sceInterParaCPU,
			5 * sizeof(double));

	cudaMemcpyToSymbol(sceIntraPara, mechPara.sceIntraParaCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceIntraParaDiv, mechPara.sceIntraParaDivCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(ProfilebeginPos, &allocPara.startPosProfile,
			sizeof(uint));
	cudaMemcpyToSymbol(ECMbeginPos, &allocPara.startPosECM, sizeof(uint));
	cudaMemcpyToSymbol(cellNodeBeginPos, &allocPara.startPosCells,
			sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerECM, &allocPara.maxNodePerECM, sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerCell, &allocPara.maxNodeOfOneCell,
			sizeof(uint));
	cudaMemcpyToSymbol(sceCartPara, mechPara.sceCartParaCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceProfilePara, mechPara.sceProfileParaCPU,
			7 * sizeof(double));
	cudaMemcpyToSymbol(sceInterDiffPara, mechPara.sceInterDiffParaCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceECMPara, mechPara.sceECMParaCPU, 5 * sizeof(double));
}

void SceNodes::copyParaToGPUConstMem_M() {
	readParas_M();
	cudaMemcpyToSymbol(cellNodeBeginPos_M, &allocPara_M.bdryNodeCount,
			sizeof(uint));
	cudaMemcpyToSymbol(allNodeCountPerCell_M, &allocPara_M.maxAllNodePerCell,
			sizeof(uint));
	cudaMemcpyToSymbol(membrThreshold_M, &allocPara_M.maxMembrNodePerCell,
			sizeof(uint));
	cudaMemcpyToSymbol(bondAdhCriLen_M, &mechPara_M.bondAdhCriLenCPU_M,
			sizeof(double));

	cudaMemcpyToSymbol(bondStiff_M, &mechPara_M.bondStiffCPU_M, sizeof(double));
	cudaMemcpyToSymbol(bondStiff_Mitotic, &mechPara_M.bondStiffCPU_Mitotic, sizeof(double));//Ali June 16
	cudaMemcpyToSymbol(growthPrgrCriVal_M, &mechPara_M.growthPrgrCriValCPU_M,
			sizeof(double));
	cudaMemcpyToSymbol(maxAdhBondLen_M, &mechPara_M.maxAdhBondLenCPU_M,
			sizeof(double));
	cudaMemcpyToSymbol(minAdhBondLen_M, &mechPara_M.minAdhBondLenCPU_M,
			sizeof(double));
	cudaMemcpyToSymbol(sceInterBPara_M, mechPara_M.sceInterBParaCPU_M,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceInterBPara_Jones_On_M, &mechPara_M.sceInterBParaCPU_Jones_On_M,
			    sizeof(int)); //Ali 
	cudaMemcpyToSymbol(sceInterBPara_Jones_M,    mechPara_M.sceInterBParaCPU_Jones_M,
			3 * sizeof(double)); //Ali 
	cudaMemcpyToSymbol(sceIntnlBPara_M, mechPara_M.sceIntnlBParaCPU_M,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceIntraPara_M, mechPara_M.sceIntraParaCPU_M,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceIntraParaDiv_M, mechPara_M.sceIntraParaDivCPU_M,
			5 * sizeof(double));

}

void SceNodes::initDimension(double domainMinX, double domainMaxX,
		double domainMinY, double domainMaxY, double domainBucketSize) {
	domainPara.minX = domainMinX;
	domainPara.maxX = domainMaxX;
	domainPara.minY = domainMinY;
	domainPara.maxY = domainMaxY;
	domainPara.gridSpacing = domainBucketSize;
	domainPara.XBucketSize = (domainPara.maxX - domainPara.minX)
			/ domainPara.gridSpacing + 1;
	domainPara.YBucketSize = (domainPara.maxY - domainPara.minY)
			/ domainPara.gridSpacing + 1;
	domainPara.totalBucketCount = domainPara.XBucketSize
			* domainPara.YBucketSize;

	auxVecs.keyBegin.resize(domainPara.totalBucketCount);
	auxVecs.keyEnd.resize(domainPara.totalBucketCount);
}

std::vector<std::pair<uint, uint> > SceNodes::obtainPossibleNeighborPairs() {
	std::vector<std::pair<uint, uint> > result;
	thrust::host_vector<uint> keyBeginCPU = auxVecs.keyBegin;
	thrust::host_vector<uint> keyEndCPU = auxVecs.keyEnd;
	thrust::host_vector<uint> bucketKeysCPU = auxVecs.bucketKeys;
	thrust::host_vector<uint> bucketValuesCPU = auxVecs.bucketValues;
	thrust::host_vector<uint> bucketValuesExtendedCPU =
			auxVecs.bucketValuesIncludingNeighbor;
	uint iterationCounter = 0;
	int size = bucketKeysCPU.size();
	for (int i = 0; i < size; i++) {
		for (int j = keyBeginCPU[bucketKeysCPU[i]];
				j < keyEndCPU[bucketKeysCPU[i]]; j++) {
			int node1 = bucketValuesCPU[i];
			int node2 = bucketValuesExtendedCPU[j];
			if (node1 >= node2) {
				continue;
			} else {
				result.push_back(std::make_pair<uint, uint>(node1, node2));
			}
			iterationCounter++;
		}
	}
	return result;
}

void SceNodes::readParas_M() {


	//////////////////////
	//// Block 1 /////////
	//////////////////////
	double U0_InterB =
			globalConfigVars.getConfigValue("SceInterB_U0").toDouble();
	double V0_InterB =
			globalConfigVars.getConfigValue("SceInterB_V0").toDouble();
	double k1_InterB =
			globalConfigVars.getConfigValue("SceInterB_k1").toDouble();
	double k2_InterB =
			globalConfigVars.getConfigValue("SceInterB_k2").toDouble();
	double interBEffectiveRange = globalConfigVars.getConfigValue(
			"InterBEffectiveRange").toDouble();
	mechPara_M.sceInterBParaCPU_M[0] = U0_InterB;
	mechPara_M.sceInterBParaCPU_M[1] = V0_InterB;
	mechPara_M.sceInterBParaCPU_M[2] = k1_InterB;
	mechPara_M.sceInterBParaCPU_M[3] = k2_InterB;
	mechPara_M.sceInterBParaCPU_M[4] = interBEffectiveRange;

//Ali
	//////////////////////
	//// Block 1.5 /////////
	//////////////////////
	int     On_InterB_Jones =
			globalConfigVars.getConfigValue("SceInterB_Jones_On").toDouble();
	double eps_InterB_Jones =
			globalConfigVars.getConfigValue("SceInterB_Jones_eps").toDouble();
	double sig_InterB_Jones =
			globalConfigVars.getConfigValue("SceInterB_Jones_sig").toDouble();
	double interBEffectiveRange_Jones = globalConfigVars.getConfigValue(
			"InterBEffectiveRange_Jones").toDouble();
        
	mechPara_M.sceInterBParaCPU_Jones_On_M = On_InterB_Jones;
	mechPara_M.sceInterBParaCPU_Jones_M[0] = eps_InterB_Jones;
	mechPara_M.sceInterBParaCPU_Jones_M[1] = sig_InterB_Jones;
	mechPara_M.sceInterBParaCPU_Jones_M[2] = interBEffectiveRange_Jones;
//Ali



	//////////////////////
	//// Block 2 /////////
	//////////////////////
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
	mechPara_M.sceIntnlBParaCPU_M[0] = U0_IntnlB;
	mechPara_M.sceIntnlBParaCPU_M[1] = V0_IntnlB;
	mechPara_M.sceIntnlBParaCPU_M[2] = k1_IntnlB;
	mechPara_M.sceIntnlBParaCPU_M[3] = k2_IntnlB;
	mechPara_M.sceIntnlBParaCPU_M[4] = intnlBEffectiveRange;

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
	mechPara_M.sceIntraParaCPU_M[0] = U0_Intra;
	mechPara_M.sceIntraParaCPU_M[1] = V0_Intra;
	mechPara_M.sceIntraParaCPU_M[2] = k1_Intra;
	mechPara_M.sceIntraParaCPU_M[3] = k2_Intra;
	mechPara_M.sceIntraParaCPU_M[4] = intraLinkEffectiveRange;

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
	mechPara_M.sceIntraParaDivCPU_M[0] = U0_Intra_Div;
	mechPara_M.sceIntraParaDivCPU_M[1] = V0_Intra_Div;
	mechPara_M.sceIntraParaDivCPU_M[2] = k1_Intra_Div;
	mechPara_M.sceIntraParaDivCPU_M[3] = k2_Intra_Div;
	mechPara_M.sceIntraParaDivCPU_M[4] = intraDivEffectiveRange;

	//////////////////////
	//// Block 5 /////////
	//////////////////////
	double bondAdhCriLen =
			globalConfigVars.getConfigValue("BondAdhCriLen").toDouble();
	mechPara_M.bondAdhCriLenCPU_M = bondAdhCriLen;

	double bondStiff = globalConfigVars.getConfigValue("BondStiff").toDouble();
	mechPara_M.bondStiffCPU_M = bondStiff;

	//Ali June 16
	double bondStiff_Mitotic = globalConfigVars.getConfigValue("BondStiff_Mitotic").toDouble();
	mechPara_M.bondStiffCPU_Mitotic = bondStiff_Mitotic;

	double growthPrgrCriVal = globalConfigVars.getConfigValue(
			"GrowthPrgrCriVal").toDouble();
	mechPara_M.growthPrgrCriValCPU_M = growthPrgrCriVal;
	double maxAdhBondLen =
			globalConfigVars.getConfigValue("MaxAdhBondLen").toDouble();
	mechPara_M.maxAdhBondLenCPU_M = maxAdhBondLen;
	double minAdhBondLen =
			globalConfigVars.getConfigValue("MinAdhBondLen").toDouble();
	mechPara_M.minAdhBondLenCPU_M = minAdhBondLen;
}

void SceNodes::debugNAN() {
	uint totalActiveNodeC = allocPara_M.currentActiveCellCount
			* allocPara_M.maxAllNodePerCell;
	double res = thrust::reduce(infoVecs.nodeLocX.begin(),
			infoVecs.nodeLocX.begin() + totalActiveNodeC);

	if (isnan(res)) {
		std::cout << "fatal error! NAN found" << std::endl;
		std::cout.flush();
		exit(0);
	}
}

std::vector<std::pair<uint, uint> > SceNodes::obtainPossibleNeighborPairs_M() {
	std::vector<std::pair<uint, uint> > result;
	thrust::host_vector<uint> keyBeginCPU = auxVecs.keyBegin;
	thrust::host_vector<uint> keyEndCPU = auxVecs.keyEnd;
	thrust::host_vector<uint> bucketKeysCPU = auxVecs.bucketKeys;
	thrust::host_vector<uint> bucketValuesCPU = auxVecs.bucketValues;
	thrust::host_vector<uint> bucketValuesExtendedCPU =
			auxVecs.bucketValuesIncludingNeighbor;
	uint iterationCounter = 0;

	uint maxNodePerCell = allocPara_M.maxAllNodePerCell;
	uint offSet = allocPara_M.bdryNodeCount;
	uint memThreshold = allocPara_M.maxMembrNodePerCell;
	int size = bucketKeysCPU.size();

	int node1, node2, cellRank1, cellRank2, nodeRank1, nodeRank2;
	for (int i = 0; i < size; i++) {
		for (int j = keyBeginCPU[bucketKeysCPU[i]];
				j < keyEndCPU[bucketKeysCPU[i]]; j++) {
			node1 = bucketValuesCPU[i];
			node2 = bucketValuesExtendedCPU[j];
			if (node1 >= node2) {
				continue;
			} else {
				cellRank1 = (node1 - offSet) / maxNodePerCell;
				nodeRank1 = (node1 - offSet) % maxNodePerCell;
				cellRank2 = (node2 - offSet) / maxNodePerCell;
				nodeRank2 = (node2 - offSet) % maxNodePerCell;
				if (nodeRank1 >= memThreshold && nodeRank2 >= memThreshold
						&& cellRank1 == cellRank2) {
					result.push_back(std::make_pair<uint, uint>(node1, node2));
				}
			}
			iterationCounter++;
		}
	}
	return result;
}

void SceNodes::initValues(std::vector<CVector>& initBdryCellNodePos,
		std::vector<CVector>& initProfileNodePos,
		std::vector<CVector>& initCartNodePos,
		std::vector<CVector>& initECMNodePos,
		std::vector<CVector>& initFNMCellNodePos,
		std::vector<CVector>& initMXCellNodePos) {

	uint FNMNodeCount = initFNMCellNodePos.size();
	uint MXNodeCount = initMXCellNodePos.size();

	uint beginAddressOfProfile = allocPara.startPosProfile;
	uint beginAddressOfCart = allocPara.startPosCart;
// find the begining position of ECM.
	uint beginAddressOfECM = allocPara.startPosECM;
// find the begining position of FNM cells.
	uint beginAddressOfFNM = allocPara.startPosCells;
// find the begining position of MX cells.
	uint beginAddressOfMX = beginAddressOfFNM + FNMNodeCount;

	std::vector<double> initBdryCellNodePosX = getArrayXComp(
			initBdryCellNodePos);
	thrust::copy(initBdryCellNodePosX.begin(), initBdryCellNodePosX.end(),
			infoVecs.nodeLocX.begin());
	std::vector<double> initBdryCellNodePosY = getArrayYComp(
			initBdryCellNodePos);
	thrust::copy(initBdryCellNodePosY.begin(), initBdryCellNodePosY.end(),
			infoVecs.nodeLocY.begin());

// copy x and y position of nodes of Profile to actual node position.
	std::vector<double> initProfileNodePosX = getArrayXComp(initProfileNodePos);
	thrust::copy(initProfileNodePosX.begin(), initProfileNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfProfile);
	std::vector<double> initProfileNodePosY = getArrayYComp(initProfileNodePos);
	thrust::copy(initProfileNodePosY.begin(), initProfileNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfProfile);

// copy x and y position of nodes of Profile to actual node position.
	std::vector<double> initCartNodePosX = getArrayXComp(initCartNodePos);
	thrust::copy(initCartNodePosX.begin(), initCartNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfCart);
	std::vector<double> initCartNodePosY = getArrayYComp(initCartNodePos);
	thrust::copy(initCartNodePosY.begin(), initCartNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfCart);

// copy x and y position of nodes of ECM to actual node position.
	std::vector<double> initECMNodePosX = getArrayXComp(initECMNodePos);
	thrust::copy(initECMNodePosX.begin(), initECMNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfECM);
	std::vector<double> initECMNodePosY = getArrayYComp(initECMNodePos);
	thrust::copy(initECMNodePosY.begin(), initECMNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfECM);

	for (int i = 0; i < initECMNodePosX.size(); i++) {
		assert(infoVecs.nodeLocX[i + beginAddressOfECM] == initECMNodePosX[i]);
		assert(!isnan(initECMNodePosX[i]));
	}

// copy x and y position of nodes of FNM cells to actual node position.
	std::vector<double> initFNMCellNodePosX = getArrayXComp(initFNMCellNodePos);
	thrust::copy(initFNMCellNodePosX.begin(), initFNMCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfFNM);
	std::vector<double> initFNMCellNodePosY = getArrayYComp(initFNMCellNodePos);
	thrust::copy(initFNMCellNodePosY.begin(), initFNMCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfFNM);

	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfFNM,
			infoVecs.nodeCellType.begin() + beginAddressOfMX, FNM);

// copy x and y position of nodes of MX cells to actual node position.
	std::vector<double> initMXCellNodePosX = getArrayXComp(initMXCellNodePos);
	thrust::copy(initMXCellNodePosX.begin(), initMXCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfMX);
	std::vector<double> initMXCellNodePosY = getArrayYComp(initMXCellNodePos);
	thrust::copy(initMXCellNodePosY.begin(), initMXCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfMX);

	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfMX,
			infoVecs.nodeCellType.begin() + beginAddressOfMX + MXNodeCount, MX);
}
// It copies the information of node locations from CPU to GPU
void SceNodes::initValues_M(std::vector<bool>& initIsActive,
		std::vector<CVector>& initCellNodePos,
		std::vector<CVector>& initCellNodeMultip_actomyo,
		std::vector<CVector>& initCellNodeMultip_integrin,
		std::vector<SceNodeType>& nodeTypes,
		std::vector<double>& mDppV,
		std::vector<MembraneType1>& mTypeV) {

	std::vector<double> initCellNodePosX = getArrayXComp(initCellNodePos);

	std::vector<double> initCellNodePosY = getArrayYComp(initCellNodePos);

	// std::vector<double> initCellNodePosZ = getArrayZComp(initCellNodePos);

	std::vector<double> initCellNodeActomyosin_basal = getArrayXActomyosin(initCellNodeMultip_actomyo);
	std::vector<double> initCellNodeActomyosin_apical = getArrayYActomyosin(initCellNodeMultip_actomyo);
	std::cout<<"size of initCellNodeActomyosin = "<<initCellNodeActomyosin_basal.size()<< std::endl;
	std::vector<double> initCellNodeIntegrin = getArrayXIntegrin(initCellNodeMultip_integrin);
	// std::vector<double> initCellNodeIntegrin_apical = getArrayYIntegrin(initCellNodeMultip_integrin);
	std::cout<<"size of initCellNodeIntegrin = "<<initCellNodeIntegrin.size()<<std::endl;
	thrust::copy(initCellNodePosX.begin(), initCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + allocPara_M.bdryNodeCount);
	thrust::copy(initCellNodePosY.begin(), initCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + allocPara_M.bdryNodeCount);
	// thrust::copy(initCellNodePosZ.begin(), initCellNodePosZ.end(),
	// 		infoVecs.nodeLocZ.begin() + allocPara_M.bdryNodeCount);
	
	thrust::copy(initCellNodeActomyosin_basal.begin(), initCellNodeActomyosin_basal.end(),
			infoVecs.nodeActomyosinMultip_basal.begin() + allocPara_M.bdryNodeCount);
	thrust::copy(initCellNodeActomyosin_apical.begin(), initCellNodeActomyosin_apical.end(),
			infoVecs.nodeActomyosinMultip_apical.begin() + allocPara_M.bdryNodeCount);
			
	thrust::copy(initCellNodeIntegrin.begin(), initCellNodeIntegrin.end(),
			infoVecs.nodeIntegrinMultip.begin() + allocPara_M.bdryNodeCount);
	// thrust::copy(initCellNodeIntegrin_apical.begin(), initCellNodeIntegrin_apical.end(),
			// infoVecs.nodeIntegrinMultip_apical.begin() + allocPara_M.bdryNodeCount);
			
	thrust::copy(nodeTypes.begin(), nodeTypes.end(),
			infoVecs.nodeCellType.begin() + allocPara_M.bdryNodeCount);
	thrust::copy(mDppV.begin(), mDppV.end(),
			infoVecs.dppLevel.begin() ); // Ali 
	thrust::copy(mTypeV.begin(), mTypeV.end(),
			infoVecs.memNodeType1.begin() ); // Ali 
	thrust::copy(initIsActive.begin(), initIsActive.end(),
			infoVecs.nodeIsActive.begin() + allocPara_M.bdryNodeCount);

}

VtkAnimationData SceNodes::obtainAnimationData(AnimationCriteria aniCri) {
	VtkAnimationData vtkData;
	std::vector<std::pair<uint, uint> > pairs = obtainPossibleNeighborPairs();
	cout << "size of potential pairs = " << pairs.size() << endl;

// unordered_map is more efficient than map, but it is a c++ 11 feature
// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

// Doesn't have to copy the entire nodeLocX array.
// Only copy the first half will be sufficient
	thrust::host_vector<double> hostTmpVectorLocX = infoVecs.nodeLocX;
	thrust::host_vector<double> hostTmpVectorLocY = infoVecs.nodeLocY;
	thrust::host_vector<double> hostTmpVectorLocZ = infoVecs.nodeLocZ;

	thrust::host_vector<double> hostTmpVectorForceX;
	thrust::host_vector<double> hostTmpVectorForceY;
	thrust::host_vector<double> hostTmpVectorForceZ;
	thrust::host_vector<double> hostTmpVectorVelVal;

	assert(hostTmpVectorLocX.size() == hostTmpVectorLocY.size());
	assert(hostTmpVectorLocY.size() == hostTmpVectorLocZ.size());

	thrust::host_vector<SceNodeType> hostTmpVectorNodeType =
			infoVecs.nodeCellType;
	thrust::host_vector<uint> hostTmpVectorNodeRank = infoVecs.nodeCellRank;
	thrust::host_vector<double> hostTmpVectorNodeStress;

	if (aniCri.animationType != CellType) {
		hostTmpVectorForceX = infoVecs.nodeInterForceX;
		hostTmpVectorForceY = infoVecs.nodeInterForceY;
		hostTmpVectorForceZ = infoVecs.nodeInterForceZ;

		assert(hostTmpVectorForceX.size() == hostTmpVectorLocX.size());
		assert(hostTmpVectorForceX.size() == hostTmpVectorForceY.size());
		assert(hostTmpVectorForceX.size() == hostTmpVectorForceZ.size());

		uint vecSize = hostTmpVectorForceX.size();
		hostTmpVectorVelVal.resize(vecSize);
		for (uint i = 0; i < vecSize; i++) {
			hostTmpVectorVelVal[i] = sqrt(
					hostTmpVectorForceX[i] * hostTmpVectorForceX[i]
							+ hostTmpVectorForceY[i] * hostTmpVectorForceY[i]
							+ hostTmpVectorForceZ[i] * hostTmpVectorForceZ[i]);
		}

	}
	if (aniCri.animationType == Force) {
		vtkData.isArrowIncluded = true;
	} else {
		vtkData.isArrowIncluded = false;
	}

	uint curIndex = 0;
	for (uint i = 0; i < pairs.size(); i++) {
		uint node1Index = pairs[i].first;
		uint node2Index = pairs[i].second;
		double node1X = hostTmpVectorLocX[node1Index];
		double node1Y = hostTmpVectorLocY[node1Index];
		double node1Z = hostTmpVectorLocZ[node1Index];
		SceNodeType node1T = hostTmpVectorNodeType[node1Index];
		uint node1R = hostTmpVectorNodeRank[node1Index];
		double node2X = hostTmpVectorLocX[node2Index];
		double node2Y = hostTmpVectorLocY[node2Index];
		double node2Z = hostTmpVectorLocZ[node2Index];
		SceNodeType node2T = hostTmpVectorNodeType[node2Index];
		uint node2R = hostTmpVectorNodeRank[node2Index];

		if (aniCri.isPairQualify(node1Index, node2Index, node1X, node1Y, node1Z,
				node1T, node1R, node2X, node2Y, node2Z, node2T, node2R)) {
			IndexMap::iterator it = locIndexToAniIndexMap.find(pairs[i].first);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].first, curIndex));
				curIndex++;
				PointAniData ptAniData;
				if (aniCri.animationType == ForceAbsVal) {
					ptAniData.colorScale = hostTmpVectorVelVal[node1Index];
				} else if (aniCri.animationType == Force) {
					ptAniData.colorScale = hostTmpVectorVelVal[node1Index];
					if (hostTmpVectorVelVal[node1Index] > aniCri.threshold) {
						ptAniData.dir.x = hostTmpVectorForceX[node1Index]
								/ hostTmpVectorVelVal[node1Index]
								* aniCri.arrowLength;
						ptAniData.dir.y = hostTmpVectorForceY[node1Index]
								/ hostTmpVectorVelVal[node1Index]
								* aniCri.arrowLength;
						ptAniData.dir.z = hostTmpVectorForceZ[node1Index]
								/ hostTmpVectorVelVal[node1Index]
								* aniCri.arrowLength;
					} else {
						ptAniData.dir.x = 0;
						ptAniData.dir.y = 0;
						ptAniData.dir.z = 0;
					}
				} else {
					ptAniData.colorScale = nodeTypeToScale(node1T);
				}
				ptAniData.pos = CVector(node1X, node1Y, node1Z);
				vtkData.pointsAniData.push_back(ptAniData);
			}
			it = locIndexToAniIndexMap.find(pairs[i].second);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].second, curIndex));
				curIndex++;
				PointAniData ptAniData;
				if (aniCri.animationType == ForceAbsVal) {
					ptAniData.colorScale = hostTmpVectorVelVal[node2Index];
				} else if (aniCri.animationType == Force) {
					ptAniData.colorScale = hostTmpVectorVelVal[node2Index];
					if (hostTmpVectorVelVal[node2Index] > aniCri.threshold) {
						ptAniData.dir.x = hostTmpVectorForceX[node2Index]
								/ hostTmpVectorVelVal[node2Index]
								* aniCri.arrowLength;
						ptAniData.dir.y = hostTmpVectorForceY[node2Index]
								/ hostTmpVectorVelVal[node2Index]
								* aniCri.arrowLength;
						ptAniData.dir.z = hostTmpVectorForceZ[node2Index]
								/ hostTmpVectorVelVal[node2Index]
								* aniCri.arrowLength;
					} else {
						ptAniData.dir.x = 0;
						ptAniData.dir.y = 0;
						ptAniData.dir.z = 0;
					}
				} else {
					ptAniData.colorScale = nodeTypeToScale(node2T);
				}
				ptAniData.pos = CVector(node2X, node2Y, node2Z);
				vtkData.pointsAniData.push_back(ptAniData);
			}

			it = locIndexToAniIndexMap.find(pairs[i].first);
			uint aniIndex1 = it->second;
			it = locIndexToAniIndexMap.find(pairs[i].second);
			uint aniIndex2 = it->second;

			LinkAniData linkData;
			linkData.node1Index = aniIndex1;
			linkData.node2Index = aniIndex2;
			vtkData.linksAniData.push_back(linkData);
		}
	}

	uint profileStartIndex = allocPara.startPosProfile;
	uint profileEndIndex = profileStartIndex
			+ allocPara.currentActiveProfileNodeCount;

	for (uint i = profileStartIndex; i < profileEndIndex; i++) {
		PointAniData ptAniData;
		ptAniData.pos = CVector(hostTmpVectorLocX[i], hostTmpVectorLocY[i],
				hostTmpVectorLocZ[i]);

		if (aniCri.animationType == ForceAbsVal) {
			ptAniData.colorScale = hostTmpVectorVelVal[i];
		} else if (aniCri.animationType == Force) {
			ptAniData.colorScale = hostTmpVectorVelVal[i];
			if (hostTmpVectorVelVal[i] > aniCri.threshold) {
				ptAniData.dir.x = hostTmpVectorForceX[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
				ptAniData.dir.y = hostTmpVectorForceY[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
				ptAniData.dir.z = hostTmpVectorForceZ[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
			}
		} else {
			ptAniData.colorScale = nodeTypeToScale(hostTmpVectorNodeType[i]);
		}
		vtkData.pointsAniData.push_back(ptAniData);
		LinkAniData linkData;
		linkData.node1Index = curIndex;
		linkData.node2Index = curIndex + 1;
		if (i != profileEndIndex - 1) {
			vtkData.linksAniData.push_back(linkData);
		}
		curIndex++;
	}

	uint cartStartIndex = allocPara.startPosCart;
	uint cartEndIndex = cartStartIndex + allocPara.maxCartNodeCount;
	for (uint i = cartStartIndex; i < cartEndIndex; i++) {
		bool isActive = infoVecs.nodeIsActive[i];
		if (!isActive) {
			continue;
		}

		PointAniData ptAniData;
		ptAniData.pos = CVector(hostTmpVectorLocX[i], hostTmpVectorLocY[i],
				hostTmpVectorLocZ[i]);
		if (aniCri.animationType == ForceAbsVal) {
			ptAniData.colorScale = hostTmpVectorVelVal[i];
		} else if (aniCri.animationType == Force) {
			ptAniData.colorScale = hostTmpVectorVelVal[i];
			if (hostTmpVectorVelVal[i] > aniCri.threshold) {
				ptAniData.dir.x = hostTmpVectorForceX[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
				ptAniData.dir.y = hostTmpVectorForceY[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
				ptAniData.dir.z = hostTmpVectorForceZ[i]
						/ hostTmpVectorVelVal[i] * aniCri.arrowLength;
			}
		} else {
			ptAniData.colorScale = nodeTypeToScale(hostTmpVectorNodeType[i]);
		}
		vtkData.pointsAniData.push_back(ptAniData);

		bool isNextActive;
		if (i == cartEndIndex - 1) {
			isNextActive = false;
		} else {
			isNextActive = infoVecs.nodeIsActive[i + 1];
		}
		if (isNextActive) {
			LinkAniData linkData;
			linkData.node1Index = curIndex;
			linkData.node2Index = curIndex + 1;
			vtkData.linksAniData.push_back(linkData);
		}
		curIndex++;
	}

	return vtkData;
}

// TODO
VtkAnimationData SceNodes::obtainAnimationData_M(AnimationCriteria aniCri) {
	VtkAnimationData vtkData;
	std::vector<std::pair<uint, uint> > pairs = obtainPossibleNeighborPairs_M();
	cout << "size of potential pairs = " << pairs.size() << endl;

// unordered_map is more efficient than map, but it is a c++ 11 feature
// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

// Doesn't have to copy the entire nodeLocX array.
// Only copy the first half will be sufficient
	thrust::host_vector<double> hostTmpVectorLocX = infoVecs.nodeLocX;
	thrust::host_vector<double> hostTmpVectorLocY = infoVecs.nodeLocY;
	thrust::host_vector<bool> hostIsActiveVec = infoVecs.nodeIsActive;
	thrust::host_vector<int> hostBondVec = infoVecs.nodeAdhereIndex;
	thrust::host_vector<double> hostMembrTenMag = infoVecs.membrTensionMag;
	thrust::host_vector<SceNodeType> hostTmpVectorNodeType =
			infoVecs.nodeCellType;

	uint activeCellCount = allocPara_M.currentActiveCellCount;
	uint maxNodePerCell = allocPara_M.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_M.maxMembrNodePerCell;
	uint beginIndx = allocPara_M.bdryNodeCount;
//uint endIndx = beginIndx + activeCellCount * maxNodePerCell;

//uint cellRank1, nodeRank1, cellRank2, nodeRank2;
	uint index1;
	int index2;
	std::vector<BondInfo> bondInfoVec;

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
	vtkData.bondsInfo = bondInfoVec;

	uint curIndex = 0;
	for (uint i = 0; i < pairs.size(); i++) {
		uint node1Index = pairs[i].first;
		uint node2Index = pairs[i].second;
		double node1X = hostTmpVectorLocX[node1Index];
		double node1Y = hostTmpVectorLocY[node1Index];

		double node2X = hostTmpVectorLocX[node2Index];
		double node2Y = hostTmpVectorLocY[node2Index];

		if (aniCri.isPairQualify_M(node1X, node1Y, node2X, node2Y)) {
			IndexMap::iterator it = locIndexToAniIndexMap.find(pairs[i].first);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].first, curIndex));
				curIndex++;
				PointAniData ptAniData;
				//ptAniData.colorScale = nodeTypeToScale(
				//		hostTmpVectorNodeType[node1Index]);
				ptAniData.colorScale = -1;
				ptAniData.colorScale2 = -1;//AAMIRI
				ptAniData.pos = CVector(node1X, node1Y, 0);
				vtkData.pointsAniData.push_back(ptAniData);
			}
			it = locIndexToAniIndexMap.find(pairs[i].second);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].second, curIndex));
				curIndex++;
				PointAniData ptAniData;
				//ptAniData.colorScale = nodeTypeToScale(
				//		hostTmpVectorNodeType[node1Index]);
				ptAniData.colorScale = -1;
				ptAniData.colorScale2 = -1;//AAMIRI
				ptAniData.pos = CVector(node2X, node2Y, 0);
				vtkData.pointsAniData.push_back(ptAniData);
			}

			it = locIndexToAniIndexMap.find(pairs[i].first);
			uint aniIndex1 = it->second;
			it = locIndexToAniIndexMap.find(pairs[i].second);
			uint aniIndex2 = it->second;

			LinkAniData linkData;
			linkData.node1Index = aniIndex1;
			linkData.node2Index = aniIndex2;
			vtkData.linksAniData.push_back(linkData);
		}
	}

	return vtkData;
}

void SceNodes::findBucketBounds() {
	thrust::counting_iterator<unsigned int> search_begin(0);
	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), search_begin,
			search_begin + domainPara.totalBucketCount,
			auxVecs.keyBegin.begin());
	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), search_begin,
			search_begin + domainPara.totalBucketCount, auxVecs.keyEnd.begin());
}

void SceNodes::findBucketBounds_M() {
	thrust::counting_iterator<uint> search_begin(0);
	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExtProc_M, search_begin,
			search_begin + domainPara.totalBucketCount,
			auxVecs.keyBegin.begin());
	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExtProc_M, search_begin,
			search_begin + domainPara.totalBucketCount, auxVecs.keyEnd.begin());
}

void SceNodes::findBucketBounds3D() {
	thrust::counting_iterator<uint> search_begin(0);
	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExtProc_M, search_begin,
			search_begin + domainPara.totalBucketCount,
			auxVecs.keyBegin.begin());
	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExtProc_M, search_begin,
			search_begin + domainPara.totalBucketCount, auxVecs.keyEnd.begin());
}

void SceNodes::prepareSceForceComputation() {
	buildBuckets2D();
	extendBuckets2D();
	findBucketBounds();
}

void SceNodes::prepareSceForceComputation_M() {
	buildBuckets2D_M();
	extendBuckets2D_M();
	findBucketBounds_M();
}

void SceNodes::prepareSceForceComputation3D() {
	buildBuckets3D();
	extendBuckets3D();
	findBucketBounds3D();
}

void SceNodes::addNewlyDividedCells(
		thrust::device_vector<double> &nodeLocXNewCell,
		thrust::device_vector<double> &nodeLocYNewCell,
		thrust::device_vector<double> &nodeLocZNewCell,
		thrust::device_vector<bool> &nodeIsActiveNewCell,
		thrust::device_vector<SceNodeType> &nodeCellTypeNewCell) {

// data validation
	uint nodesSize = nodeLocXNewCell.size();
	assert(nodesSize % allocPara.maxNodeOfOneCell == 0);
	uint addCellCount = nodesSize / allocPara.maxNodeOfOneCell;

// position that we will add newly divided cells.
	uint shiftStartPosNewCell = allocPara.startPosCells
			+ allocPara.currentActiveCellCount * allocPara.maxNodeOfOneCell;

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodeLocXNewCell.begin(),
							nodeLocYNewCell.begin(), nodeLocZNewCell.begin(),
							nodeIsActiveNewCell.begin(),
							nodeCellTypeNewCell.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodeLocXNewCell.end(),
							nodeLocYNewCell.end(), nodeLocZNewCell.end(),
							nodeIsActiveNewCell.end(),
							nodeCellTypeNewCell.end())),
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(),
							infoVecs.nodeCellType.begin()))
					+ shiftStartPosNewCell);

	// total number of cells has increased.
	allocPara.currentActiveCellCount = allocPara.currentActiveCellCount
			+ addCellCount;
}

void SceNodes::buildBuckets2D() {
	int totalActiveNodes;
	if (controlPara.simuType != Disc_M) {
		totalActiveNodes = allocPara.startPosCells
				+ allocPara.currentActiveCellCount * allocPara.maxNodeOfOneCell;
	} else {
		totalActiveNodes = allocPara_M.bdryNodeCount
				+ allocPara_M.currentActiveCellCount
						* allocPara_M.maxAllNodePerCell;
	}

	auxVecs.bucketKeys.resize(totalActiveNodes);
	auxVecs.bucketValues.resize(totalActiveNodes);
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(totalActiveNodes);

// takes counting iterator and coordinates
// return tuple of keys and values
// transform the points to their bucket indices
	thrust::transform(
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), countingIterBegin)),
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), countingIterBegin))
					+ totalActiveNodes,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			pointToBucketIndex2D(domainPara.minX, domainPara.maxX,
					domainPara.minY, domainPara.maxY, domainPara.gridSpacing));

// sort the points by their bucket index
	thrust::sort_by_key(auxVecs.bucketKeys.begin(), auxVecs.bucketKeys.end(),
			auxVecs.bucketValues.begin());
// for those nodes that are inactive, key value of UINT_MAX will be returned.
// we need to removed those keys along with their values.
	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.end(), UINT_MAX);

	auxVecs.bucketKeys.erase(auxVecs.bucketKeys.end() - numberOfOutOfRange,
			auxVecs.bucketKeys.end());
	auxVecs.bucketValues.erase(auxVecs.bucketValues.end() - numberOfOutOfRange,
			auxVecs.bucketValues.end());
}

void SceNodes::buildBuckets2D_M() {
	int totalActiveNodes = allocPara_M.bdryNodeCount
			+ allocPara_M.currentActiveCellCount
					* allocPara_M.maxAllNodePerCell;

	thrust::counting_iterator<uint> iBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices
	thrust::transform(
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), iBegin)),
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), iBegin))
					+ totalActiveNodes,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			pointToBucketIndex2D(domainPara.minX, domainPara.maxX,
					domainPara.minY, domainPara.maxY, domainPara.gridSpacing));

	// sort the points by their bucket index
	thrust::sort_by_key(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + totalActiveNodes,
			auxVecs.bucketValues.begin());
	// for those nodes that are inactive, key value of UINT_MAX will be returned.
	// we need to removed those keys along with their values.
	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + totalActiveNodes, UINT_MAX);

	endIndx_M = totalActiveNodes - numberOfOutOfRange;
}

void SceNodes::buildBuckets3D() {
	int totalActiveNodes = allocPara_M.bdryNodeCount
			+ allocPara_M.currentActiveCellCount
					* allocPara_M.maxAllNodePerCell;

	thrust::counting_iterator<uint> iBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices
	thrust::transform(
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), iBegin)),
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), iBegin))
					+ totalActiveNodes,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			BucketIndexer3D(domainPara.minX, domainPara.maxX, domainPara.minY,
					domainPara.maxY, domainPara.minZ, domainPara.maxZ,
					domainPara.gridSpacing));

	// sort the points by their bucket index
	thrust::sort_by_key(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + totalActiveNodes,
			auxVecs.bucketValues.begin());
	// for those nodes that are inactive, key value of UINT_MAX will be returned.
	// we need to removed those keys along with their values.
	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + totalActiveNodes, UINT_MAX);

	endIndx_M = totalActiveNodes - numberOfOutOfRange;

}

__device__
double computeDist(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2)
					+ (zPos - zPos2) * (zPos - zPos2));
}

__device__
double computeDist2D(double &xPos, double &yPos, double &xPos2, double &yPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2));
}

__device__
void calculateAndAddECMForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {

	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceECMPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceECMPara[0] / sceECMPara[2]
				* exp(-linkLength / sceECMPara[2])
				+ sceECMPara[1] / sceECMPara[3]
						* exp(-linkLength / sceECMPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.3;
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}
__device__
void calculateAndAddProfileForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	forceValue = -sceProfilePara[5] * (linkLength - sceProfilePara[6]);

	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}

__device__
void calculateAndAddIntraForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue;
	if (linkLength > sceIntraPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceIntraPara[0] / sceIntraPara[2]
				* exp(-linkLength / sceIntraPara[2])
				+ sceIntraPara[1] / sceIntraPara[3]
						* exp(-linkLength / sceIntraPara[3]);
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}

__device__
void calAndAddIntraForceDiv(double& xPos, double& yPos, double& zPos,
		double& xPos2, double& yPos2, double& zPos2, double& growPro,
		double& xRes, double& yRes, double& zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue;
	if (linkLength > sceIntraPara[4]) {
		forceValue = 0;
	} else {
		if (growPro > sceIntraParaDiv[4]) {
			double intraPara0 = growPro * (sceIntraParaDiv[0])
					+ (1.0 - growPro) * sceIntraPara[0];
			double intraPara1 = growPro * (sceIntraParaDiv[1])
					+ (1.0 - growPro) * sceIntraPara[1];
			double intraPara2 = growPro * (sceIntraParaDiv[2])
					+ (1.0 - growPro) * sceIntraPara[2];
			double intraPara3 = growPro * (sceIntraParaDiv[3])
					+ (1.0 - growPro) * sceIntraPara[3];
			forceValue = -intraPara0 / intraPara2
					* exp(-linkLength / intraPara2)
					+ intraPara1 / intraPara3 * exp(-linkLength / intraPara3);
		} else {
			forceValue = -sceIntraPara[0] / sceIntraPara[2]
					* exp(-linkLength / sceIntraPara[2])
					+ sceIntraPara[1] / sceIntraPara[3]
							* exp(-linkLength / sceIntraPara[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}

__device__
void calAndAddIntraDiv_M(double& xPos, double& yPos, double& xPos2,
		double& yPos2, double& growPro, double& xRes, double& yRes) {
	double linkLength = computeDist2D(xPos, yPos, xPos2, yPos2);
	double forceValue;
	if (growPro > growthPrgrCriVal_M) {
		if (linkLength > sceIntraParaDiv_M[4]) {
			forceValue = 0;
		} else {
			double percent = (growPro - growthPrgrCriVal_M)
					/ (1.0 - growthPrgrCriVal_M);
			double intraPara0 = percent * (sceIntraParaDiv_M[0])
					+ (1.0 - percent) * sceIntraPara_M[0];
			double intraPara1 = percent * (sceIntraParaDiv_M[1])
					+ (1.0 - percent) * sceIntraPara_M[1];
			double intraPara2 = percent * (sceIntraParaDiv_M[2])
					+ (1.0 - percent) * sceIntraPara_M[2];
			double intraPara3 = percent * (sceIntraParaDiv_M[3])
					+ (1.0 - percent) * sceIntraPara_M[3];
			forceValue = -intraPara0 / intraPara2
					* exp(-linkLength / intraPara2)
					+ intraPara1 / intraPara3 * exp(-linkLength / intraPara3);
		}
	} else {
		if (linkLength > sceIntraPara_M[4]) {
			forceValue = 0;
		} else {
			forceValue = -sceIntraPara_M[0] / sceIntraPara_M[2]
					* exp(-linkLength / sceIntraPara_M[2])
					+ sceIntraPara_M[1] / sceIntraPara_M[3]
							* exp(-linkLength / sceIntraPara_M[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddIntraB_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes) {
	double linkLength = computeDist2D(xPos, yPos, xPos2, yPos2);
	double forceValue;
	if (linkLength > sceIntnlBPara_M[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceIntnlBPara_M[0] / sceIntnlBPara_M[2]
				* exp(-linkLength / sceIntnlBPara_M[2])
				+ sceIntnlBPara_M[1] / sceIntnlBPara_M[3]
						* exp(-linkLength / sceIntnlBPara_M[3]);
	}
	//if (forceValue > 0) {
	//	forceValue = 0;
	//}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddInter_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes) {
	double linkLength = computeDist2D(xPos, yPos, xPos2, yPos2);
	double forceValue;
	if (linkLength > sceInterBPara_M[4]) {
		forceValue = 0;
	} else {
		
			forceValue = -sceInterBPara_M[0] / sceInterBPara_M[2]
					* exp(-linkLength / sceInterBPara_M[2])
					+ sceInterBPara_M[1] / sceInterBPara_M[3]
							* exp(-linkLength / sceInterBPara_M[3]);
		
	//	if (forceValue > 0) {  //Ali 
	//		forceValue = 0;
	//	}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddInter_M_mitotic(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes, double scaling) {
	double linkLength = computeDist2D(xPos, yPos, xPos2, yPos2);
	double forceValue;
	if (linkLength > sceInterBPara_M[4]) {
		forceValue = 0;
	} else {
				forceValue = -(1.0+(2.5-1.0)*scaling)*sceInterBPara_M[0] / (/*(1*(1-scaling)+scaling*(1.0/2.0))**/sceInterBPara_M[2])
					* exp(-linkLength / (/*(1*(1-scaling)+scaling*(1.0/2.0))**/sceInterBPara_M[2]))
					+ (1.0+(2.5-1.0)*scaling)*sceInterBPara_M[1] / (/*(1*(1-scaling)+scaling*(1.0/2.0))**/sceInterBPara_M[3])
							* exp(-linkLength / (/*(1*(1-scaling)+scaling*(1.0/2.0))**/sceInterBPara_M[3]));
				// forceValue = -sceInterBPara_M[0] / (sceInterBPara_M[2])
				// 	* exp(-linkLength / (sceInterBPara_M[2]))
				// 	+ sceInterBPara_M[1] / (sceInterBPara_M[3])
				// 			* exp(-linkLength / (sceInterBPara_M[3]));
		
	//	if (forceValue > 0) {  //Ali 
	//		forceValue = 0;
	//	}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddInter_M2(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes) {
	double linkLength = computeDist2D(xPos, yPos, xPos2, yPos2);
	double forceValue;
	if (linkLength > sceInterBPara_Jones_M[2]) {
		forceValue = 0;
	} else {
		forceValue =24*sceInterBPara_Jones_M[0]/linkLength*pow(sceInterBPara_Jones_M[1]/linkLength,6)*
                                                          ( 1.0-2 *pow(sceInterBPara_Jones_M[1]/linkLength,6) ) ;
			
				
						
		if (forceValue > 0) {
			forceValue = 0;
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}
//Ali

__device__
void calculateAndAddInterForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}

__device__
void calAndAddInterForceDisc(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes, double& interForceX, double& interForceY,
		double& interForceZ) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
	}
	double fX = forceValue * (xPos2 - xPos) / linkLength;
	double fY = forceValue * (yPos2 - yPos) / linkLength;
	double fZ = forceValue * (zPos2 - zPos) / linkLength;
	xRes = xRes + fX;
	yRes = yRes + fY;
	zRes = zRes + fZ;
	interForceX = interForceX + fX;
	interForceY = interForceY + fY;
	interForceZ = interForceZ + fZ;
}

__device__
void calculateAndAddCartForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceCartPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceCartPara[0] / sceCartPara[2]
				* exp(-linkLength / sceCartPara[2])
				+ sceCartPara[1] / sceCartPara[3]
						* exp(-linkLength / sceCartPara[3]);
		if (linkLength > 1.0e-12) {
			//double dotProduct = (xPos2 - xPos) / linkLength * cartGrowDirVec[0]
			//		+ (yPos2 - yPos) / linkLength * cartGrowDirVec[1]
			//		+ (zPos2 - zPos) / linkLength * cartGrowDirVec[2];
			//forceValue = forceValue * dotProduct;
			// this is just a temperary solution -- the direction should not be fixed.
			xRes = xRes - forceValue * cartGrowDirVec[0];
			yRes = yRes - forceValue * cartGrowDirVec[1];
			zRes = zRes - forceValue * cartGrowDirVec[2];
			//xRes = xRes + forceValue * (xPos2 - xPos);
			//yRes = yRes + forceValue * (yPos2 - yPos);
			//zRes = zRes + forceValue * (zPos2 - zPos);
		}
		if (forceValue > 0) {
			//forceValue = forceValue * 0.01;
			forceValue = 0;
			//xRes = xRes + forceValue * (xPos2 - xPos);
			//yRes = yRes + forceValue * (yPos2 - yPos);
			//zRes = zRes + forceValue * (zPos2 - zPos);
		}
	}

}

__device__
void calculateAndAddDiffInterCellForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterDiffPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterDiffPara[0] / sceInterDiffPara[2]
				* exp(-linkLength / sceInterDiffPara[2])
				+ sceInterDiffPara[1] / sceInterDiffPara[3]
						* exp(-linkLength / sceInterDiffPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.2;
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}
__device__
void calculateAndAddInterForceDiffType(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.3;
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}

__device__ bool bothNodesCellNode(uint nodeGlobalRank1, uint nodeGlobalRank2,
		uint cellNodesThreshold) {
	if (nodeGlobalRank1 < cellNodesThreshold
			&& nodeGlobalRank2 < cellNodesThreshold) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isSameCell(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos
			|| nodeGlobalRank2 < cellNodeBeginPos) {
		return false;
	}
	if ((nodeGlobalRank1 - cellNodeBeginPos) / nodeCountPerCell
			== (nodeGlobalRank2 - cellNodeBeginPos) / nodeCountPerCell) {
		return true;
	} else {
		return false;
	}
}

//Ali
__device__

bool Is_Lennard_Jones() {
            if (sceInterBPara_Jones_On_M==1) {
                 return true ; 
                     }
            else    {
                  return false ; 
                    }
 }

__device__
bool isSameCell_m(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos_M
			|| nodeGlobalRank2 < cellNodeBeginPos_M) {
		return false;
	}
	if ((nodeGlobalRank1 - cellNodeBeginPos_M) / allNodeCountPerCell_M
			== (nodeGlobalRank2 - cellNodeBeginPos_M) / allNodeCountPerCell_M) {
		return true;
	} else {
		return false;
	}
}

__device__
bool bothInternal(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos_M
			|| nodeGlobalRank2 < cellNodeBeginPos_M) {
		return false;
	}
	uint nodeRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	uint nodeRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	if (nodeRank1 >= membrThreshold_M && nodeRank2 >= membrThreshold_M) {
		return true;
	} else {
		return false;
	}
}

__device__
bool bothMembr(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos_M
			|| nodeGlobalRank2 < cellNodeBeginPos_M) {
		return false;
	}
	uint nodeRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	uint nodeRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	if (nodeRank1 < membrThreshold_M && nodeRank2 < membrThreshold_M) {
		return true;
	} else {
		return false;
	}
}

__device__
bool bothMembrDiffCell(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos_M
			|| nodeGlobalRank2 < cellNodeBeginPos_M) {
		return false;
	}
	uint cellRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			/ allNodeCountPerCell_M;
	uint cellRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			/ allNodeCountPerCell_M;
	if (cellRank1 == cellRank2) {
		return false;
	}
	uint nodeRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	uint nodeRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	if (nodeRank1 < membrThreshold_M && nodeRank2 < membrThreshold_M) {
		return true;
	} else {
		return false;
	}
}

//AAMIRI
/*
__device__
bool isNodeOnMembrane(uint nodeGlobalRank) {

	uint nodeRank = (nodeGlobalRank - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;

	if (nodeGlobalRank >= cellNodeBeginPos_M && nodeRank < membrThreshold_M){
		return true;
	} else{
		return false;
	}

}
*/

__device__
bool sameCellMemIntnl(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos_M
			|| nodeGlobalRank2 < cellNodeBeginPos_M) {
		return false;
	}
	uint cellRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			/ allNodeCountPerCell_M;
	uint cellRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			/ allNodeCountPerCell_M;
	if (cellRank1 != cellRank2) {
		return false;
	}
	uint nodeRank1 = (nodeGlobalRank1 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	uint nodeRank2 = (nodeGlobalRank2 - cellNodeBeginPos_M)
			% allNodeCountPerCell_M;
	if ((nodeRank1 < membrThreshold_M && nodeRank2 >= membrThreshold_M)
			|| (nodeRank2 < membrThreshold_M && nodeRank1 >= membrThreshold_M)) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isSameECM(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if ((nodeGlobalRank1 - ECMbeginPos) / nodeCountPerECM
			== (nodeGlobalRank2 - ECMbeginPos) / nodeCountPerECM) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isNeighborECMNodes(uint nodeGlobalRank1, uint nodeGlobalRank2) {
// this means that two nodes are from the same ECM
	if ((nodeGlobalRank1 - ECMbeginPos) / nodeCountPerECM
			== (nodeGlobalRank2 - ECMbeginPos) / nodeCountPerECM) {
		// this means that two nodes are actually close to each other
		// seems to be strange because of unsigned int.
		if ((nodeGlobalRank1 > nodeGlobalRank2
				&& nodeGlobalRank1 - nodeGlobalRank2 == 1)
				|| (nodeGlobalRank2 > nodeGlobalRank1
						&& nodeGlobalRank2 - nodeGlobalRank1 == 1)) {
			return true;
		}
	}
	return false;
}

__device__ bool isNeighborProfileNodes(uint nodeGlobalRank1,
		uint nodeGlobalRank2) {
	if ((nodeGlobalRank1 > nodeGlobalRank2
			&& nodeGlobalRank1 - nodeGlobalRank2 == 1)
			|| (nodeGlobalRank2 > nodeGlobalRank1
					&& nodeGlobalRank2 - nodeGlobalRank1 == 1)) {
		return true;
	}
	return false;
}

__device__ bool ofSameType(uint cellType1, uint cellType2) {
	if (cellType1 == cellType2) {
		return true;
	} else {
		return false;
	}
}

__device__ bool bothCellNodes(SceNodeType &type1, SceNodeType &type2) {
	if ((type1 == MX || type1 == FNM) && (type2 == MX || type2 == FNM)) {
		return true;
	} else {
		return false;
	}
}

__device__
void attemptToAdhere(bool& isSuccess, uint& index, double& dist,
		uint& nodeRank2, double& xPos1, double& yPos1, double& xPos2,
		double& yPos2) {
	double length = computeDist2D(xPos1, yPos1, xPos2, yPos2);
	if (length <= bondAdhCriLen_M) {
		if (isSuccess) {
			if (length < dist) {
				dist = length;
				index = nodeRank2;
			}
		} else {
			isSuccess = true;
			index = nodeRank2;
			dist = length;
		}
	}
}

__device__
void handleAdhesionForce_M(int& adhereIndex, double& xPos, double& yPos,
		double& curAdherePosX, double& curAdherePosY, double& xRes,
		double& yRes, double& alpha, double & beta) {
	double curLen = computeDist2D(xPos, yPos, curAdherePosX, curAdherePosY);
	//if (curLen > maxAdhBondLen_M) {
	//	adhereIndex = -1;
	//	return;
//	} else {
		if (curLen > minAdhBondLen_M) {
			double forceValue = beta*(curLen - minAdhBondLen_M) * (bondStiff_M * alpha + bondStiff_Mitotic * (1.0-alpha) );
			xRes = xRes + forceValue * (curAdherePosX - xPos) / curLen;
			yRes = yRes + forceValue * (curAdherePosY - yPos) / curLen;
		}

//	}
}

//Ali for reaction force

__device__
void handleAdhesionForce_M2(double& xPos, double& yPos,
		double& curAdherePosX, double& curAdherePosY, double& xRes,
		double& yRes, double& alpha) {
		double curLen = computeDist2D(xPos, yPos, curAdherePosX, curAdherePosY);
		if (curLen > minAdhBondLen_M ) {
			double forceValue = (curLen - minAdhBondLen_M) * (bondStiff_M * alpha + bondStiff_Mitotic * (1.0-alpha) );
			xRes = forceValue * (curAdherePosX - xPos) / curLen;
			yRes = forceValue * (curAdherePosY - yPos) / curLen;
		}
		else {
			xRes=0 ;
			yRes=0 ; 
		}

	}




//Ali June 16
__device__
double getMitoticAdhCoef(double& growProg, double& growProgNeigh){
	double alpha = 1.0;


	if (growProg > growthPrgrCriVal_M && growProgNeigh > growthPrgrCriVal_M){
			alpha = 1.0 - ( 0.5*(growProg+growProgNeigh)-growthPrgrCriVal_M )/(1.0 - growthPrgrCriVal_M);
		//	adhSkipped = true;
		}
	else if (growProg > growthPrgrCriVal_M){
			alpha = 1.0 - (growProg-growthPrgrCriVal_M)/(1.0 - growthPrgrCriVal_M);
		//	adhSkipped = true;
		}
	else if (growProgNeigh > growthPrgrCriVal_M){
			alpha = 1.0 - (growProgNeigh-growthPrgrCriVal_M)/(1.0 - growthPrgrCriVal_M);
		//	adhSkipped = true;
		}


	return alpha;
}

__device__
void calculateForceBetweenLinkNodes(double &xLoc, double &yLoc, double &zLoc,
		double &xLocLeft, double &yLocLeft, double &zLocLeft, double &xLocRight,
		double &yLocRight, double &zLocRight, double &xVel, double &yVel,
		double &zVel) {
	double linkLengthLeft = computeDist(xLoc, yLoc, zLoc, xLocLeft, yLocLeft,
			zLocLeft);
	double forceValueLeft = sceProfilePara[5]
			* (linkLengthLeft - sceProfilePara[6]);
	xVel = xVel + forceValueLeft * (xLocLeft - xLoc) / linkLengthLeft;
	yVel = yVel + forceValueLeft * (yLocLeft - yLoc) / linkLengthLeft;
	zVel = zVel + forceValueLeft * (zLocLeft - zLoc) / linkLengthLeft;

	double linkLengthRight = computeDist(xLoc, yLoc, zLoc, xLocRight, yLocRight,
			zLocRight);
	double forceValueRight = sceProfilePara[5]
			* (linkLengthRight - sceProfilePara[6]);
	xVel = xVel + forceValueRight * (xLocRight - xLoc) / linkLengthRight;
	yVel = yVel + forceValueRight * (yLocRight - yLoc) / linkLengthRight;
	zVel = zVel + forceValueRight * (zLocRight - zLoc) / linkLengthRight;

}

__device__
void handleSceForceNodesBasic(uint& nodeRank1, uint& nodeRank2, double& xPos,
		double& yPos, double& zPos, double& xPos2, double& yPos2, double& zPos2,
		double& xRes, double& yRes, double& zRes, double* _nodeLocXAddress,
		double* _nodeLocYAddress, double* _nodeLocZAddress) {
	if (isSameCell(nodeRank1, nodeRank2)) {
		calculateAndAddIntraForce(xPos, yPos, zPos, _nodeLocXAddress[nodeRank2],
				_nodeLocYAddress[nodeRank2], _nodeLocZAddress[nodeRank2], xRes,
				yRes, zRes);
	} else {
		calculateAndAddInterForce(xPos, yPos, zPos, _nodeLocXAddress[nodeRank2],
				_nodeLocYAddress[nodeRank2], _nodeLocZAddress[nodeRank2], xRes,
				yRes, zRes);
	}
}

__device__
void handleSceForceNodesDisc(uint& nodeRank1, uint& nodeRank2, double& xPos,
		double& yPos, double& zPos, double& xPos2, double& yPos2, double& zPos2,
		double& xRes, double& yRes, double& zRes, double& interForceX,
		double& interForceY, double& interForceZ, double* _nodeLocXAddress,
		double* _nodeLocYAddress, double* _nodeLocZAddress,
		double* _nodeGrowProAddr) {
	if (isSameCell(nodeRank1, nodeRank2)) {
		calAndAddIntraForceDiv(xPos, yPos, zPos, _nodeLocXAddress[nodeRank2],
				_nodeLocYAddress[nodeRank2], _nodeLocZAddress[nodeRank2],
				_nodeGrowProAddr[nodeRank2], xRes, yRes, zRes);
	} else {
		calAndAddInterForceDisc(xPos, yPos, zPos, _nodeLocXAddress[nodeRank2],
				_nodeLocYAddress[nodeRank2], _nodeLocZAddress[nodeRank2], xRes,
				yRes, zRes, interForceX, interForceY, interForceZ);
	}
}

// __device__
// void handleSceForceNodesDisc_M(uint& nodeRank1, uint& nodeRank2, double& xPos,
// 		double& yPos, double& xPos2, double& yPos2, double& xRes, double& yRes,
// 		double* _nodeLocXAddress, double* _nodeLocYAddress,
// 		double* _nodeGrowProAddr) {

// 	if (isSameCell_m(nodeRank1, nodeRank2)) {
// 		if (bothInternal(nodeRank1, nodeRank2)) {
// 			// both nodes are internal type.
// 			calAndAddIntraDiv_M(xPos, yPos, _nodeLocXAddress[nodeRank2],
// 					_nodeLocYAddress[nodeRank2], _nodeGrowProAddr[nodeRank2],
// 					xRes, yRes);
// 		} else if (bothMembr(nodeRank1, nodeRank2)) {
// 			// both nodes epithilium type. no sce force applied.
// 			// nothing to do here.
// 		} else {
// 			// one node is epithilium type the other is internal type.
// 			calAndAddIntraB_M(xPos, yPos, _nodeLocXAddress[nodeRank2],
// 					_nodeLocYAddress[nodeRank2], xRes, yRes);
// 		}
// 	} else {
// 		if (bothMembr(nodeRank1, nodeRank2)) {
// 			calAndAddInter_M(xPos, yPos, _nodeLocXAddress[nodeRank2],
// 					_nodeLocYAddress[nodeRank2], xRes, yRes);
// 		}
// 	}
// }

void SceNodes::extendBuckets2D() {
	static const uint extensionFactor2D = 9;
	uint valuesCount = auxVecs.bucketValues.size();
	auxVecs.bucketKeysExpanded.resize(valuesCount * extensionFactor2D);
	auxVecs.bucketValuesIncludingNeighbor.resize(
			valuesCount * extensionFactor2D);

	/**
	 * beginning of constant iterator
	 */
	thrust::constant_iterator<uint> first(extensionFactor2D);
	/**
	 * end of constant iterator.
	 * the plus sign only indicate movement of position, not value.
	 * e.g. movement is 5 and first iterator is initialized as 9
	 * result array is [9,9,9,9,9];
	 */
	thrust::constant_iterator<uint> last = first + valuesCount;

	expand(first, last,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							auxVecs.bucketValuesIncludingNeighbor.begin())));

	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd = countingBegin
			+ valuesCount * extensionFactor2D;

	thrust::transform(
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.end(), countingEnd)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			NeighborFunctor2D(domainPara.XBucketSize, domainPara.YBucketSize));

	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), UINT_MAX);

	int sizeBeforeShrink = auxVecs.bucketKeysExpanded.size();
	int numberInsideRange = sizeBeforeShrink - numberOfOutOfRange;
	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(),
			auxVecs.bucketValuesIncludingNeighbor.begin());
	auxVecs.bucketKeysExpanded.erase(
			auxVecs.bucketKeysExpanded.begin() + numberInsideRange,
			auxVecs.bucketKeysExpanded.end());
	auxVecs.bucketValuesIncludingNeighbor.erase(
			auxVecs.bucketValuesIncludingNeighbor.begin() + numberInsideRange,
			auxVecs.bucketValuesIncludingNeighbor.end());
}

void SceNodes::extendBuckets2D_M() {
	endIndxExt_M = endIndx_M * 9;
	/**
	 * beginning of constant iterator
	 */
	thrust::constant_iterator<uint> first(9);
	/**
	 * end of constant iterator.
	 * the plus sign only indicate movement of position, not value.
	 * e.g. movement is 5 and first iterator is initialized as 9
	 * result array is [9,9,9,9,9];
	 */
	thrust::constant_iterator<uint> last = first + endIndx_M;

	expand(first, last,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							auxVecs.bucketValuesIncludingNeighbor.begin())));

	thrust::counting_iterator<uint> countingBegin(0);

	thrust::transform(
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)) + endIndxExt_M,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			NeighborFunctor2D(domainPara.XBucketSize, domainPara.YBucketSize));

	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExt_M, UINT_MAX);

	endIndxExtProc_M = endIndxExt_M - numberOfOutOfRange;
	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExt_M,
			auxVecs.bucketValuesIncludingNeighbor.begin());
}

void SceNodes::extendBuckets3D() {
	endIndxExt_M = endIndx_M * 27;
	/**
	 * beginning of constant iterator
	 */
	thrust::constant_iterator<uint> first(27);
	/**
	 * end of constant iterator.
	 * the plus sign only indicate movement of position, not value.
	 * e.g. movement is 5 and first iterator is initialized as 9
	 * result array is [9,9,9,9,9];
	 */
	thrust::constant_iterator<uint> last = first + endIndx_M; // this is NOT numerical addition!

	expand(first, last,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							auxVecs.bucketValuesIncludingNeighbor.begin())));

	thrust::counting_iterator<uint> countingBegin(0);

	thrust::transform(
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)) + endIndxExt_M,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			NgbrFunc3D(domainPara.XBucketSize, domainPara.YBucketSize,
					domainPara.ZBucketSize));

	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExt_M, UINT_MAX);

	endIndxExtProc_M = endIndxExt_M - numberOfOutOfRange;
	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.begin() + endIndxExt_M,
			auxVecs.bucketValuesIncludingNeighbor.begin());
}

void SceNodes::applySceForcesBasic() {
	uint* valueAddress = thrust::raw_pointer_cast(
			&auxVecs.bucketValuesIncludingNeighbor[0]);
	double* nodeLocXAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocY[0]);
	double* nodeLocZAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocZ[0]);

	thrust::transform(
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.begin()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.begin()),
							auxVecs.bucketValues.begin(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.begin()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.end()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.end()),
							auxVecs.bucketValues.end(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.end()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(infoVecs.nodeVelX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelZ.begin(),
									auxVecs.bucketValues.begin()))),
			AddSceForceBasic(valueAddress, nodeLocXAddress, nodeLocYAddress,
					nodeLocZAddress));
}

void SceNodes::applySceForcesDisc() {
	uint* valueAddress = thrust::raw_pointer_cast(
			&auxVecs.bucketValuesIncludingNeighbor[0]);
	double* nodeLocXAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocY[0]);
	double* nodeLocZAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocZ[0]);
	double* nodeGrowProAddr = thrust::raw_pointer_cast(
			&infoVecs.nodeGrowPro[0]);

	thrust::transform(
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.begin()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.begin()),
							auxVecs.bucketValues.begin(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.begin()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.end()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.end()),
							auxVecs.bucketValues.end(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.end()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(infoVecs.nodeVelX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelZ.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(
									infoVecs.nodeInterForceX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(
									infoVecs.nodeInterForceY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(
									infoVecs.nodeInterForceZ.begin(),
									auxVecs.bucketValues.begin()))),
			AddSceForceDisc(valueAddress, nodeLocXAddress, nodeLocYAddress,
					nodeLocZAddress, nodeGrowProAddr));
}

void SceNodes::applySceForcesDisc_M(double mitoticThreshold) {

	ECellType eCellTypeTmp ; 
	thrust::host_vector <ECellType> eCellTypeVHost ;

	eCellTypeVHost.resize(allocPara_M.currentActiveCellCount, notActive) ; 
	//for (int i= 0 ; i<allocPara_M.currentActiveCellCount;  i++) {
	//	eCellTypeTmp=cellsSceNodes->getCellInfoVecs().eCellTypeV2[i];
	//	cout << "Epithelial cell type is="<<eCellTypeTmp <<endl ; 
		//eCellTypeVHost.push_back(eCellTypeTmp) ; 
//	}
	if (adhUpdate == true){
		// std::cout<<"Entering adhesion pair between neighboring cells update step!"<<std::endl;
	}
	if (adhUpdate) {
		adhUpdate=false ;
		int maxNumAdh=180 ;
		//vector <ECellType> eCellTypeV2Host ;

     	thrust :: copy (infoVecs.nodeLocX.begin(),infoVecs.nodeLocX.end(),infoVecs.nodeLocXHost.begin()) ; // Ali	
		
     	thrust :: copy (infoVecs.nodeLocY.begin(),infoVecs.nodeLocY.end(),infoVecs.nodeLocYHost.begin()) ; // Ali 	
		
     	thrust :: copy (infoVecs.nodeIsActive.begin(),infoVecs.nodeIsActive.end(),infoVecs.nodeIsActiveHost.begin()) ; // Ali 	
		
     	thrust :: copy (infoVecs.nodeCellRankFront.begin() ,infoVecs.nodeCellRankFront.end() ,infoVecs.nodeCellRankFrontHost.begin()) ; // Ali 	
		
     	thrust :: copy (infoVecs.nodeCellRankBehind.begin(),infoVecs.nodeCellRankBehind.end(),infoVecs.nodeCellRankBehindHost.begin()) ; // Ali 	
		
     	thrust :: copy (infoVecs.memNodeType1.begin(),infoVecs.memNodeType1.end(),infoVecs.memNodeType1Host.begin()) ; // Ali 
		
		thrust :: copy (cellsSceNodes->getCellInfoVecs().eCellTypeV2.begin(),cellsSceNodes->getCellInfoVecs().eCellTypeV2.begin()+allocPara_M.currentActiveCellCount,eCellTypeVHost.begin()) ;
	 	
		thrust::fill(infoVecs.nodeAdhereIndexHost.begin(),infoVecs.nodeAdhereIndexHost.end(), -1) ;
		
	 	thrust::fill(infoVecs.nodeMemMirrorIndexHost.begin(),infoVecs.nodeMemMirrorIndexHost.end(), -1) ;  //Ali it is important to reset the values
		
	 	//thrust::fill(infoVecs.nodeIsLateralMemHost.begin(),infoVecs.nodeIsLateralMemHost.end(), false) ;  //Ali
	 	thrust::fill(infoVecs.nodeAdhMinDist.begin(),infoVecs.nodeAdhMinDist.end(), 10000) ;  //Ali
		
		// std::cout<<"AdhUpdate ERROR 1"<<std::endl;

        int totalActiveNodes = allocPara_M.currentActiveCellCount* allocPara_M.maxAllNodePerCell; // Ali
	  	int maxMembNode=    allocPara_M.maxMembrNodePerCell ; 
		
	  	int maxNodePerCell= allocPara_M.maxAllNodePerCell ; 
      	double  distMinP2,distP2 ;
	  	int indexAdhNode ; 
		bool findAnyNode ;
	  	double maxAdhLen= mechPara_M.bondAdhCriLenCPU_M; 
	  	int cellRankTmp1, cellRankTmp2 ; 
	  	int deactiveIdMyPair, deactiveIdAdhPair ;
		int activeMemCount [ allocPara_M.currentActiveCellCount] ;
		int firstApiLat [ allocPara_M.currentActiveCellCount] ; 
		int secondApiLat [ allocPara_M.currentActiveCellCount] ; 
		int cellRank, iNext, jJunction ;
		std::vector <SubApicalInfoEachCell> subApicalInfo ; 
		
		// std::cout<<"AdhUpdate ERROR 2"<<std::endl;

		//setup required basic parameters 
        for (int i=0 ; i< allocPara_M.currentActiveCellCount ; i++ ){
			activeMemCount[i] = 0 ; 
		}
		// std::cout<<"AdhUpdate ERROR 3"<<std::endl;

		for (int i=0 ; i<totalActiveNodes ;  i++) {
			infoVecs.isSubApicalJunctionHost[i]=false ; 
		}
		// std::cout<<"AdhUpdate ERROR 3"<<std::endl;
		for (int i=0 ; i<totalActiveNodes ;  i++) {
			if (infoVecs.nodeIsActiveHost[i]==true && (i%maxNodePerCell)<maxMembNode){
				cellRank=i/maxNodePerCell ; 
				activeMemCount [cellRank]=activeMemCount [cellRank]+1 ; 
			}
		}
	// std::cout<<"AdhUpdate ERROR 4"<<std::endl;
		subApicalInfo.clear() ; 
		

		//Find the subapical nodes in front of the cell
		int cellRankOld=-1 ; 
	
		for (int i=0 ; i<totalActiveNodes ;  i++) {
			if (i == totalActiveNodes-1){
				// std::cout<<"I am at the last step of this loop 1"<<std::endl;
			}
			// if (i/maxNodePerCell == 31 || i/maxNodePerCell==86){
			// 	std::cout<<"i/maxNodePerCell for subapical front = "<<i/maxNodePerCell<<std::endl;
			// 	std::cout<<"infoVecs.nodeIsActiveHost["<<i<<"] = "<<infoVecs.nodeIsActiveHost[i]<<std::endl;
			// 	std::cout<<"i%maxNodePerCell = "<< (i%maxNodePerCell)<<std::endl;
			// }
				if (infoVecs.nodeIsActiveHost[i]==true && (i%maxNodePerCell)<maxMembNode){ // check active and membrane
					// std::cout<<"IF "<<i<<std::endl;
					cellRank=i/maxNodePerCell ;

					eCellTypeTmp=eCellTypeVHost[cellRank];
					iNext=i+1 ; 
					if ( (i%maxNodePerCell)==(activeMemCount[cellRank]-1)) {  // if the node is the last node of cell's membrane
						iNext=iNext-activeMemCount [cellRank] ;
					}
					if ( (infoVecs.memNodeType1Host[i]==lateralA && infoVecs.memNodeType1Host[iNext]==apical1) ) { // find the apical junction
						firstApiLat[cellRank]=i ; // lateral node 
						// std::cout<<"1cellRank["<<cellRank<<"], condition met for 1st adhesion joints"<<std::endl;
						for (int j=0 ; j<NumAdhAfter(cellRank,eCellTypeTmp) ; j++) {   //find junction nodes // 
							jJunction=firstApiLat[cellRank]-j ; 
							if (jJunction <(cellRank*maxNodePerCell)) {
								jJunction=jJunction + activeMemCount [cellRank] ;
								//cout << " The subApicalNodes of cell rank " << cellRank << " passed the first node ID" << endl ; 
							}
		 					infoVecs.isSubApicalJunctionHost[jJunction]=true ;

							if (cellRank !=cellRankOld) {

						//		cout << " for cell rank= " << cellRank << " subapicalInfo has been created." << endl ;
								SubApicalInfoEachCell subApicalInfoEachCell(maxNumAdh); 
								subApicalInfo.push_back(subApicalInfoEachCell); 
								cellRankOld=cellRank ; 
							}
							subApicalInfo[cellRank].nodeIdFront[j]=jJunction ; 
							// if (cellRank >= 85){
							// 	std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
							// }

						}

					}
					else if (infoVecs.memNodeType1Host[iNext]==lateralA && infoVecs.memNodeType1Host[i]==apical1){
						// std::cout<<"ELSE IF "<<i<<std::endl;
						// std::cout<<"2cellRank["<<cellRank<<"], condition met for 1st adhesion joints"<<std::endl;
						// firstApiLat[cellRank]=i ; // lateral node
						firstApiLat[cellRank]=iNext ; // lateral node 
						for (int j=0 ; j<NumAdhAfter(cellRank,eCellTypeTmp) ; j++) {   //find junction nodes // 
							// jJunction=firstApiLat[cellRank]-j ;
							jJunction=firstApiLat[cellRank]+j ; 
							// if (jJunction <(cellRank*maxNodePerCell)) {
							// if (jJunction >= (cellRank*maxNodePerCell + activeMemCount[cellRank])) {
							if (jJunction >= ((cellRank)*maxNodePerCell + activeMemCount[cellRank])) {
								// jJunction=jJunction + activeMemCount [cellRank] ;
								jJunction=jJunction - activeMemCount [cellRank] ;
								//cout << " The subApicalNodes of cell rank " << cellRank << " passed the first node ID" << endl ; 
							}
						
		 					infoVecs.isSubApicalJunctionHost[jJunction]=true ;

							if (cellRank !=cellRankOld) {

						//		cout << " for cell rank= " << cellRank << " subapicalInfo has been created." << endl ;
								SubApicalInfoEachCell subApicalInfoEachCell(maxNumAdh); 
								subApicalInfo.push_back(subApicalInfoEachCell); 
								cellRankOld=cellRank ; 
							}
							subApicalInfo[cellRank].nodeIdFront[j]=jJunction ; 
							// if (cellRank >= 85){
							// 	std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
							// }

						}
					}
				
				}
		}
		// std::cout << "     >>> 14 <<<" << endl;
		// cout << "first set of adhesion joints are found" << endl ; 
		// if (cellRank == 86){
		// 	for (int j = 0; j < maxNumAdh; j++){
		// 		std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdFront["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdFront[j]<<std::endl;
		// 	}
		// 	for (int j = 0; j < maxNumAdh; j++){
		// 		std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
		// 	}
		// }
		
		//Find the subapical nodes supposingly behind (Before) the cell
		
		for (int i=0 ; i<totalActiveNodes ;  i++) {
			if (i == totalActiveNodes-1){
				// std::cout<<"I am at the last step of this loop 2"<<std::endl;
			}
			// if (i/maxNodePerCell > 86){
			// 	std::cout<<"cellRank index exceed 86, that's it"<<std::endl;
			// }
				if (infoVecs.nodeIsActiveHost[i]==true && (i%maxNodePerCell)<maxMembNode){
					cellRank=i/maxNodePerCell ;
					// if (cellRank == 86){
					// 	std::cout<<"i/maxNodePerCell = "<< i/maxNodePerCell<<" cellRank = "<<cellRank<<std::endl;
					// }
					//eCellType= eCellTypeV2Host[cellRank];
					eCellTypeTmp= eCellTypeVHost [cellRank];
					// if (cellRank == 31 || cellRank == 86){
					// 	std::cout<<"2nd eCellTypeTmp["<<cellRank<<"] = "<<eCellTypeTmp<<std::endl;
					// }
					iNext=i+1 ; 
					// if (cellRank == 86){// || cellRank == 31){
					// 	std::cout<<"i = "<<i<<", iNext = "<<iNext<<std::endl;
					// 	std::cout<<"i%maxNodePerCell = "<<i%maxNodePerCell<<std::endl;
					// 	std::cout<<"activeMemCount["<<cellRank<<"] = "<<activeMemCount[cellRank]<<std::endl;
					// }
					if ( (i%maxNodePerCell)==(activeMemCount [cellRank]-1)) {
						iNext=iNext-activeMemCount [cellRank]  ; 
						// if (cellRank >= 85){
						// 	std::cout<<"activeMemCount[cellRank] = "<<activeMemCount[cellRank]<<" altered iNext = "<<iNext<<std::endl;
						// }
					}
					// if (cellRank == 86){
					// 	std::cout<<" ^^^ 1 ^^^"<<std::endl;
					// 	for (int j = 0; j < maxNumAdh; j++){
					// 		std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdFront["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdFront[j]<<std::endl;
					// 	}
					// 	for (int j = 0; j < maxNumAdh; j++){
					// 		std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
					// 	}
					// }
					// if (cellRank == 86|| cellRank == 31){
					// 	std::cout<<"cellRank["<<cellRank<<"], infoVecs.memNodeType1Host["<<i<<"] = "<<infoVecs.memNodeType1Host[i]<<std::endl;
					// 	std::cout<<"cellRank["<<cellRank<<"], infoVecs.memNodeType1Host["<<iNext<<"] = "<<infoVecs.memNodeType1Host[iNext]<<std::endl;
					// 	// if (infoVecs.memNodeType1Host[i]==apical1 && infoVecs.memNodeType1Host[iNext]==lateralB
					// 	// 	|| (infoVecs.memNodeType1Host[iNext]==apical1 && ( infoVecs.memNodeType1Host[i]==lateralB ))){
					// 	// 	std::cout<<"condition met (secondApiLat)"<<std::endl;
					// 	// } 
					// 	// else{
					// 	// 	std::cout<<"condition not met (secondApiLat)"<<std::endl;
					// 	// }
					// }
					// if (infoVecs.memNodeType1Host[i]==apical1 && ( infoVecs.memNodeType1Host[iNext]==lateralB ) ) {

						
					if ((infoVecs.memNodeType1Host[i]==apical1 && ( infoVecs.memNodeType1Host[iNext]==lateralB )) ) {
						// std::cout<<"Cell Rank = "<<cellRank<<" This If?"<<std::endl;
						// std::cout<<"1cellRank["<<cellRank<<"], condition met for 2nd adhesion joints"<<std::endl;
						secondApiLat[cellRank]=iNext ; 
						// if (cellRank >= 85){
						// 	std::cout<<"secondApiLat["<<cellRank<<"] = "<<secondApiLat[cellRank]<<std::endl;
						// 	std::cout<<"NumAdhBefore = "<<NumAdhBefore(cellRank, eCellTypeTmp)<<std::endl;
						// }
						for (int j=0 ; j<NumAdhBefore(cellRank,eCellTypeTmp) ; j++) {   //find junction nodes
							jJunction=secondApiLat[cellRank]+j ; 
							// if (cellRank >= 85){
							// 	std::cout<<"j = "<<j<<" jJunction - activeMemCount["<<cellRank<<"] = "<<jJunction - activeMemCount[cellRank]<<std::endl;
							// }
							if (jJunction>=(cellRank*maxNodePerCell+activeMemCount [cellRank]) ) {
								jJunction=jJunction - activeMemCount [cellRank];

								//cout << " The subApicalNodes of cell rank " << cellRank << " passed the last node ID" << endl ; 
							}
		 					infoVecs.isSubApicalJunctionHost[jJunction]=true ;
							//  if (cellRank >= 85){
							//  	std::cout<<"infoVecs.isSubApicalJunctionHost["<<jJunction<<"] = "<<infoVecs.isSubApicalJunctionHost[jJunction]<<std::endl;
							//  }
							//  if (cellRank >= 85){
							// 	std::cout<<"subApicalInfo dimension = "<<subApicalInfo.size()<<std::endl;
							// 	std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
							// }
							subApicalInfo[cellRank].nodeIdBehind[j]=jJunction ; // the vector of structures for active cells has already been generated.
							// std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
	
						}
				
					}
					else if ( (infoVecs.memNodeType1Host[iNext]==apical1 && ( infoVecs.memNodeType1Host[i]==lateralB ))){
						// std::cout<<"Cell Rank = "<<cellRank<<" Or This Else If?"<<std::endl;
						// secondApiLat[cellRank]=iNext ;
						// std::cout<<"2cellRank["<<cellRank<<"], condition met for 2nd adhesion joints"<<std::endl;
						secondApiLat[cellRank]=i ; 
						// if (cellRank >= 85){
						// 	std::cout<<"secondApiLat["<<cellRank<<"] = "<<secondApiLat[cellRank]<<std::endl;
						// 	std::cout<<"NumAdhBefore = "<<NumAdhBefore(cellRank, eCellTypeTmp)<<std::endl;
						// }
						for (int j=0 ; j<NumAdhBefore(cellRank,eCellTypeTmp) ; j++) {   //find junction nodes
							jJunction=secondApiLat[cellRank]-j ; 
							// if (cellRank >= 85){
							// 	std::cout<<"j = "<<j<<" jJunction - activeMemCount["<<cellRank<<"] = "<<jJunction - activeMemCount[cellRank]<<std::endl;
							// }
							if (jJunction < (cellRank*maxNodePerCell) ) {
								jJunction=jJunction + activeMemCount [cellRank];

								//cout << " The subApicalNodes of cell rank " << cellRank << " passed the last node ID" << endl ; 
							}
							// std::cout<<"jJunction = "<<jJunction<<std::endl;
		 					infoVecs.isSubApicalJunctionHost[jJunction]=true ;
							//  if (cellRank == 86){
							// 	std::cout<<"j = "<<j<<std::endl;
							//  	std::cout<<"infoVecs.isSubApicalJunctionHost["<<jJunction<<"] = "<<infoVecs.isSubApicalJunctionHost[jJunction]<<std::endl;
							// 	std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
							//  }
							//  if (cellRank >= 85){
							// 	std::cout<<"subApicalInfo dimension = "<<subApicalInfo.size()<<std::endl;
							// 	std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
							// }
							subApicalInfo[cellRank].nodeIdBehind[j]=jJunction ; // the vector of structures for active cells has already been generated.
							// std::cout<<"subApicalInfo["<<cellRank<<"].nodeIdBehind["<<j<<"] = "<<subApicalInfo[cellRank].nodeIdBehind[j]<<std::endl;
	
						}
					}
				}
				
		}
		// std::cout << "     >>> 15 <<<" << endl;

		// cout << "Second set of adhesion joints are found" << endl ; 
		//for (int i=0 ; i<totalActiveNodes ;  i++) {
		//	if (infoVecs.isSubApicalJunctionHost[i]) {
		//		cout << "for cell with rank  "	<<int(i/maxNodePerCell) << "node rank of subApical junction is " << i << endl ;  
		 //	}
	//	}

		// cout << " size of vector storing information of apical junctions is " << subApicalInfo.size() << endl ; 
		if (subApicalInfo.size() != 0 ) {  // to pass the first time step in which the membrane node type is not defined.
			for ( int i= 0 ; i<allocPara_M.currentActiveCellCount ; i++) {
				
				for ( int j=0 ; j<maxNumAdh ; j++) {
			
					int idFront=subApicalInfo[i].nodeIdFront[j] ;
					// if (i == 31){
						// std::cout<<"ERROR 1"<<std::endl;
					// }	
					int idBehind=subApicalInfo[i].nodeIdBehind[j] ;
					// if (i == 31){
						// std::cout<<"ERROR 2"<<std::endl;
					// }	
					int cellRankFront=infoVecs.nodeCellRankFrontHost[i] ;
					// if (i == 31){
						// std::cout<<"ERROR 3"<<std::endl;
					// }	 
					int cellRankBehind=infoVecs.nodeCellRankBehindHost[i] ;
					// if (i == 31){
						// std::cout<<"ERROR 4"<<std::endl;
					// 	std::cout<<subApicalInfo[cellRankFront].nodeIdBehind.size()<<std::endl;
					// }	
					if (idFront != -1 && cellRankFront != -1) {
						infoVecs.nodeAdhereIndexHost[idFront]=subApicalInfo[cellRankFront].nodeIdBehind[j] ;
					}
					// if (i == 31){
						// std::cout<<"ERROR 5"<<std::endl;
					// }	
					if (idBehind !=-1 && cellRankBehind != -1) {
						infoVecs.nodeAdhereIndexHost[idBehind]=subApicalInfo[cellRankBehind].nodeIdFront[j] ;
					}
					// if (i == 31){
						// std::cout<<"ERROR 6"<<std::endl;
					// }	
					if ( eCellTypeVHost[i]==pouch && NumAdhBefore(i,pouch)==NumAdhAfter(i,pouch) ) {
						infoVecs.nodeMemMirrorIndexHost[idFront]=idBehind ;
						infoVecs.nodeMemMirrorIndexHost[idBehind]=idFront ;
					};
					// if (i == 31){
					// 	std::cout<<"ERROR 7"<<std::endl;
					// }	
					// if (i == 86 || i == 31){
					// 	std::cout<<"idFront["<<j<<"] = "<<idFront<<", loc = "<<infoVecs.nodeLocX[idFront]<<" "<<infoVecs.nodeLocY[idFront]<<", nodeAdhereIndexHost loc = "<<infoVecs.nodeLocX[infoVecs.nodeAdhereIndexHost[idFront]]
					// 			<<" "<<infoVecs.nodeLocY[infoVecs.nodeAdhereIndexHost[idFront]]<<" "<<infoVecs.nodeLocZ[infoVecs.nodeAdhereIndexHost[idFront]]<<std::endl;
					// 	std::cout<<"idBehind["<<j<<"] = "<<idBehind<<", loc = "<<infoVecs.nodeLocX[idBehind]<<" "<<infoVecs.nodeLocY[idBehind]<<", nodeAdhereIndexHost loc = "<<infoVecs.nodeLocX[infoVecs.nodeAdhereIndexHost[idBehind]]
					// 			<<" "<<infoVecs.nodeLocY[infoVecs.nodeAdhereIndexHost[idBehind]]<<" "<<infoVecs.nodeLocZ[infoVecs.nodeAdhereIndexHost[idBehind]]<<std::endl;
					// }
				}

				
			}

// std::cout << "     >>> 16 <<<" << endl;
	
/////////////////////////////////// start adhesion for apical nodes of pouch cells with apical nodes of peripodial cells ///////////////////////
			
	 		for (int i=0 ; i<totalActiveNodes ;  i++) {
				cellRankTmp1=i/maxNodePerCell ;
		 		distMinP2=10000 ; // large number
	  			findAnyNode=false ; 
				if (eCellTypeVHost[cellRankTmp1]==pouch && infoVecs.memNodeType1Host[i]==apical1) { 
		 			for (int j=0 ; j<totalActiveNodes ; j++) {
						cellRankTmp2=j/maxNodePerCell ; 
					//	if ( cellRankTmp2>=74 && cellRankTmp2<=76 && infoVecs.memNodeType1Host[i]==apical1) { 
						if (eCellTypeVHost[cellRankTmp2]==peri && infoVecs.memNodeType1Host[j]==apical1   ) {
							distP2=pow( infoVecs.nodeLocXHost[i]-infoVecs.nodeLocXHost[j],2)+
			         	    	   pow( infoVecs.nodeLocYHost[i]-infoVecs.nodeLocYHost[j],2) ;
							if (distP2<distMinP2   && distP2<maxAdhLen*maxAdhLen) {
							//	cout << " I am inside a function where there is one apical pouch and one apical perip node and it is min" << endl ; 
								distMinP2=distP2 ;
								indexAdhNode=j ; 
								findAnyNode=true ;
							}
		  				}
					}
		 	   	}
                
			   	if ( findAnyNode && sqrt(distMinP2)<infoVecs.nodeAdhMinDist[indexAdhNode]){
				//	cout << " I am inside apical adhesion" << endl ; 
	     			deactiveIdAdhPair=infoVecs.nodeAdhereIndexHost[indexAdhNode] ;
					if (deactiveIdAdhPair != -1){	
				//		cout << " I am inside deactiving one perip adhesion" << endl ; 
	     				infoVecs.nodeAdhereIndexHost[deactiveIdAdhPair]=-1 ;
	     				infoVecs.nodeAdhMinDist[deactiveIdAdhPair]=10000 ;
					}
	     			infoVecs.nodeAdhereIndexHost[i]=indexAdhNode ;
	     			infoVecs.nodeAdhereIndexHost[indexAdhNode]=i ; 
					infoVecs.nodeAdhMinDist[indexAdhNode]=sqrt(distMinP2) ; 
					infoVecs.nodeAdhMinDist[i]=sqrt(distMinP2) ;
				}

		 	}
		
		  	// cout << " I am ready to copy the data in adhesion function to the GPU " << endl ; 
	
/////////////////////////////////// start adhesion for apical nodes of pouch cells with apical nodes of peripodial cells ///////////////////////
  		} // finish if of bypassing the first time
		// copy back to GPU 
		
		thrust::copy(infoVecs.nodeAdhereIndexHost.begin(),infoVecs.nodeAdhereIndexHost.end(), infoVecs.nodeAdhereIndex.begin()) ;  //Ali
		
		thrust::copy(infoVecs.nodeMemMirrorIndexHost.begin(),infoVecs.nodeMemMirrorIndexHost.end(), infoVecs.nodeMemMirrorIndex.begin()) ;  //Ali
		
    	thrust::copy(infoVecs.isSubApicalJunctionHost.begin(),infoVecs.isSubApicalJunctionHost.end(), infoVecs.isSubApicalJunction.begin()) ;  //Ali
		// std::cout << "     >>> 17 <<<" << endl;// std::cout << "     >>> 16 <<<" << endl;
		for (int i = 0; i < infoVecs.nodeAdhereIndex.size(); i++){
			uint cellRank=i/maxNodePerCell ;
			if (cellRank >= allocPara.currentActiveCellCount){
				continue;
			}
			// std::cout<<"i = "<<i<<" node loc = "<<infoVecs.nodeLocX[i]<<" "<<infoVecs.nodeLocY[i]<<" "<<cellRank<<" "<<infoVecs.nodeMemMirrorIndex[i]<<" "<<infoVecs.nodeAdhereIndex[i]<<" "<<infoVecs.nodeCellType[i]<<" " <<infoVecs.memNodeType1[i]<<" "<<infoVecs.nodeIsActive[i]<< std::endl;
		}

	} // finish the if condition for finding the pair node

	uint* valueAddress = thrust::raw_pointer_cast(
			&auxVecs.bucketValuesIncludingNeighbor[0]);
	double* nodeLocXAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocY[0]);
	int* nodeAdhIdxAddress = thrust::raw_pointer_cast(
			&infoVecs.nodeAdhereIndex[0]);
	int* membrIntnlAddress = thrust::raw_pointer_cast(
			&infoVecs.membrIntnlIndex[0]);
	double* nodeGrowProAddr = thrust::raw_pointer_cast(
			&infoVecs.nodeGrowPro[0]);
	double* quiescencePerNodeAddr = thrust::raw_pointer_cast(
			&infoVecs.quiescencePerNode[0]);
	
	thrust::transform(
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.begin()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.begin()),
							auxVecs.bucketValues.begin(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.begin()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.begin() + endIndx_M),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.begin() + endIndx_M),
							auxVecs.bucketValues.end(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.begin() + endIndx_M),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.begin() + endIndx_M))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(infoVecs.nodeVelX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelY.begin(),
									auxVecs.bucketValues.begin()))),
			AddForceDisc_M(valueAddress, nodeLocXAddress, nodeLocYAddress,
					nodeAdhIdxAddress, membrIntnlAddress,
					nodeGrowProAddr, mitoticThreshold, quiescencePerNodeAddr));
	
}

const SceDomainPara& SceNodes::getDomainPara() const {
	return domainPara;
}

void SceNodes::setDomainPara(const SceDomainPara& domainPara) {
	this->domainPara = domainPara;
}

const NodeAllocPara& SceNodes::getAllocPara() const {
	return allocPara;
}

void SceNodes::setAllocPara(const NodeAllocPara& allocPara) {
	this->allocPara = allocPara;
}

const NodeAuxVecs& SceNodes::getAuxVecs() const {
	return auxVecs;
}

void SceNodes::setAuxVecs(const NodeAuxVecs& auxVecs) {
	this->auxVecs = auxVecs;
}

NodeInfoVecs& SceNodes::getInfoVecs() {
	return infoVecs;
}

std::vector<std::vector<int> > SceNodes::obtainLabelMatrix(
		PixelizePara& pixelPara) {
	std::vector<std::vector<int> > result;
	std::vector<NodeWithLabel> nodeLabels;
	ResAnalysisHelper resHelper;
	resHelper.setPixelPara(pixelPara);

	thrust::host_vector<double> hostTmpVectorLocX = infoVecs.nodeLocX;
	thrust::host_vector<double> hostTmpVectorLocY = infoVecs.nodeLocY;
	thrust::host_vector<double> hostTmpVectorLocZ = infoVecs.nodeLocZ;
	thrust::host_vector<SceNodeType> hostTmpVectorNodeType =
			infoVecs.nodeCellType;
	thrust::host_vector<uint> hostTmpVectorNodeRank = infoVecs.nodeCellRank;
	thrust::host_vector<uint> hostTmpVectorIsActive = infoVecs.nodeIsActive;

	uint startIndex = allocPara.startPosCells;
	uint endIndex = startIndex
			+ allocPara.currentActiveCellCount * allocPara.maxNodeOfOneCell;
	for (uint i = startIndex; i < endIndex; i++) {
		if (hostTmpVectorIsActive[i] == true) {
			NodeWithLabel nodeLabel;
			nodeLabel.cellRank = hostTmpVectorNodeRank[i];
			nodeLabel.position = CVector(hostTmpVectorLocX[i],
					hostTmpVectorLocY[i], hostTmpVectorLocZ[i]);
			nodeLabels.push_back(nodeLabel);
		}
	}

	result = resHelper.outputLabelMatrix(nodeLabels);
	return result;
}

void SceNodes::initControlPara(bool isStab) {
	int simuTypeConfigValue =
			globalConfigVars.getConfigValue("SimulationType").toInt();
	controlPara.simuType = parseTypeFromConfig(simuTypeConfigValue);
	controlPara.controlSwitchs.outputBmpImg = globalConfigVars.getSwitchState(
			"Switch_OutputBMP");
	controlPara.controlSwitchs.outputLabelMatrix =
			globalConfigVars.getSwitchState("Switch_OutputLabelMatrix");
	controlPara.controlSwitchs.outputStat = globalConfigVars.getSwitchState(
			"Switch_OutputStat");
	controlPara.controlSwitchs.outputVtkFile = globalConfigVars.getSwitchState(
			"Switch_OutputVtk");
	if (isStab) {
		controlPara.controlSwitchs.stab = ON;
	} else {
		controlPara.controlSwitchs.stab = OFF;
	}

}

void SceNodes::sceForcesPerfTesting() {
	prepareSceForceComputation();
	applySceForcesBasic();
}

void SceNodes::sceForcesPerfTesting_M() {
	prepareSceForceComputation_M();
	applySceForcesBasic_M();
}

void SceNodes::applySceForcesBasic_M() {
}

void SceNodes::sceForcesDisc() {
	prepareSceForceComputation();
	applySceForcesDisc();
}

//Node velocity is reset back to zero in this function at the beginning of every new time step
void SceNodes::sceForcesDisc_M(double timeRatio, double timeRatio_Crit_Division, int cycle,
								double contractActomyo_multip_perCell1, double contractActomyo_multip_perCell2, double contractActomyo_multip_perCell3,
								double contractActomyo_multip_perCell_apical1, double contractActomyo_multip_perCell_apical2, double contractActomyo_multip_perCell_apical3, double mitoticThreshold) {
#ifdef DebugMode
	cudaEvent_t start1, start2, start3, stop;
	float elapsedTime1, elapsedTime2, elapsedTime3;
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&start3);
	cudaEventCreate(&stop);
	cudaEventRecord(start1, 0);
#endif
	// cout << " confirm--- 1 ---" << endl;
	cout.flush();
	prepareSceForceComputation_M(); //buckets for optimization of searching algorithm

#ifdef DebugMode
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif
	if (timeRatio == 0 && cycle <= 0){
		// std::cout<<"Setting the anisotropic contractility profile"<<std::endl;
		// std::cout<<"This message should only appear once, if more than once is shown, something is wrong!"<<std::endl;
		// std::vector<double> contractActomyo_multip;
		// //Starting from cell 0 to cell 62
		// contractActomyo_multip.push_back(0.0);//CellID0
		// contractActomyo_multip.push_back(0.0);
		// contractActomyo_multip.push_back(0.536452492);
		// contractActomyo_multip.push_back(0.52731934 );
		// contractActomyo_multip.push_back(0.586284249);
		// contractActomyo_multip.push_back(0.634834161);
		// contractActomyo_multip.push_back(0.639320622);
		// contractActomyo_multip.push_back(0.690594456);
		// contractActomyo_multip.push_back(0.704374299);
		// contractActomyo_multip.push_back(0.817977888);
		// contractActomyo_multip.push_back(0.825348502);
		// contractActomyo_multip.push_back(0.831277039);
		// contractActomyo_multip.push_back(0.791539817);
		// contractActomyo_multip.push_back(0.845217113);
		// contractActomyo_multip.push_back(0.963467393);
		// contractActomyo_multip.push_back(0.949206858);
		// contractActomyo_multip.push_back(0.777279282);
		// contractActomyo_multip.push_back(0.796026278);
		// contractActomyo_multip.push_back(0.782887358);
		// contractActomyo_multip.push_back(0.810767505);
		// contractActomyo_multip.push_back(0.783528281);
		// contractActomyo_multip.push_back(0.654382311);
		// contractActomyo_multip.push_back(0.631950008);
		// contractActomyo_multip.push_back(0.734657907);
		// contractActomyo_multip.push_back(0.653260695);
		// contractActomyo_multip.push_back(0.667040538);
		// contractActomyo_multip.push_back(0.653581157);
		// contractActomyo_multip.push_back(0.686909149);
		// contractActomyo_multip.push_back(0.59509694 );
		// contractActomyo_multip.push_back(0.636596699);
		// contractActomyo_multip.push_back(0.56449287 );
		// contractActomyo_multip.push_back(0.655343695);
		// contractActomyo_multip.push_back(0.650536773);
		// contractActomyo_multip.push_back(0.794584201);
		// contractActomyo_multip.push_back(0.664156385);
		// contractActomyo_multip.push_back(0.692997917);
		// contractActomyo_multip.push_back(0.675212306);
		// contractActomyo_multip.push_back(0.623297548);
		// contractActomyo_multip.push_back(0.662714309);
		// contractActomyo_multip.push_back(0.536612722);
		// contractActomyo_multip.push_back(0.561288255);
		// contractActomyo_multip.push_back(0.701650377);
		// contractActomyo_multip.push_back(0.613683705);
		// contractActomyo_multip.push_back(0.647332158);
		// contractActomyo_multip.push_back(0.602467553);
		// contractActomyo_multip.push_back(0.696042301);
		// contractActomyo_multip.push_back(0.661592693);
		// contractActomyo_multip.push_back(0.697003685);
		// contractActomyo_multip.push_back(0.732574908);
		// contractActomyo_multip.push_back(0.669123538);
		// contractActomyo_multip.push_back(0.635475084);
		// contractActomyo_multip.push_back(0.62586124 );
		// contractActomyo_multip.push_back(0.632590931);
		// contractActomyo_multip.push_back(0.701490146);
		// contractActomyo_multip.push_back(0.712065374);
		// contractActomyo_multip.push_back(0.795385355);
		// contractActomyo_multip.push_back(0.798910431);
		// contractActomyo_multip.push_back(0.838647653);
		// contractActomyo_multip.push_back(0.777279282);
		// contractActomyo_multip.push_back(0.693158148);
		// contractActomyo_multip.push_back(0.647492389);
		// contractActomyo_multip.push_back(0.631950008);
		// contractActomyo_multip.push_back(0.611760936);


		// std::vector<double> contractActomyo_multip_apical;
		// contractActomyo_multip_apical.push_back(0.0);//CellID0
		// contractActomyo_multip_apical.push_back(0.0);
		// contractActomyo_multip_apical.push_back(0.563851947);
		// contractActomyo_multip_apical.push_back(0.617529242);
		// contractActomyo_multip_apical.push_back(0.626502163);
		// contractActomyo_multip_apical.push_back(0.623938471);
		// contractActomyo_multip_apical.push_back(0.623137318);
		// contractActomyo_multip_apical.push_back(0.669764461);
		// contractActomyo_multip_apical.push_back(0.63675693 );
		// contractActomyo_multip_apical.push_back(0.672808845);
		// contractActomyo_multip_apical.push_back(0.686268226);
		// contractActomyo_multip_apical.push_back(0.679858997);
		// contractActomyo_multip_apical.push_back(0.708059606);
		// contractActomyo_multip_apical.push_back(0.7002083  );
		// contractActomyo_multip_apical.push_back(0.791059125);
		// contractActomyo_multip_apical.push_back(0.752123057);
		// contractActomyo_multip_apical.push_back(0.770549591);
		// contractActomyo_multip_apical.push_back(0.756609518);
		// contractActomyo_multip_apical.push_back(0.786412434);
		// contractActomyo_multip_apical.push_back(0.76814613 );
		// contractActomyo_multip_apical.push_back(0.683544304);
		// contractActomyo_multip_apical.push_back(0.701329915);
		// contractActomyo_multip_apical.push_back(0.721038295);
		// contractActomyo_multip_apical.push_back(0.798910431);
		// contractActomyo_multip_apical.push_back(0.745874059);
		// contractActomyo_multip_apical.push_back(0.767505207);
		// contractActomyo_multip_apical.push_back(0.787534049);
		// contractActomyo_multip_apical.push_back(0.823746195);
		// contractActomyo_multip_apical.push_back(0.835603269);
		// contractActomyo_multip_apical.push_back(0.873257491);
		// contractActomyo_multip_apical.push_back(0.822945041);
		// contractActomyo_multip_apical.push_back(1          );
		// contractActomyo_multip_apical.push_back(0.834161192);
		// contractActomyo_multip_apical.push_back(0.783528281);
		// contractActomyo_multip_apical.push_back(0.830315655);
		// contractActomyo_multip_apical.push_back(0.778240667);
		// contractActomyo_multip_apical.push_back(0.795064893);
		// contractActomyo_multip_apical.push_back(0.775516744);
		// contractActomyo_multip_apical.push_back(0.756289056);
		// contractActomyo_multip_apical.push_back(0.714148374);
		// contractActomyo_multip_apical.push_back(0.692517225);
		// contractActomyo_multip_apical.push_back(0.750841211);
		// contractActomyo_multip_apical.push_back(0.785130588);
		// contractActomyo_multip_apical.push_back(0.714949527);
		// contractActomyo_multip_apical.push_back(0.675372536);
		// contractActomyo_multip_apical.push_back(0.71126422 );
		// contractActomyo_multip_apical.push_back(0.656785772);
		// contractActomyo_multip_apical.push_back(0.658227848);
		// contractActomyo_multip_apical.push_back(0.655023233);
		// contractActomyo_multip_apical.push_back(0.616567858);
		// contractActomyo_multip_apical.push_back(0.62810447 );
		// contractActomyo_multip_apical.push_back(0.629065855);
		// contractActomyo_multip_apical.push_back(0.638679699);
		// contractActomyo_multip_apical.push_back(0.660471078);
		// contractActomyo_multip_apical.push_back(0.667521231);
		// contractActomyo_multip_apical.push_back(0.568017946);
		// contractActomyo_multip_apical.push_back(0.633231854);
		// contractActomyo_multip_apical.push_back(0.564332639);
		// contractActomyo_multip_apical.push_back(0.548309566);
		// contractActomyo_multip_apical.push_back(0.584842173);
		// contractActomyo_multip_apical.push_back(0.550552796);
		// contractActomyo_multip_apical.push_back(0.539817337);
		// contractActomyo_multip_apical.push_back(0.529722801);

		for (int i = 0; i < allocPara_M.maxTotalNodeCount; i++){
			uint cellRank = i/allocPara_M.maxAllNodePerCell;
			if (cellRank >= 2 && cellRank <= 26){
			// if (cellRank >= 2 && cellRank <= 15){
				// infoVecs.contractActomyo_multip[i] = pow(contractActomyo_multip[cellRank], 1.0);////1.0;
				// infoVecs.contractActomyo_multip_apical[i] = pow(contractActomyo_multip_apical[cellRank], 1.0);//0.25;
				infoVecs.contractActomyo_multip[i] = contractActomyo_multip_perCell1;//pow(contractActomyo_multip[cellRank], 3.0);////1.0;
				infoVecs.contractActomyo_multip_apical[i] = contractActomyo_multip_perCell_apical1;//0.25;
				if (i == 2*allocPara_M.maxAllNodePerCell || i == 26*allocPara_M.maxAllNodePerCell){
					std::cout<<"infoVecs.contractActomyo_multip["<<i<<"] = "<<infoVecs.contractActomyo_multip[i]<<std::endl;
					std::cout<<"infoVecs.contractActomyo_multip_apical["<<i<<"] = "<<infoVecs.contractActomyo_multip_apical[i]<<std::endl;
				}
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = pow(contractActomyo_multip[cellRank], 1.0);//infoVecs.contractActomyo_multip[i];
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = pow(contractActomyo_multip_apical[cellRank], 1.0);//0.25;//infoVecs.contractActomyo_multip_apical[i];
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = contractActomyo_multip_perCell1;//infoVecs.contractActomyo_multip[i];
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = contractActomyo_multip_perCell_apical1;//0.25;//infoVecs.contractActomyo_multip_apical[i];
			}
			else if (cellRank >= 27 && cellRank <= 37){
			// else if (cellRank >= 16 && cellRank <= 48){
				// infoVecs.contractActomyo_multip[i] = pow(contractActomyo_multip[cellRank], 1.0);////1.0;
				// infoVecs.contractActomyo_multip_apical[i] = pow(contractActomyo_multip_apical[cellRank], 1.0);//0.25;
				infoVecs.contractActomyo_multip[i] = contractActomyo_multip_perCell2;//1.5;//pow(contractActomyo_multip[cellRank], 3.0);//0.25;//1.0;
				infoVecs.contractActomyo_multip_apical[i] = contractActomyo_multip_perCell_apical2;//1.5;
				// if (cellRank >=28 && cellRank <=34){
				// 	infoVecs.contractActomyo_multip_apical[i] = 2.0;
				// }
				if (i == 27*allocPara_M.maxAllNodePerCell || i == 37*allocPara_M.maxAllNodePerCell){
					std::cout<<"infoVecs.contractActomyo_multip["<<i<<"] = "<<infoVecs.contractActomyo_multip[i]<<std::endl;
					std::cout<<"infoVecs.contractActomyo_multip_apical["<<i<<"] = "<<infoVecs.contractActomyo_multip_apical[i]<<std::endl;
				}
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = pow(contractActomyo_multip[cellRank], 1.0);//infoVecs.contractActomyo_multip[i];
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = pow(contractActomyo_multip_apical[cellRank], 1.0);
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = contractActomyo_multip_perCell2;//1.5;//infoVecs.contractActomyo_multip[i];
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = contractActomyo_multip_perCell_apical2;//1.5;//infoVecs.contractActomyo_multip_apical[i];
				// if (cellRank >=28 && cellRank <=34){
				// 	cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = contractActomyo_multip_perCell2*1.125;
				// 	cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = contractActomyo_multip_perCell_apical2*1.0;
				// }
			}
			else if (cellRank >= 38 && cellRank <= 62){
			// else if (cellRank >= 49 && cellRank <= 62){
				// infoVecs.contractActomyo_multip[i] = pow(contractActomyo_multip[cellRank], 1.0);////1.0;
				// infoVecs.contractActomyo_multip_apical[i] = pow(contractActomyo_multip_apical[cellRank], 1.0);//0.25;
				infoVecs.contractActomyo_multip[i] = contractActomyo_multip_perCell3;//pow(contractActomyo_multip[cellRank], 3.0);//1.0;
				infoVecs.contractActomyo_multip_apical[i] = contractActomyo_multip_perCell_apical3;//0.25;
				if (i == 38*allocPara_M.maxAllNodePerCell || i == 62*allocPara_M.maxAllNodePerCell){
					std::cout<<"infoVecs.contractActomyo_multip["<<i<<"] = "<<infoVecs.contractActomyo_multip[i]<<std::endl;
					std::cout<<"infoVecs.contractActomyo_multip_apical["<<i<<"] = "<<infoVecs.contractActomyo_multip_apical[i]<<std::endl;
				}
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = pow(contractActomyo_multip[cellRank], 1.0);//infoVecs.contractActomyo_multip[i];
				// cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = pow(contractActomyo_multip_apical[cellRank], 1.0);
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = contractActomyo_multip_perCell3;//infoVecs.contractActomyo_multip[i];
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = contractActomyo_multip_perCell_apical3;//0.25;//infoVecs.contractActomyo_multip_apical[i];
			}
			else{
				infoVecs.contractActomyo_multip[i] = 0.0;
				infoVecs.contractActomyo_multip_apical[i] = 0.0;
				// std::cout<<"infoVecs.contractActomyo_multip["<<i<<"] = "<<infoVecs.contractActomyo_multip[i]<<std::endl;
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[cellRank] = infoVecs.contractActomyo_multip[i];
				cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[cellRank] = infoVecs.contractActomyo_multip_apical[i];
			}
		}
		for (int j = 0; j < allocPara_M.maxCellCount; j++){
			std::cout<<"contractActomyo_multip_perCell["<<j<<"] = "<<cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[j]<<", and multip_apical_perCell = "<<cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[j]<<std::endl;
		}

		for (int i = 0; i < allocPara_M.currentActiveCellCount; i++){
			std::cout<<"basal actomyo contractility for cell["<<i<<"] : "<<cellsSceNodes->getCellInfoVecs().contractActomyo_multip_perCell[i] <<std::endl;
			std::cout<<"apical actomyo contractility for cell["<<i<<"] : "<<cellsSceNodes->getCellInfoVecs().contractActomyo_multip_apical_perCell[i]<<std::endl;
		}
	}
	// cout << " confirm--- 2 ---" << endl;
	cout.flush();
	applySceForcesDisc_M(mitoticThreshold); // compate the MMD forces and also finds the nearset neighbor for applying the adhesion


#ifdef DebugMode
	cudaEventRecord(start3, 0);
	cudaEventSynchronize(start3);
	cudaEventElapsedTime(&elapsedTime2, start2, start3);
#endif
	// cout << " confirm--- 3 ---" << endl;
	cout.flush();
	processMembrAdh_M(timeRatio, timeRatio_Crit_Division, cycle); //applying the adhesion force 

	// cout << " confirm--- 4 ---" << endl;
	cout.flush();

	copyInterCellForces_M();//AAMIRI	

// cout << " confirm--- 5 ---" << endl;

#ifdef DebugMode
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime3, start3, stop);

	std::cout << "time spent in Node logic: " << elapsedTime1 << " "
	<< elapsedTime2 << " " << elapsedTime3 << std::endl;
#endif
}

double SceNodes::getMaxEffectiveRange() {
	int simuTypeConfigValue =
			globalConfigVars.getConfigValue("SimulationType").toInt();
	SimulationType type = parseTypeFromConfig(simuTypeConfigValue);
	if (type != Disc_M) {
		double interLinkEffectiveRange = globalConfigVars.getConfigValue(
				"InterCellLinkEffectRange").toDouble();
		double maxEffectiveRange = interLinkEffectiveRange;

		double intraLinkEffectiveRange = globalConfigVars.getConfigValue(
				"IntraCellLinkEffectRange").toDouble();
		if (intraLinkEffectiveRange > maxEffectiveRange) {
			maxEffectiveRange = intraLinkEffectiveRange;
		}

		double cartEffectiveRange = 0;
		// cartilage effective range does not apply for other types of simulation.
		try {
			cartEffectiveRange = globalConfigVars.getConfigValue(
					"CartForceEffectiveRange").toDouble();
		} catch (SceException &exce) {

		}
		if (cartEffectiveRange > maxEffectiveRange) {
			maxEffectiveRange = cartEffectiveRange;
		}
		return maxEffectiveRange;
	} else {
		double membrMembrEffRange = globalConfigVars.getConfigValue(
				"InterBEffectiveRange").toDouble();
		double membrIntnlEffRange = globalConfigVars.getConfigValue(
				"IntnlBEffectRange").toDouble();
		double intnlIntnlEffRange = globalConfigVars.getConfigValue(
				"IntraEffectRange").toDouble();
		double intnlDivEffRange = globalConfigVars.getConfigValue(
				"IntraDivEffectRange").toDouble();
		double maxEffRange = 0;
		std::vector<double> ranges;
		ranges.push_back(membrMembrEffRange);
		// all these are now
		//ranges.push_back(membrIntnlEffRange);
		//ranges.push_back(intnlIntnlEffRange);
		//ranges.push_back(intnlDivEffRange);
		maxEffRange = *std::max_element(ranges.begin(), ranges.end());
		return maxEffRange;
	}
}

void SceNodes::setInfoVecs(const NodeInfoVecs& infoVecs) {
	this->infoVecs = infoVecs;
}

void SceNodes::allocSpaceForNodes(uint maxTotalNodeCount,uint maxNumCells, uint currentActiveCellCount) {
	infoVecs.nodeLocX.resize(maxTotalNodeCount);
	infoVecs.nodeLocXHost.resize(maxTotalNodeCount); //Ali 
	infoVecs.nodeLocY.resize(maxTotalNodeCount);
	infoVecs.nodeLocYHost.resize(maxTotalNodeCount); // Ali
	infoVecs.nodeLocZ.resize(maxTotalNodeCount);
	infoVecs.nodeVelX.resize(maxTotalNodeCount);
	infoVecs.nodeVelY.resize(maxTotalNodeCount);
	infoVecs.nodeVelZ.resize(maxTotalNodeCount);
	infoVecs.basalContractPair.resize(maxTotalNodeCount,-1); // Ali

	infoVecs.nodeActomyosinMultip_basal.resize(maxTotalNodeCount);
	infoVecs.nodeActomyosinMultip_apical.resize(maxTotalNodeCount);
	infoVecs.nodeIntegrinMultip.resize(maxTotalNodeCount);
	infoVecs.quiescencePerNode.resize(maxTotalNodeCount, -100.0);

	//infoVecs.nodeContractLevel.resize(maxTotalNodeCount,0.0);// Ali
	infoVecs.nodeF_MM_C_X.resize(maxTotalNodeCount,0.0);// Ali
	infoVecs.nodeF_MM_C_Y.resize(maxTotalNodeCount,0.0);// Ali
	infoVecs.nodeContractEnergyT.resize(maxTotalNodeCount,0.0);// Ali

	infoVecs.nodeF_MI_M_x.resize(maxTotalNodeCount);  //Ali
	infoVecs.nodeF_MI_M_y.resize(maxTotalNodeCount);  //Ali
	infoVecs.nodeF_MI_M_T.resize(maxTotalNodeCount);  //Ali
	infoVecs.nodeF_MI_M_N.resize(maxTotalNodeCount,0.0);  //Ali
	infoVecs.nodeVelTangent.resize(maxTotalNodeCount);//AAMIRI
	infoVecs.nodeVelNormal.resize(maxTotalNodeCount);//AAMIRI
	infoVecs.nodeCurvature.resize(maxTotalNodeCount, 0.0);//AAMIRI
	infoVecs.nodeActinLevel.resize(maxTotalNodeCount, 0.0);//Ali
	infoVecs.nodeExtForceX.resize(maxTotalNodeCount,0.0);//AAMIRI-Ali
	infoVecs.nodeExtForceY.resize(maxTotalNodeCount,0.0);//AAMIRI-Ali
	infoVecs.nodeInterCellForceX.resize(maxTotalNodeCount,0.0);//AAMIRI-Ali
	infoVecs.nodeInterCellForceY.resize(maxTotalNodeCount,0.0);//AAMIRI-Ali
	infoVecs.nodeInterCellForceTangent.resize(maxTotalNodeCount);//AAMIRI
	infoVecs.nodeInterCellForceNormal.resize(maxTotalNodeCount);//AAMIRI
	infoVecs.nodeMaxForce.resize(maxTotalNodeCount);
	//infoVecs.nodeIsBasalMem.resize(maxTotalNodeCount,false); //Ali
	//infoVecs.nodeIsLateralMem.resize(maxTotalNodeCount,false); //Ali
	infoVecs.nodeIsApicalMem.resize(maxTotalNodeCount,0); //Ali
	infoVecs.nodeIsBasalMem.resize(maxTotalNodeCount,0); //Ali
	//infoVecs.nodeIsLateralMemHost.resize(maxTotalNodeCount,false); //Ali
	infoVecs.nodeCellType.resize(maxTotalNodeCount);
	infoVecs.nodeCellRank.resize(maxTotalNodeCount);
	infoVecs.contractActomyo_multip.resize(maxTotalNodeCount, 1.0);
	infoVecs.contractActomyo_multip_apical.resize(maxTotalNodeCount, 0.0);
	infoVecs.nodesIsEnteringMitotic.resize(maxTotalNodeCount, false);
	infoVecs.nodesContractileSpringGrowthProgress.resize(maxTotalNodeCount, 0.0);
	infoVecs.nodeIsActive.resize(maxTotalNodeCount);
	infoVecs.nodeIsActiveHost.resize(maxTotalNodeCount); // Ali
	infoVecs.nodeAdhMinDist.resize(maxTotalNodeCount); // Ali
	infoVecs.nodeCellRankFront.resize(maxNumCells,-1); // Ali
	infoVecs.nodeCellRankBehind.resize(maxNumCells,-1); // Ali
	infoVecs.nodeCellRankFrontHost.resize(maxNumCells,-1); // Ali
	infoVecs.nodeCellRankBehindHost.resize(maxNumCells,-1); // Ali
	infoVecs.nodeCellRankFrontOld.resize(maxNumCells,-1); // Ali
	infoVecs.nodeCellRankBehindOld.resize(maxNumCells,-1); // Ali
	if (controlPara.simuType == Disc
			|| controlPara.simuType == SingleCellTest) {
		infoVecs.nodeGrowPro.resize(maxTotalNodeCount);
		infoVecs.nodeInterForceX.resize(maxTotalNodeCount);
		infoVecs.nodeInterForceY.resize(maxTotalNodeCount);
		infoVecs.nodeInterForceZ.resize(maxTotalNodeCount);

	}
	if (controlPara.simuType == Disc_M) {
		infoVecs.nodeAdhereIndex.resize(maxTotalNodeCount);
		infoVecs.nodeMemMirrorIndex.resize(maxTotalNodeCount);  //Ali 
		infoVecs.nodeAdhIndxHostCopy.resize(maxTotalNodeCount);
		infoVecs.nodeAdhereIndexHost.resize(maxTotalNodeCount); //Ali 
		infoVecs.nodeMemMirrorIndexHost.resize(maxTotalNodeCount); //Ali 
		infoVecs.membrIntnlIndex.resize(maxTotalNodeCount);
		infoVecs.nodeGrowPro.resize(maxTotalNodeCount, false);
		infoVecs.nodeMitotic.resize(maxTotalNodeCount);
		infoVecs.membrTensionMag.resize(maxTotalNodeCount, 0);
		infoVecs.membrTenMagRi.resize(maxTotalNodeCount, 0);
		infoVecs.membrDistToRi.resize(maxTotalNodeCount, 0);//AAMIRI
		infoVecs.membrLinkRiMidX.resize(maxTotalNodeCount, 0);
		infoVecs.membrLinkRiMidY.resize(maxTotalNodeCount, 0);
		infoVecs.membrBendLeftX.resize(maxTotalNodeCount, 0);
		infoVecs.membrBendSpringEnergy.resize(maxTotalNodeCount, 0.0);
		infoVecs.membrLinSpringEnergy.resize(maxTotalNodeCount, 0.0);
		infoVecs.nodeIIEnergy.resize(maxTotalNodeCount, 0.0);
		infoVecs.nodeIMEnergy.resize(maxTotalNodeCount, 0.0);
		infoVecs.lagrangeFX.resize(maxTotalNodeCount, 0.0);
		infoVecs.lagrangeFY.resize(maxTotalNodeCount, 0.0);
		infoVecs.lagrangeFN.resize(maxTotalNodeCount, 0.0);
		infoVecs.membrBendLeftY.resize(maxTotalNodeCount, 0);
		infoVecs.membrBendRightX.resize(maxTotalNodeCount, 0);
		infoVecs.membrBendRightY.resize(maxTotalNodeCount, 0);
		infoVecs.dppLevel.resize(maxTotalNodeCount, 0.0); //Ali 
		infoVecs.memNodeType1.resize(maxTotalNodeCount, notAssigned1); //Ali 
		infoVecs.memNodeType1Host.resize(maxTotalNodeCount, notAssigned1); //Ali 
		infoVecs.isSubApicalJunction.resize(maxTotalNodeCount, false); //Ali

		//host vectors
	    infoVecs.nodeIsActiveH.resize(maxTotalNodeCount, 0.0);  ; //for solver
	    infoVecs.locXOldHost.resize(maxTotalNodeCount, 0.0);   ; //for solver
	    infoVecs.locYOldHost.resize(maxTotalNodeCount, 0.0);  ; //for solver
	    infoVecs.rHSXHost.resize(maxTotalNodeCount, 0.0);  ; //for solver
	    infoVecs.rHSYHost.resize(maxTotalNodeCount, 0.0);  ; //for solver
	    infoVecs.hCoefD.resize(maxTotalNodeCount, 0.0); 
		infoVecs.hCoefUd.resize(maxTotalNodeCount, 0.0); 
		infoVecs.hCoefLd.resize(maxTotalNodeCount, 0.0);  ; // for solver

	    infoVecs.isSubApicalJunctionHost.resize(maxTotalNodeCount, false); //Ali 

		auxVecs.bucketKeys.resize(maxTotalNodeCount);
		auxVecs.bucketValues.resize(maxTotalNodeCount);
		auxVecs.bucketKeysExpanded.resize(maxTotalNodeCount * 9);
		auxVecs.bucketValuesIncludingNeighbor.resize(maxTotalNodeCount * 9);
	}
	if (currentActiveCellCount != 1) {
		thrust:: sequence (infoVecs.nodeCellRankFront.begin() ,infoVecs.nodeCellRankFront.begin() +currentActiveCellCount) ; //Ali
		thrust:: sequence (infoVecs.nodeCellRankBehind.begin(),infoVecs.nodeCellRankBehind.begin()+currentActiveCellCount) ; //Ali

		thrust:: device_vector<int>  tmp1 ; 
		thrust:: device_vector<int>  tmp2 ; 

    	tmp1.resize(currentActiveCellCount,1) ; 
    	tmp2.resize(currentActiveCellCount,-1) ;

		thrust:: transform(tmp1.begin(),tmp1.begin()+currentActiveCellCount,
						   infoVecs.nodeCellRankFront.begin(),infoVecs.nodeCellRankFront.begin(), thrust::plus<int>()) ; //Ali
		thrust:: transform(tmp2.begin(),tmp2.begin()+currentActiveCellCount,
	        	           infoVecs.nodeCellRankBehind.begin(),infoVecs.nodeCellRankBehind.begin(),thrust::plus<int>()) ; //Ali
					  
		infoVecs.nodeCellRankBehind[0]=currentActiveCellCount-1 ; 
		infoVecs.nodeCellRankFront[currentActiveCellCount-1]=0 ;


	}
}
void SceNodes::initNodeAllocPara(uint totalBdryNodeCount,
		uint maxProfileNodeCount, uint maxCartNodeCount, uint maxTotalECMCount,
		uint maxNodeInECM, uint maxTotalCellCount, uint maxNodeInCell) {
	allocPara.maxCellCount = maxTotalCellCount;
	allocPara.maxNodeOfOneCell = maxNodeInCell;
	allocPara.maxNodePerECM = maxNodeInECM;
	allocPara.maxECMCount = maxTotalECMCount;
	allocPara.maxProfileNodeCount = maxProfileNodeCount;
	allocPara.maxCartNodeCount = maxCartNodeCount;

	allocPara.currentActiveProfileNodeCount = 0;
	allocPara.currentActiveCartNodeCount = 0;
	allocPara.BdryNodeCount = totalBdryNodeCount;
	allocPara.currentActiveCellCount = 0;
	allocPara.maxTotalECMNodeCount = allocPara.maxECMCount
			* allocPara.maxNodePerECM;
	allocPara.currentActiveECM = 0;

	allocPara.maxTotalCellNodeCount = maxTotalCellCount
			* allocPara.maxNodeOfOneCell;

	allocPara.startPosProfile = totalBdryNodeCount;
	allocPara.startPosCart = allocPara.startPosProfile
			+ allocPara.maxProfileNodeCount;
	allocPara.startPosECM = allocPara.startPosCart + allocPara.maxCartNodeCount;
	allocPara.startPosCells = allocPara.startPosECM
			+ allocPara.maxTotalECMNodeCount;
}

void SceNodes::initNodeAllocPara_M(uint totalBdryNodeCount,
		uint maxTotalCellCount, uint maxEpiNodePerCell,
		uint maxInternalNodePerCell) {
	allocPara_M.bdryNodeCount = totalBdryNodeCount;
	allocPara_M.currentActiveCellCount = 0;
	allocPara_M.maxCellCount = maxTotalCellCount;
	allocPara_M.maxAllNodePerCell = maxEpiNodePerCell + maxInternalNodePerCell;
	allocPara_M.maxMembrNodePerCell = maxEpiNodePerCell;
	allocPara_M.maxIntnlNodePerCell = maxInternalNodePerCell;
	allocPara_M.maxTotalNodeCount = allocPara_M.maxAllNodePerCell
			* allocPara_M.maxCellCount;
}

void SceNodes::removeNodes(int cellRank, vector<uint> &removeSeq) {
	uint cellBeginIndex = allocPara.startPosCells
			+ cellRank * allocPara.maxNodeOfOneCell;
	uint cellEndIndex = cellBeginIndex + allocPara.maxNodeOfOneCell;
	thrust::host_vector<double> cellXCoords(allocPara.maxNodeOfOneCell);
	thrust::host_vector<double> cellYCoords(allocPara.maxNodeOfOneCell);
	thrust::copy(infoVecs.nodeLocX.begin() + cellBeginIndex,
			infoVecs.nodeLocX.begin() + cellEndIndex, cellXCoords.begin());
	thrust::copy(infoVecs.nodeLocY.begin() + cellBeginIndex,
			infoVecs.nodeLocY.begin() + cellEndIndex, cellYCoords.begin());
	vector<bool> isRemove(allocPara.maxNodeOfOneCell, false);

	/*
	 std::cout << "before, X: [";
	 for (uint i = 0; i < allocPara.maxNodeOfOneCell; i++) {
	 std::cout << cellXCoords[i] << " ";
	 }
	 std::cout << "]" << endl;
	 std::cout << "before, Y: [";
	 for (uint i = 0; i < allocPara.maxNodeOfOneCell; i++) {
	 std::cout << cellYCoords[i] << " ";
	 }
	 std::cout << "]" << endl;
	 */

	for (uint i = 0; i < removeSeq.size(); i++) {
		isRemove[removeSeq[i]] = true;
	}
	thrust::host_vector<double> cellXRemoved(allocPara.maxNodeOfOneCell);
	thrust::host_vector<double> cellYRemoved(allocPara.maxNodeOfOneCell);
	uint curIndex = 0;
	for (uint i = 0; i < allocPara.maxNodeOfOneCell; i++) {
		if (isRemove[i] == false) {
			cellXRemoved[curIndex] = cellXCoords[i];
			cellYRemoved[curIndex] = cellYCoords[i];
			curIndex++;
		}
	}

	/*
	 std::cout << "after, X: [";
	 for (uint i = 0; i < allocPara.maxNodeOfOneCell; i++) {
	 std::cout << cellXRemoved[i] << " ";
	 }
	 std::cout << "]" << endl;
	 std::cout << "after, Y: [";
	 for (uint i = 0; i < allocPara.maxNodeOfOneCell; i++) {
	 std::cout << cellYRemoved[i] << " ";
	 }
	 std::cout << "]" << endl;
	 */
	thrust::copy(cellXRemoved.begin(), cellXRemoved.end(),
			infoVecs.nodeLocX.begin() + cellBeginIndex);
	thrust::copy(cellYRemoved.begin(), cellYRemoved.end(),
			infoVecs.nodeLocY.begin() + cellBeginIndex);
}

void SceNodes::processMembrAdh_M(double timeRatio, double timeRatio_Crit_Division, int cycle) {
	keepAdhIndxCopyInHost_M(timeRatio, timeRatio_Crit_Division, cycle);
	// cout << " I am before applyMembrAdh_M" << endl ; 
	applyMembrAdh_M();
	//removeInvalidPairs_M();  //Ali changed position 
}

void SceNodes::keepAdhIndxCopyInHost_M(double timeRatio, double timeRatio_Crit_Division, int cycle) {
	uint maxTotalNode = allocPara_M.currentActiveCellCount
			* allocPara_M.maxAllNodePerCell;
			// std::cout<<"maxTotalNode = "<<maxTotalNode<<std::endl;
			// std::cout<<infoVecs.nodeAdhereIndex.size()<<std::endl;
			// std::cout<<infoVecs.nodeAdhIndxHostCopy.size()<<std::endl;
			// if (cycle == 1 && timeRatio == timeRatio_Crit_Division){
			// 	for (int i = 0 ; i < infoVecs.nodeAdhIndxHostCopy.size(); i++){
			// 		std::cout<<"nodeAdhIndxHostCopy at 2nd division = "<<infoVecs.nodeAdhIndxHostCopy[i]<<std::endl;
			// 	}
			// }
	// thrust::copy(infoVecs.nodeAdhereIndex.begin(),
	// 		infoVecs.nodeAdhereIndex.begin() + maxTotalNode,
	// 		infoVecs.nodeAdhIndxHostCopy.begin());
}

void SceNodes::removeInvalidPairs_M() {
	int* nodeAdhIdxAddress = thrust::raw_pointer_cast(
			&infoVecs.nodeAdhereIndex[0]);
	uint curActiveNodeCt = allocPara_M.currentActiveCellCount
			* allocPara_M.maxAllNodePerCell;
	thrust::counting_iterator<int> iBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							infoVecs.nodeAdhereIndex.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							infoVecs.nodeAdhereIndex.begin()))
					+ curActiveNodeCt, infoVecs.nodeAdhereIndex.begin(),
			AdjustAdh(nodeAdhIdxAddress));
}

void SceNodes::applyMembrAdh_M() {
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iBegin2(0);
	uint maxTotalNode = allocPara_M.currentActiveCellCount
			* allocPara_M.maxAllNodePerCell;
	double* nodeLocXAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocY[0]);

	double* nodeGrowProAddr = thrust::raw_pointer_cast(&infoVecs.nodeGrowPro[0]);
	int* nodeAdhAddr = thrust::raw_pointer_cast(&infoVecs.nodeAdhereIndex[0]);
	double* nodedppLevelAddr = thrust::raw_pointer_cast(&infoVecs.dppLevel[0]);
	//thrust::counting_iterator<uint> iBegin_node(0); 
    // cout << " Apical adhesion is " << isApicalAdhPresent << endl ;
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeIsActive.begin(),
							infoVecs.nodeAdhereIndex.begin(), iBegin,
							infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin(),
							infoVecs.memNodeType1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeIsActive.begin(),
							infoVecs.nodeAdhereIndex.begin(), iBegin,
							infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin(),
							infoVecs.memNodeType1.begin())) + maxTotalNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin())),
			ApplyAdh(nodeLocXAddress, nodeLocYAddress, nodeGrowProAddr,nodeAdhAddr,nodedppLevelAddr,isApicalAdhPresent));
		
		//for (int i=0 ; i<140 ; i++){
		//	cout <<"adhesion index for "<<i << " is "<<infoVecs.nodeAdhereIndex[i]<< endl ; 
//		}

/*	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeIsActive.begin(),
							iBegin2,
							infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeIsActive.begin(),
							iBegin2,
							infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin())) + maxTotalNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeVelX.begin(),
									   infoVecs.nodeVelY.begin())),
			ApplyAdhReaction(nodeLocXAddress, nodeLocYAddress, nodeGrowProAddr,nodeAdhAddr,maxTotalNode));
*/
}

//AAMIRI-Ali 
//up to here nodeVelX and nodeVelY are the summation of Adhesion force between cells and Membrane-membrane forces between differen cells
void SceNodes::copyInterCellForces_M(){

	thrust::copy(infoVecs.nodeVelX.begin(), infoVecs.nodeVelX.end(),
			infoVecs.nodeInterCellForceX.begin());

	thrust::copy(infoVecs.nodeVelY.begin(), infoVecs.nodeVelY.end(),
			infoVecs.nodeInterCellForceY.begin());

}
