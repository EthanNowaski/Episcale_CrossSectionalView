/**
 * @file SimulationDomainGPU.cu
 * @brief this file contains domain level logic.
 * @author Wenzhao Sun wsun2@nd.edu Ali Nematbakhsh nematba@ucr.edu
 * @bug no know bugs
 */

#include "SimulationDomainGPU.h"

using namespace std;

#define DebugModeDomain

/**
 * Constructor.
 * reads values from config file.
 */
SimulationDomainGPU::SimulationDomainGPU() {
	readAllParameters();
}

void SimulationDomainGPU::initializeNodes_M(std::vector<SceNodeType> &nodeTypes,
		std::vector<bool> &nodeIsActive, std::vector<CVector> &initNodesVec, std::vector<CVector> &initNodeMultip_actomyo, std::vector<CVector> &initNodeMultip_integrin,
		std::vector<uint> &initActiveMembrNodeCounts,
		std::vector<uint> &initActiveIntnlNodeCounts,
		std::vector<double> &initGrowProgVec, 
		std::vector<ECellType> & eCellTypeV1 ,
		std::vector<double> & mDppV,
		std::vector<MembraneType1> & mTypeV,double InitTimeStage) {  //Ali
	/*
	 * Initialize SceNodes by constructor. first two parameters come from input parameters
	 * while the last four parameters come from Config file.
	 */
	std::cout << "Initializing nodes ...... " << std::endl;
	uint initTmpNumActiveCells ; //Ali  
    initTmpNumActiveCells=initActiveMembrNodeCounts.size() ;   //Ali size of this vector is the initial number of active cells 
	nodes = SceNodes(memPara.maxCellInDomain, memPara.maxAllNodePerCell,initTmpNumActiveCells);  // Ali // this function includes giving initial size to GPU vectors for node values 
	//nodes = SceNodes(memPara.maxCellInDomain, memPara.maxAllNodePerCell);  // Ali 

	// array size of cell type array
	uint nodeTypeSize = nodeTypes.size();
	// array size of initial active node count of cells array.
	uint initMembrNodeCountSize = initActiveMembrNodeCounts.size();
	uint initIntnlNodeCountSize = initActiveIntnlNodeCounts.size();
	// two sizes must match.
	assert(initMembrNodeCountSize == initIntnlNodeCountSize);
	assert(
			memPara.maxCellInDomain * memPara.maxAllNodePerCell
					== nodeTypes.size());

	/*
	 * second part: actual initialization
	 * copy data from main system memory to GPU memory
	 */
	NodeAllocPara_M para = nodes.getAllocParaM();
	para.currentActiveCellCount = initMembrNodeCountSize;
	assert(
			initNodesVec.size() / para.maxAllNodePerCell
					== initMembrNodeCountSize);
	nodes.setAllocParaM(para);

	cout << " I am above initValues_M " << endl ; 
	nodes.initValues_M(nodeIsActive, initNodesVec, initNodeMultip_actomyo, initNodeMultip_integrin, nodeTypes, mDppV,mTypeV);  // it copies the infomration of nodes such as locations from CPU to GPU
	cout << " I paased initValues_M " << endl ; 

	double simulationTotalTime =
			globalConfigVars.getConfigValue("SimulationTotalTime").toDouble();

	double simulationTimeStep =
			globalConfigVars.getConfigValue("SimulationTimeStep").toDouble();
	int TotalNumOfOutputFrames =
			globalConfigVars.getConfigValue("TotalNumOfOutputFrames").toInt();
	string uniqueSymbolOutput =
			globalConfigVars.getConfigValue("UniqueSymbol").toString();

	int 	freqPlotData=int ( (simulationTotalTime-InitTimeStage)/(simulationTimeStep*TotalNumOfOutputFrames) ) ;

	eCM.Initialize(memPara.maxAllNodePerCell, memPara.maxMembrNodePerCell,memPara.maxAllNodePerCell*memPara.maxCellInDomain, freqPlotData, uniqueSymbolOutput);

	cells = SceCells(&nodes, & eCM, & solver, initActiveMembrNodeCounts,
			initActiveIntnlNodeCounts, initGrowProgVec, eCellTypeV1, InitTimeStage);  //Ali

	nodes.Initialize_SceNodes  ( &cells) ;
	eCM.Initialize_SceECM(& nodes, & cells, & solver) ; 
}


// This function is called by the main function of the code, discMain_M.cpp
void SimulationDomainGPU::initialize_v2_M(SimulationInitData_V2_M& initData, double  InitTimeStage) {   //Ali 
	std::cout << "Start initializing simulation domain ......" << std::endl;
	memPara.isStab = initData.isStab;
	initializeNodes_M(initData.nodeTypes, initData.initIsActive,
			initData.initNodeVec,initData.initNodeMultip_actomyo, initData.initNodeMultip_integrin, initData.initActiveMembrNodeCounts,
			initData.initActiveIntnlNodeCounts, initData.initGrowProgVec, 
			initData.eCellTypeV1,initData.mDppV,initData.mTypeV, InitTimeStage);  // Ali
	std::cout << "Finished initializing nodes positions" << std::endl;
	nodes.initDimension(domainPara.minX, domainPara.maxX, domainPara.minY,
			domainPara.maxY, domainPara.gridSpacing);
	//std::cout << "finished init nodes dimension" << std::endl;
	// The domain task is not stabilization unless specified in the next steps.
	stabPara.isProcessStab = false;
	std::cout << "Finished initializing simulation domain" << std::endl;
}

/**
 * Highest level logic of domain.
 *
 */
void SimulationDomainGPU::runAllLogic(double dt) {

	if (memPara.simuType == Disc) {
		nodes.sceForcesDisc();
	}

	// This function applies velocity so nodes actually move inside this function.
	if (memPara.simuType == Disc) {
		cells.runAllCellLevelLogicsDisc(dt);
	}

	if (memPara.simuType == SingleCellTest) {
		nodes.sceForcesDisc();
		cells.runStretchTest(dt);
	}
}

void SimulationDomainGPU::runAllLogic_M(double & dt, double Damp_Coef, double InitTimeStage, 
											double timeRatio, double timeRatio_Crit_actomyo, double timeRatio_Crit_ECM, double timeRatio_Crit_Division,
												double volume_Increase_Target_Ratio, double volume_Increase_Scale, double postDivision_restorationRateScale, int cycle,
												double distFromNucleus_max, double distFromNucleus_min, double distFromNucleus_normalMax1, double distFromNucleus_normalMax2, double distFromNucleus_normalMax3, 
												double distFromNucleus_normalMax_apical1, double distFromNucleus_normalMax_apical2, double distFromNucleus_normalMax_apical3, 
												double percentage_before_timeRatio_Crit_Division_scaling,
												double growthProgressSpeed, int maxApicalBasalNodeNum, double maxLengthToAddMemNodes, double mitoRndActomyoStrengthScaling, double thresholdToIntroduceNewCell, 
												double contractActomyo_multip_perCell1, double contractActomyo_multip_perCell2, double contractActomyo_multip_perCell3,
												double contractActomyo_multip_perCell_apical1, double contractActomyo_multip_perCell_apical2, double contractActomyo_multip_perCell_apical3, double mitoticThreshold) {   

#ifdef DebugModeDomain
	cudaEvent_t start1, start2, stop;
	float elapsedTime1, elapsedTime2;
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&stop);
	cudaEventRecord(start1, 0);
#endif
	// cout << "--- 1 ---" << endl;
	cout.flush();
	
	nodes.sceForcesDisc_M(timeRatio, timeRatio_Crit_Division, cycle,
								contractActomyo_multip_perCell1, contractActomyo_multip_perCell2, contractActomyo_multip_perCell3,
								contractActomyo_multip_perCell_apical1, contractActomyo_multip_perCell_apical2, contractActomyo_multip_perCell_apical3,
								mitoticThreshold);
	// cout << "--- 2 ---" << endl;
	cout.flush();
#ifdef DebugModeDomain
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif
	// cout << "--- 3 ---" << endl;
	cout.flush();

	cells.runAllCellLogicsDisc_M(dt,Damp_Coef,InitTimeStage, timeRatio, timeRatio_Crit_actomyo, timeRatio_Crit_ECM, timeRatio_Crit_Division, volume_Increase_Target_Ratio, volume_Increase_Scale, postDivision_restorationRateScale, cycle,
									distFromNucleus_max, distFromNucleus_min, distFromNucleus_normalMax1, distFromNucleus_normalMax2,distFromNucleus_normalMax3,
									distFromNucleus_normalMax_apical1,distFromNucleus_normalMax_apical2,distFromNucleus_normalMax_apical3, 
									percentage_before_timeRatio_Crit_Division_scaling, growthProgressSpeed, maxApicalBasalNodeNum,
									 maxLengthToAddMemNodes, mitoRndActomyoStrengthScaling, thresholdToIntroduceNewCell, mitoticThreshold);
	// cout << "--- 4 ---" << endl;
	cout.flush();
#ifdef DebugModeDomain
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime2, start2, stop);
	//std::cout << "time spent in Simu Domain logic: " << elapsedTime1 << " "
	//<< elapsedTime2 << std::endl;
#endif
}

void SimulationDomainGPU::readMemPara() {
	int simuTypeConfigValue =
			globalConfigVars.getConfigValue("SimulationType").toInt();

	memPara.simuType = parseTypeFromConfig(simuTypeConfigValue);

	memPara.maxCellInDomain =
			globalConfigVars.getConfigValue("MaxCellInDomain").toInt();
	if (memPara.simuType != Disc_M) {
		memPara.maxNodePerCell = globalConfigVars.getConfigValue(
				"MaxNodePerCell").toInt();
	}

	memPara.maxECMInDomain = 0;
	memPara.maxNodePerECM = 0;
	memPara.FinalToInitProfileNodeCountRatio = 0;

	if (memPara.simuType == Disc_M) {
		memPara.maxMembrNodePerCell = globalConfigVars.getConfigValue(
				"MaxMembrNodeCountPerCell").toInt();
		memPara.maxIntnlNodePerCell = globalConfigVars.getConfigValue(
				"MaxIntnlNodeCountPerCell").toInt();
		memPara.maxAllNodePerCell = globalConfigVars.getConfigValue(
				"MaxAllNodeCountPerCell").toInt();
		assert(
				memPara.maxMembrNodePerCell + memPara.maxIntnlNodePerCell
						== memPara.maxAllNodePerCell);
	}

     cout << "I am in readMemPara"<< endl ; 
}

void SimulationDomainGPU::readDomainPara() {
	domainPara.minX = globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble();
	domainPara.maxX = globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble();
	domainPara.minY = globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble();
	domainPara.maxY = globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble();
	//domainPara.minZ = globalConfigVars.getConfigValue("DOMAIN_ZMIN").toDouble();
	//domainPara.maxZ = globalConfigVars.getConfigValue("DOMAIN_ZMAX").toDouble();
	domainPara.gridSpacing = nodes.getMaxEffectiveRange();
	domainPara.XBucketSize = (domainPara.maxX - domainPara.minX)
			/ domainPara.gridSpacing + 1;
	domainPara.YBucketSize = (domainPara.maxY - domainPara.minY)
			/ domainPara.gridSpacing + 1;

     cout << "I am in readDomainPara"<< endl ; 
}

void SimulationDomainGPU::readAllParameters() {
	readMemPara();
	readDomainPara();
cout<< "I am in SimulationDomainGPU constructor" << endl ; 
}

void SimulationDomainGPU::outputVtkFilesWithCri(std::string scriptNameBase,
		int rank, AnimationCriteria aniCri) {
	nodes.prepareSceForceComputation();
	VtkAnimationData aniData = nodes.obtainAnimationData(aniCri);
	aniData.outputVtkAni(scriptNameBase, rank);
}

void SimulationDomainGPU::outputVtkFilesWithCri_M(std::string scriptNameBase,
		int rank, AnimationCriteria aniCri) {
	nodes.prepareSceForceComputation();
	//std::cout << "started generate raw data" << std::endl;
	AniRawData rawAni = cells.obtainAniRawData(aniCri);
	//std::cout << "finished generate raw data" << std::endl;
	VtkAnimationData aniData = cells.outputVtkData(rawAni, aniCri);
	//std::cout << "finished generate vtk data" << std::endl;
	aniData.outputVtkAni(scriptNameBase, rank);
	//std::cout << "finished generate vtk file" << std::endl;
}

void SimulationDomainGPU::outputVtkGivenCellColor(std::string scriptNameBase,
		int rank, AnimationCriteria aniCri, std::vector<double>& cellColorVal,std::vector<double>& cellsPerimeter) {
	nodes.prepareSceForceComputation();
	AniRawData rawAni = cells.obtainAniRawDataGivenCellColor(cellColorVal,
			aniCri,cellsPerimeter); //AliE
	VtkAnimationData aniData = cells.outputVtkData(rawAni, aniCri);
	aniData.outputVtkAni(scriptNameBase, rank);
}

void SimulationDomainGPU::outputVtkColorByCell_T1(std::string scriptNameBase,
		int rank, AnimationCriteria aniCri) {
	assert(aniCri.animationType == T1Tran);
	std::vector<double> t1ColorVec = processT1Color();
	//outputVtkGivenCellColor(scriptNameBase, rank, aniCri, t1ColorVec);  //Ali comment node uing now so no need for modification
}

std::vector<double> SimulationDomainGPU::processT1Color() {
	std::vector<double> result;
	result.resize(cells.getAllocParaM().currentActiveCellCount);
	for (int i = 0; i < int(result.size()); i++) {
		if (t1CellSet.find(i) == t1CellSet.end()) {
			result[i] = 0;
		} else {
			result[i] = 1;
		}
	}
	return result;
}

void SimulationDomainGPU::outputVtkColorByCell_polySide(
		std::string scriptNameBase, int rank, AnimationCriteria aniCri) {

        std::cout  <<"I am in modified function" <<std:: endl; 
	assert(aniCri.animationType == PolySide);
        std:: vector<double> cellsPerimeter ; 
	std::vector<double> polySideColorVec = processPolySideColor(cellsPerimeter);
	outputVtkGivenCellColor(scriptNameBase, rank, aniCri, polySideColorVec,cellsPerimeter);
        cellsPerimeter.clear(); 
}
void SimulationDomainGPU::outputResumeData(uint frame) {
	WriteResumeData writeResumeData ; 
    std::cout  <<"I am writing Resume Data file" <<std:: endl;

	//Gather information
	std::vector<AniResumeData> aniResumeDatas= cells.obtainResumeData() ;
	aniResumeDatas.push_back                    (eCM.obtainResumeData()); 
    
	//Fetch the input parameters
	std::string uniqueSymbol        = globalConfigVars.getConfigValue(
	                                  "UniqueSymbol").toString() ; 
	std::string membFileNameResume  = globalConfigVars.getConfigValue(
	                        		  "MembraneNodes_FileName_Resume").toString() ;
	std::string intnlFileNameResume = globalConfigVars.getConfigValue(
			                          "IntnlNodes_FileName_Resume").toString() ;
	
	//Write the gathered information 
	//0 is for membrane nodes and 1 is for internal node, 2 for cells, 	3 is for ECM nodes
	writeResumeData.writeForMembAndIntnl(aniResumeDatas.at(0),aniResumeDatas.at(1), membFileNameResume, intnlFileNameResume, uniqueSymbol) ;  
	writeResumeData.writeForCells       (aniResumeDatas.at(2), uniqueSymbol) ; 
	writeResumeData.writeForECM         (aniResumeDatas.at(3), uniqueSymbol) ; 
}

std::vector<double> SimulationDomainGPU::processPolySideColor(std:: vector<double> & cellsPerimeter) {
	CellsStatsData cellStatsVec = cells.outputPolyCountData();
	
        for (int i=0; i< int (cellStatsVec.cellsStats.size());  i++) {
        //std::cout << cellStatsVec.cellsStats.size() <<"cellsStatsvec vector size is" <<std:: endl; 
       	   cellsPerimeter.push_back(cellStatsVec.cellsStats[i].cellPerim) ;
        }

         //AliE        
	std::vector<double> result = cellStatsVec.outputPolySides(); 
	return result;
}

vector<vector<int> > SimulationDomainGPU::outputLabelMatrix(
		std::string resultNameBase, int rank, PixelizePara& pixelPara) {
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string resultNameRank = ss.str();
	std::string matrixFileName = resultNameBase + resultNameRank + ".dat";
	vector<vector<int> > matrix = nodes.obtainLabelMatrix(pixelPara);
	printMatrixToFile(matrix, matrixFileName);
	return matrix;
}

void SimulationDomainGPU::outputGrowthProgressAuxFile(int step) {
	static bool isFirstTime = true;
	std::string auxDataFileName = globalConfigVars.getConfigValue(
			"DataOutputFolder").toString()
			+ globalConfigVars.getConfigValue("GrowthAuxFileName").toString();
	if (isFirstTime) {
		std::remove(auxDataFileName.c_str());
		isFirstTime = false;
	}
	std::cout << "Updating growth progress file" << std::endl;
	ofstream ofs;
	ofs.open(auxDataFileName.c_str(), ios::app);
	ofs << step << " ";
	std::vector<double> growProVec = cells.getGrowthProgressVec();
	for (std::vector<double>::iterator it = growProVec.begin();
			it != growProVec.end(); ++it) {
		ofs << *it << " ";
	}
	ofs << std::endl;
	ofs.close();
}

void SimulationDomainGPU::analyzeLabelMatrix(vector<vector<int> > &labelMatrix,
		int step, std::string &imageFileNameBase, std::string &statFileName) {
	ResAnalysisHelper resHelper;

	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << step;
	std::string imgNameRank = ss.str();
	std::string imgFileName = imageFileNameBase + imgNameRank + ".bmp";

	resHelper.outputImg_formatBMP(imgFileName, labelMatrix);
	std::vector<double> growthProVec = cells.getGrowthProgressVec();
	if (memPara.simuType == Disc) {
		resHelper.outputStat_PolygonCounting(statFileName, step, labelMatrix,
				growthProVec);
		outputGrowthProgressAuxFile(step);
	} else {
		resHelper.outputStat_PolygonCounting(statFileName, step, labelMatrix);
	}
}

bool SimulationDomainGPU::isDividing_ForAni() {
	if (cells.aniDebug) {
		cells.aniDebug = false;
		return true;
	}
	return false;
}

void SimulationDomainGPU::performAblation(AblationEvent& ablEvent) {
	thrust::host_vector<double> xCoord = nodes.getInfoVecs().nodeLocX;
	thrust::host_vector<double> yCoord = nodes.getInfoVecs().nodeLocY;

	AblationEvent aa;

	for (uint i = 0; i < xCoord.size(); i++) {
		double xDiff = xCoord[i] - 25.3;
		double yDiff = yCoord[i] - 25.2;
		if (xDiff * xDiff + yDiff * yDiff < 0.04) {
			uint cellRank = i / 90;
			uint nodeRank = i % 90;
			std::cout << "cell : " << cellRank << ", node: " << nodeRank
					<< "pos: (" << xCoord[i] << "," << yCoord[i] << ")"
					<< std::endl;
			bool found = false;
			for (uint j = 0; j < aa.ablationCells.size(); j++) {
				if (aa.ablationCells[j].cellNum == cellRank) {
					found = true;
					aa.ablationCells[j].nodeNums.push_back(nodeRank);
				}
			}
			if (!found) {
				AblaInfo cellNew;
				cellNew.cellNum = cellRank;
				cellNew.nodeNums.push_back(nodeRank);
				aa.ablationCells.push_back(cellNew);
			}
		}
	}

	aa.printInfo();
	int jj;
	cin >> jj;

	cells.runAblationTest(aa);
}

//We may be able to remove this function by repositioning the location of calling output data
CellsStatsData SimulationDomainGPU::outputPolyCountData() {
	// this step is necessary for obtaining correct neighbors because new cells might have been created in previous step.
	//nodes.sceForcesDisc_M(); // Ali commented this. it will intefere with logics of the code.
	return cells.outputPolyCountData();
}
SingleCellData  SimulationDomainGPU::OutputStressStrain() {
	// this step is necessary for obtaining correct neighbors because new cells might have been created in previous step.
	//nodes.sceForcesDisc_M(); // Ali commented this. it will intefere with logics of the code.
	return cells.OutputStressStrain();
}


NetworkInfo SimulationDomainGPU::buildNetInfo(CellsStatsData &polyData) {
	std::vector<NetworkNode> netNodes;
	for (uint i = 0; i < polyData.cellsStats.size(); i++) {
		NetworkNode netNode;
		netNode.setGrowP(polyData.cellsStats[i].cellGrowthProgress);
		netNode.setNodeRank(polyData.cellsStats[i].cellRank);
		netNode.setPos(polyData.cellsStats[i].cellCenter);

		std::vector<int> ngbrVec;
		std::set<int>::iterator it;
		for (it = polyData.cellsStats[i].neighborVec.begin();
				it != polyData.cellsStats[i].neighborVec.end(); ++it) {
			ngbrVec.push_back(*it);
		}

		netNode.setNgbrList(ngbrVec);
		netNodes.push_back(netNode);
	}
	NetworkInfo result(netNodes);
	return result;
}


std::set<int> SimulationDomainGPU::findT1Transition() {
	std::set<int> result;
	for (uint i = 0; i < preT1Vec.size(); i++) {
		for (uint j = 0; j < preT1Vec[i].size(); j++) {
			if (netInfo.isT1Tran(preT1Vec[i][j])) {
				result.insert(preT1Vec[i][j].nodeRank);
				result.insert(preT1Vec[i][j].centerNgbr);
				result.insert(preT1Vec[i][j].sideNgbrs[0]);
				result.insert(preT1Vec[i][j].sideNgbrs[1]);
			}
		}
	}
	if (result.size() != 0) {
		std::cout << "found T1 transition!" << std::endl;
	}
	return result;
}

void SimulationDomainGPU::processT1Info(int maxStepTraceBack,
		CellsStatsData &polyData) {
	// first, construct network info
	netInfo = buildNetInfo(polyData);

	// second, find all of the previous pre-t1 states matches
	// has make t1 transition under current network info. output
	// these cell numbers.
	t1CellSet = findT1Transition();

	// finally, update the pre-T1 info vector by remove old one
	// and add new one.
	if (preT1Vec.size() >= maxStepTraceBack) {
		int eraseSize = preT1Vec.size() - maxStepTraceBack + 1;
		preT1Vec.erase(preT1Vec.begin(), preT1Vec.begin() + eraseSize);
	}
	std::vector<PreT1State> preT1States = netInfo.scanForPreT1States();
	preT1Vec.push_back(preT1States);
}

