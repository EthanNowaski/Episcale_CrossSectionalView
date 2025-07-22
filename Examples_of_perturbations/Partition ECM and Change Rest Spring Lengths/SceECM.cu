#include "SceECM.h"
#include "SceCells.h" // Because of forward declaration
//# define debugModeECM 
// bending stiffness is given inside the code. It should be given as in input from a txt file.
//isInitPhase bool variable is not active anymore.
//Right now it is assumed that ECM stiffness is the same everywhere.
__constant__ double sceInterCell_ECM[5]; 
//__constant__ double wLCPara_ECM[4]; 
__constant__ double restLenECMAdhSpringGPU  ;  
__constant__ double maxLenECMAdhSpringGPU ;
__constant__ double kAdhECMGPU ;
__constant__ double stiffnessECMBasalGPU ;
__constant__ double stiffnessECMBCGPU ;
__constant__ double stiffnessECMPeripGPU ;
__constant__ double lknotECMBasalGPU ;
__constant__ double lknotECMBCGPU ;
__constant__ double lknotECMPeripGPU ;

const double smallNumber=.000001 ; 



namespace patch{
	template <typename  T> std::string to_string (const T& n) 
	{
	std:: ostringstream stm ; 
	stm << n ; 
	return stm.str() ; 
	}
}



__device__
void DefineECMStiffnessAndLknot(EType nodeType, int nodeIdx, int* cellNeighborIdPtr, int* cellSubdomainIndxPtr, double& stiffness, double& sponLen) {

	if (nodeType==excm) {
		// stiffness=stiffnessECMBasalGPU ;   
		// sponLen=lknotECMBasalGPU  ; 
		if (cellSubdomainIndxPtr[cellNeighborIdPtr[nodeIdx]] == 0){
			stiffness=stiffnessECMBasalGPU ;   
			sponLen=lknotECMBasalGPU ;
			}
		else if (cellSubdomainIndxPtr[cellNeighborIdPtr[nodeIdx]] == 1){
			stiffness=stiffnessECMBasalGPU ;   
			sponLen=lknotECMBasalGPU * 0.5; // For ECM nodes connected to central cells the resting spring length is half
			}
		else if (cellSubdomainIndxPtr[cellNeighborIdPtr[nodeIdx]] == 2){
			stiffness=stiffnessECMBasalGPU ;   
			sponLen=lknotECMBasalGPU ;
			}
	}
	if (nodeType==perip) {
		stiffness=stiffnessECMPeripGPU ; 
		sponLen=lknotECMPeripGPU ;
	}

	if (nodeType==bc2) {
		stiffness=stiffnessECMBCGPU;  
		sponLen=lknotECMBCGPU ; 
	}

}

__device__
double calMorse_ECM(const double& linkLength ) {

	double forceValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterCell_ECM[0] / sceInterCell_ECM[2]
				* exp(-linkLength / sceInterCell_ECM[2])
				+ sceInterCell_ECM[1] / sceInterCell_ECM[3]
						* exp(-linkLength / sceInterCell_ECM[3]);
//		if (forceValue > 0) {
//			forceValue = 0;
//		}
	}

	return (forceValue) ; 
}

__device__
double calMorse_ECM_mitotic(const double& linkLength, double scaling ) {

	double forceValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		forceValue = 0;
	} else {
		forceValue = (1.0+(3.0-1.0)*scaling)*(-sceInterCell_ECM[0] / sceInterCell_ECM[2]
				* exp(-linkLength / sceInterCell_ECM[2])
				+ sceInterCell_ECM[1] / sceInterCell_ECM[3]
						* exp(-linkLength / sceInterCell_ECM[3]));
		// forceValue = -(1.0+(2.0-1.0)*scaling)*sceInterCell_ECM[0] / sceInterCell_ECM[2]
		// 		* exp(-linkLength / sceInterCell_ECM[2])
		// 		+ (1.0+(2.0-1.0)*scaling)*sceInterCell_ECM[1] / sceInterCell_ECM[3]
		// 				* exp(-linkLength / sceInterCell_ECM[3]);
//		if (forceValue > 0) {
//			forceValue = 0;
//		}
	}

	return (forceValue) ; 
}

__device__
double calMorseEnergy_ECM(const double& linkLength ) {

	double energyValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		energyValue = 0;
	} else {
		energyValue = sceInterCell_ECM[0]* exp(-linkLength / sceInterCell_ECM[2])
				    - sceInterCell_ECM[1]* exp(-linkLength / sceInterCell_ECM[3]);
	}

	return (energyValue) ; 
}

/*
__device__
double calWLC_ECM(const double& linkLength ) {

	double x=linkLength/wLCPara_ECM[0] ;
	return (wLCPara_ECM[1]*( 6*x+ ( x*x*(3.0-2*x))/( (1-x)*(1-x) ) )
	       -wLCPara_ECM[2]/pow(linkLength,wLCPara_ECM[3]) ) ; 	
}
*/
__device__
bool IsValidAdhPair(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU  && dist < maxLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}
__device__
bool IsValidAdhPairForNotInitPhase(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}


__device__
double  CalAdhECM(const double& dist ) {
	return (kAdhECMGPU*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}

__device__
double  CalAdhEnergy(const double& dist ) {
	return (0.5*kAdhECMGPU*(dist-restLenECMAdhSpringGPU)*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}


EType SceECM:: ConvertStringToEType(string eNodeRead) {
	if (eNodeRead=="perip")  {
		return perip ; 
	}
	else if (eNodeRead=="bc2") {
		return bc2 ; 
	}
	else if (eNodeRead=="excm") {
		return excm ; 
	}
	else {
		cout << "Error in defining type of external nodes" << endl ; 
		return excm ;// To just return something to avoid compiler complain 
	}
} 



SceECM::SceECM() {

	isECMNeighborSet=false ; 
	eCMRemoved=false ; 
	isECMNeighborResetPostDivision = false;
}



void SceECM::Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM, int freqPlotData, string uniqueSymbol) {

	maxAllNodePerCell=maxAllNodePerCellECM ; 
	maxMembrNodePerCell= maxMembrNodePerCellECM ; 
    // maxTotalNodes=maxTotalNodesECM ; //Ali 
    this->freqPlotData=freqPlotData ; 
    this->uniqueSymbol=uniqueSymbol ; 

	std::fstream readCoord_ECM ;
	std::fstream readInput_ECM ;
	int numberNodes_ECM ; 
	double tmpPosX_ECM,tmpPosY_ECM ; 
	vector<double> posXIni_ECM,posYIni_ECM ;
	vector <EType> eNodeVec ; 


	int resumeSimulation = globalConfigVars.getConfigValue(
	"ResumeSimulation").toInt();
	
	if (resumeSimulation==0) { 
		cout << " In the ECM module, I am in start mode" << endl ; 
		readCoord_ECM.open("./resources/coordinate_ECM21.txt") ;
	}
	else if(resumeSimulation==1) { 
		cout << " In the ECM module, I am in resume mode" << endl ;
		std::string secondInputFileName = "./resources/DataFileECM_" + uniqueSymbol + "Resume.cfg";
		readCoord_ECM.open(secondInputFileName.c_str()) ;
	}
	else{

		throw std::invalid_argument(" ResumeSimulation parameter in the input file must be either 1 or 0. Error from ECM module"); 
	}
 
if (readCoord_ECM.is_open()) {
	cout << "ECM coordinates file opened successfully" <<endl ; 
}
else {
	cout << "ECM coordinates file is not opened successfully" << endl ; 
}


string inputInfoText ; 
string eNodeRead ; 
readCoord_ECM>>numberNodes_ECM ;
for (int i=0 ; i<numberNodes_ECM ; i++){
	readCoord_ECM>>tmpPosX_ECM>>tmpPosY_ECM>>eNodeRead  ; 
	posXIni_ECM.push_back(tmpPosX_ECM) ; 
	posYIni_ECM.push_back(tmpPosY_ECM) ; 
	EType eNode=ConvertStringToEType(eNodeRead) ; 
	eNodeVec.push_back(eNode) ; 
}

readInput_ECM.open("./resources/ECM_input.txt") ;
if (readInput_ECM.is_open()) {
	cout << "ECM Mech input opened successfully" <<endl ; 
}
else {
	cout << "ECM Mech input is not opened successfully" << endl ; 
}

 
 readInput_ECM>> inputInfoText ; 
 for (int i=0 ; i<5; i++) {
 	readInput_ECM>> mechPara_ECM.sceInterCellCPU_ECM[i] ; //=39.0 ; 
 }
 
// readInput_ECM>>restLenECMSpring ;
// readInput_ECM>>eCMLinSpringStiff ;    
 readInput_ECM>>restLenECMAdhSpring  ;  
 readInput_ECM>>maxLenECMAdhSpring ;
 readInput_ECM>>kAdhECM ;
 //for ( int i=0 ; i<4 ; i++) {
//	readInput_ECM>>mechPara_ECM.wLCParaCPU_ECM[i] ;
// }    



std::fstream secondInput_ECM ; 
std:: string secondInputInfo ;  //dummy 
std::string secondInputFileName = "./resources/ECM_" + uniqueSymbol + "input.cfg";
secondInput_ECM.open(secondInputFileName.c_str()) ;
//secondInput_ECM.open("./resources/ECM_N01G00_input.cfg" ) ;
if (secondInput_ECM.is_open()) {
	cout << "Second ECM Mech input opened successfully" <<endl ; 
}
else {
	cout << "Second ECM Mech input is not opened successfully" << endl ; 
}

 secondInput_ECM>>secondInputInfo ;  // just for information no use in the code
 secondInput_ECM>>stiffnessECMBasal ;
 secondInput_ECM>>stiffnessECMBC ;
 secondInput_ECM>>stiffnessECMPerip ;
 secondInput_ECM>>lknotECMBasal ;
 secondInput_ECM>>lknotECMBC ;
 secondInput_ECM>>lknotECMPerip ;
 secondInput_ECM>>dampBasal ;
 secondInput_ECM>>dampBC ;
 secondInput_ECM>>dampApical ;

 cout <<" stiffness of ECM at the basal side is="<<stiffnessECMBasal <<endl ;   
 cout <<" stiffness of ECM at boundary is="<<stiffnessECMBC<<endl ; 
 cout <<" stiffness of ECM peripodial side is="<<stiffnessECMPerip<<endl ; 

 cout <<" rest len basal ECM is="<<lknotECMBasal<<endl ;
 cout <<" rest len boundary ECM is= "<<lknotECMBC<<endl ;
 cout << "rest len peripodial ECM is=" <<lknotECMPerip <<endl ; 

 cout << "Damping for basal ECM is="<<dampBasal<<endl ;
 cout << "Damping for boundary ECM is= "<<dampBC<<endl ;
 cout << "Damping for peripodial ECM is=" <<dampApical <<endl ; 

 cout << "number of ECM nodes is"<< numberNodes_ECM <<endl ; 
 maxTotalNodes = numberNodes_ECM*2;
 cout << "max number of ECM nodes is"<< maxTotalNodes <<endl ; 




 for (int i=0 ; i<5; i++) {
	cout <<"Morse parameter number"<<i<<" is " <<mechPara_ECM.sceInterCellCPU_ECM[i]<<endl ; 

} 

 //cout <<"rest length of ECM spring is "<<restLenECMSpring<<endl ;   

// cout <<"ECM spring stiffness is "<<eCMLinSpringStiff<<endl ; 

 cout <<"ECM Membrane neutral adhesion length is "<<restLenECMAdhSpring<<endl ;  
 cout <<"ECM Membrane max adhesion length is "<<maxLenECMAdhSpring<<endl ;
 cout <<"ECM Membrane adhesion stiffness is "<<kAdhECM<<endl ;
 cout << "ECM only applies adhesvie force" << endl ; 

//for ( int i=0 ; i<4 ; i++) {
//	cout<<"wLC parameter "<< i << " is "<<mechPara_ECM.wLCParaCPU_ECM[i]<<endl ;  ;
//}    

cudaMemcpyToSymbol(sceInterCell_ECM,mechPara_ECM.sceInterCellCPU_ECM
			,5*sizeof(double)); 
//cudaMemcpyToSymbol(wLCPara_ECM,mechPara_ECM.wLCParaCPU_ECM
//			,4*sizeof(double)); 

cudaMemcpyToSymbol(restLenECMAdhSpringGPU, &restLenECMAdhSpring,sizeof(double));

cudaMemcpyToSymbol(maxLenECMAdhSpringGPU, &maxLenECMAdhSpring,sizeof(double));
cudaMemcpyToSymbol(kAdhECMGPU, &kAdhECM,sizeof(double));

cudaMemcpyToSymbol(stiffnessECMPeripGPU, &stiffnessECMPerip,sizeof(double));
cudaMemcpyToSymbol(stiffnessECMBCGPU, &stiffnessECMBC,sizeof(double));
cudaMemcpyToSymbol(stiffnessECMBasalGPU, &stiffnessECMBasal,sizeof(double));

cudaMemcpyToSymbol(lknotECMPeripGPU, & lknotECMPerip,sizeof(double));
cudaMemcpyToSymbol(lknotECMBCGPU, & lknotECMBC,sizeof(double));
cudaMemcpyToSymbol(lknotECMBasalGPU, & lknotECMBasal,sizeof(double));


counter=100000 ; //large number
lastPrintECM=1000000 ; // large number 
outputFrameECM=0 ; 
numNodesECM= numberNodes_ECM ; //(eCMMaxX-eCMMinX)/eCMMinDist ; 



// indexECM.resize(numNodesECM,0) ;
// peripORexcm.resize(numNodesECM,perip) ;
// dampCoef.resize(numNodesECM) ; 
// nodeECMLocX.resize(numNodesECM,0.0) ;
// nodeECMLocY.resize(numNodesECM,0.0) ;


// cellNeighborId.resize(numNodesECM,-1) ;

// stiffLevel.resize(numNodesECM) ;
// sponLen.resize(numNodesECM) ;

// linSpringForceECMX.resize(numNodesECM,0.0); 
// linSpringForceECMY.resize(numNodesECM,0.0); 
// linSpringAvgTension.resize(numNodesECM,0.0); 
// linSpringEnergy.resize(numNodesECM,0.0); 
// morseEnergy.resize(numNodesECM,0.0); 
// adhEnergy.resize(numNodesECM,0.0); 


// bendSpringForceECMX.resize(numNodesECM,0.0); 
// bendSpringForceECMY.resize(numNodesECM,0.0);  
// memMorseForceECMX.resize(numNodesECM,0.0); 
// memMorseForceECMY.resize(numNodesECM,0.0);
 
// fBendCenterX.resize(numNodesECM,0.0); 
// fBendCenterY.resize(numNodesECM,0.0); 
// fBendLeftX.resize(numNodesECM,0.0); 
// fBendLeftY.resize(numNodesECM,0.0); 
// fBendRightX.resize(numNodesECM,0.0); 
// fBendRightY.resize(numNodesECM,0.0); 
 
// totalForceECMX.resize(numNodesECM,0.0); 
// totalForceECMY.resize(numNodesECM,0.0);

// totalExplicitForceECMX.resize(numNodesECM,0.0); 
// totalExplicitForceECMY.resize(numNodesECM,0.0);


// rHSX.resize(numNodesECM,0.0); 
// rHSY.resize(numNodesECM,0.0);
// //memNodeType.resize(maxTotalNodes,notAssigned1) ; 


// nodeIsActive.resize(numNodesECM,true) ; 
// thrust::sequence (indexECM.begin(),indexECM.begin()+numNodesECM);

indexECM.resize(maxTotalNodes,0) ;
peripORexcm.resize(maxTotalNodes,perip) ;
dampCoef.resize(maxTotalNodes) ; 
nodeECMLocX.resize(maxTotalNodes,0.0) ;
nodeECMLocY.resize(maxTotalNodes,0.0) ;
// isActiveECM.resize(maxTotalNodes,false);


cellNeighborId.resize(maxTotalNodes,-1) ;

stiffLevel.resize(maxTotalNodes) ;
sponLen.resize(maxTotalNodes) ;

linSpringForceECMX.resize(maxTotalNodes,0.0); 
linSpringForceECMY.resize(maxTotalNodes,0.0); 
linSpringAvgTension.resize(maxTotalNodes,0.0); 
linSpringEnergy.resize(maxTotalNodes,0.0); 
morseEnergy.resize(maxTotalNodes,0.0); 
adhEnergy.resize(maxTotalNodes,0.0); 


bendSpringForceECMX.resize(maxTotalNodes,0.0); 
bendSpringForceECMY.resize(maxTotalNodes,0.0);  
memMorseForceECMX.resize(maxTotalNodes,0.0); 
memMorseForceECMY.resize(maxTotalNodes,0.0);
 
fBendCenterX.resize(maxTotalNodes,0.0); 
fBendCenterY.resize(maxTotalNodes,0.0); 
fBendLeftX.resize(maxTotalNodes,0.0); 
fBendLeftY.resize(maxTotalNodes,0.0); 
fBendRightX.resize(maxTotalNodes,0.0); 
fBendRightY.resize(maxTotalNodes,0.0); 
 
totalForceECMX.resize(maxTotalNodes,0.0); 
totalForceECMY.resize(maxTotalNodes,0.0);

totalExplicitForceECMX.resize(maxTotalNodes,0.0); 
totalExplicitForceECMY.resize(maxTotalNodes,0.0);


rHSX.resize(maxTotalNodes,0.0); 
rHSY.resize(maxTotalNodes,0.0);
//memNodeType.resize(maxTotalNodes,notAssigned1) ; 


nodeIsActive.resize(maxTotalNodes,true) ; 
thrust::sequence (indexECM.begin(),indexECM.begin()+maxTotalNodes);
 
thrust::copy(posXIni_ECM.begin(),posXIni_ECM.end(),nodeECMLocX.begin()) ; 
thrust::copy(posYIni_ECM.begin(),posYIni_ECM.end(),nodeECMLocY.begin()) ; 
thrust::copy(eNodeVec.begin(),eNodeVec.end(),peripORexcm.begin()) ;

// for (int i = 0; i < posXIni_ECM.size(); i++){
// 	isActiveECM[i] == true;
// }

AssignDampCoef() ; 
cout << "GPU level initial coordinates and type of external nodes are: " << endl ; 
for (int i=0;  i<nodeECMLocX.size() ;  i++) {
	cout<< nodeECMLocX[i]<<", "<<nodeECMLocY[i]<<", "<<peripORexcm[i] << endl; 
}


PrintECM(0.0) ;
std::string cSVFileName = "./ECMFolder/EnergyExport_" + uniqueSymbol + ".CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str());

			EnergyExport <<"Time,"<<"TotalMorseEnergyECM," << "TotalAdhEnergyECM,"<<"TotalLinSpringEnergy,"<<"TotalEnergy, " <<"TotalEnergyDerivative"<< std::endl;


} //initilaization function finished






void SceECM:: ApplyECMConstrain(int currentActiveCellCount, int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_CoefCell, bool cellPolar, bool subCellPolar, bool isInitPhase, double mitoticThreshold){  

        bool implicit_solver_active = false ; 
	if (eCMRemoved) {
		PrintECMRemoved(curTime);
		cout << "ECM is removed" << endl ; 
		return ; 
	}

#ifdef debugModeECM 
	cudaEvent_t start1, start2, start3, start4, start5, start6, start7, start8, stop;
	float elapsedTime1, elapsedTime2, elapsedTime3, elapsedTime4, elapsedTime5, elapsedTime6,  elapsedTime7 , elapsedTime8 ; 
	cudaEventCreate(&start1);
	cudaEventCreate(&start2);
	cudaEventCreate(&start3);
	cudaEventCreate(&start4);
	cudaEventCreate(&start5);
	cudaEventCreate(&start6);
	cudaEventCreate(&start7);
	cudaEventCreate(&start8);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start1, 0);
#endif
nodeCellLocXOld.resize(totalNodeCountForActiveCellsECM) ;
nodeCellLocYOld.resize(totalNodeCountForActiveCellsECM) ;
adhPairECM_Cell.resize(totalNodeCountForActiveCellsECM,-1) ;
morseEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
adhEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
thrust::copy(nodesPointerECM->getInfoVecs().nodeLocX.begin(),nodesPointerECM->getInfoVecs().nodeLocX.begin()+totalNodeCountForActiveCellsECM,nodeCellLocXOld.begin()) ; 
thrust::copy(nodesPointerECM->getInfoVecs().nodeLocY.begin(),nodesPointerECM->getInfoVecs().nodeLocY.begin()+totalNodeCountForActiveCellsECM,nodeCellLocYOld.begin()) ; 
#ifdef debugModeECM
	cudaEventRecord(start2, 0);
	cudaEventSynchronize(start2);
	cudaEventElapsedTime(&elapsedTime1, start1, start2);
#endif
// std::cout<<"ERROR 1"<<std::endl;

counter ++ ; 
//if (counter>=100 || curTime<(100*dt) || isECMNeighborSet==false) {
if (curTime<(100*dt) || isECMNeighborSet==false) {
	isECMNeighborSet=true ; 
	counter=0 ;
	FindNeighborCandidateForCellsAndECMNodes(); 
	// int ECMAdd_indx = decideIfAddECMNode_M(numNodesECM);
	// if (ECMAdd_indx >= 0){
	// 	AddECMNode(ECMAdd_indx, numNodesECM);
	// 	numNodesECM += 1;
	// 	std::cout<<"New ecm node added!"<<std::endl;
	// 	std::cout<<"Current numECMNodes : "<<numNodesECM<<std::endl;
	// }
	// FindNeighborCandidateForCellsAndECMNodes();
	// for (int i = 0; i < numNodesECM; i++){
	// 	std::cout<<cellNeighborId[i]<<std::endl;
	// }
}
// std::cout<<"ERROR 2"<<std::endl;
if (cellsPointerECM->getCellInfoVecs().isPostDivision ==true || cellsPointerECM->getCellInfoVecs().isPostAddMembrNodes == true){
	std::cout<<"Resetting ecm and cell neighbor info! post division!"<<std::endl;
	int ECMAdd_indx = decideIfAddECMNode_M(numNodesECM);
	// for (int i = 0; i < 5; i++){
		if (ECMAdd_indx >= 0){
			AddECMNode(ECMAdd_indx, numNodesECM);
			numNodesECM += 1;
			std::cout<<"New ecm node added!"<<std::endl;
			std::cout<<"Current numECMNodes : "<<numNodesECM<<std::endl;
		}
	// }
	FindNeighborCandidateForCellsAndECMNodes();
	isECMNeighborResetPostDivision=true;
	// for (int i = 0; i < numNodesECM; i++){
	// 	std::cout<<cellNeighborId[i]<<std::endl;
	// }
	// std::cout<<"SponLen : "<<std::endl;
	// for (int i = 0; i < numNodesECM; i++){
	// 	std::cout<<sponLen[i]<<std::endl;
	// }
}

// Ethan's Edits new function that assigns stiffness and rest length to ECM nodes
AssignECMNodeStiffnessAndRestLength();

// std::cout<<"ERROR 3"<<std::endl;
#ifdef debugModeECM
	cudaEventRecord(start3, 0);
	cudaEventSynchronize(start3);
	cudaEventElapsedTime(&elapsedTime2, start2, start3);
#endif

MoveCellNodesByECMForces(totalNodeCountForActiveCellsECM,currentActiveCellCount,dt, Damp_CoefCell, mitoticThreshold) ; 

// std::cout<<"ERROR 4"<<std::endl;
/* To reduce computational cost
energyECM.totalMorseEnergyCellECM = thrust::reduce( morseEnergyCell.begin(),morseEnergyCell.begin()+totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() ); 
energyECM.totalAdhEnergyCellECM   = thrust::reduce( adhEnergyCell.begin()  ,adhEnergyCell.begin()  +totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() );
*/

CalLinSpringForce(); 
// std::cout<<"ERROR 5"<<std::endl;
CalBendSpringForce();
// std::cout<<"ERROR 6"<<std::endl;
#ifdef debugModeECM
	cudaEventRecord(start4, 0);
	cudaEventSynchronize(start4);
	cudaEventElapsedTime(&elapsedTime3, start3, start4);
#endif
CalCellForcesOnECM(mitoticThreshold) ;
// std::cout<<"ERROR 7"<<std::endl;
//energyECM.totalLinSpringEnergyECM = 0.5 * ( thrust::reduce( linSpringEnergy.begin(),linSpringEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() )); 
//to make sure it is based on the distance used for action force calculation.
/* To reduce computational cost 
energyECM.totalMorseEnergyECMCell = thrust::reduce( morseEnergy.begin(),morseEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() ); 
energyECM.totalAdhEnergyECMCell   = thrust::reduce( adhEnergy.begin()  ,adhEnergy.begin()  +numNodesECM,(double) 0.0, thrust::plus<double>() );
*/

if (!implicit_solver_active) {
// Calculate summation of all forces and move nodes if explicit solver is going to be used
    CalSumForcesOnECM() ;
    MoveNodesBySumAllForces(dt) ; 
}
// std::cout<<"ERROR 8"<<std::endl;
// if (implicit_solver_active) {
// //Calculate right hand side of implicit solver which includes explicit forces
//     CalSumOnlyExplicitForcesOnECM() ;
//     CalRHS(dt) ;
// }

#ifdef debugModeECM
	cudaEventRecord(start5, 0);
	cudaEventSynchronize(start5);
	cudaEventElapsedTime(&elapsedTime4, start4, start5);
#endif

//Create tmp CPU vectors for using in implicit solver. Declariation is not acceptable to be inisde the if condition
    // vector <double> tmpRHSX(numNodesECM); 
    // vector <double> tmpRHSY(numNodesECM); 
    // tmpHostNodeECMLocX.resize(numNodesECM); 
    // tmpHostNodeECMLocY.resize(numNodesECM); 
// if (implicit_solver_active) {
// // Copy ECM locations from GPU to CPU if implicit solver is used    
	
//     thrust::copy (rHSX.begin(), rHSX.begin()+numNodesECM, tmpRHSX.begin()); 
//     thrust::copy (rHSY.begin(), rHSY.begin()+numNodesECM, tmpRHSY.begin()); 
//     thrust::copy (nodeECMLocX.begin(), nodeECMLocX.begin()+numNodesECM, tmpHostNodeECMLocX.begin()); 
//     thrust::copy (nodeECMLocY.begin(), nodeECMLocY.begin()+numNodesECM, tmpHostNodeECMLocY.begin());
//     //cout << "max RHSX is " << *max_element(tmpRHSX.begin(), tmpRHSX.begin()+numNodesECM) << endl ;  
//     //cout << "min RHSX is " << *min_element(tmpRHSX.begin(), tmpRHSX.begin()+numNodesECM) << endl ;
//     //cout << "max RHSY is " << *max_element(tmpRHSY.begin(), tmpRHSY.begin()+numNodesECM) << endl ;  
//     //cout << "min RHSY is " << *min_element(tmpRHSY.begin(), tmpRHSY.begin()+numNodesECM) << endl ; 
// }
    #ifdef debugModeECM
	cudaEventRecord(start6, 0);
	cudaEventSynchronize(start6);
	cudaEventElapsedTime(&elapsedTime5, start5, start6);
    #endif

// if (implicit_solver_active) {
//     // setting up eqaution of motion if implicit solver is used
//     EquMotionCoef (dt); 
// }

    #ifdef debugModeECM
	cudaEventRecord(start7, 0);
	cudaEventSynchronize(start7);
	cudaEventElapsedTime(&elapsedTime6, start6, start7);
    #endif

// if (implicit_solver_active) {
//     // Fetch the implicit solver and update ECM location if implicit solver is used
//     tmpHostNodeECMLocX =solverPointer->SOR3DiagPeriodic(nodeIsActive,hCoefLd, hCoefD, hCoefUd,tmpRHSX,indexPrev, indexNext, tmpHostNodeECMLocX); 
//     tmpHostNodeECMLocY =solverPointer->SOR3DiagPeriodic(nodeIsActive,hCoefLd, hCoefD, hCoefUd,tmpRHSY,indexPrev,indexNext, tmpHostNodeECMLocY);

//     // copy ECM node locations back from CPU to GPU if implicit solver is used    
//     thrust::copy (tmpHostNodeECMLocX.begin(), tmpHostNodeECMLocX.begin()+numNodesECM, nodeECMLocX.begin()); 
//     thrust::copy (tmpHostNodeECMLocY.begin(), tmpHostNodeECMLocY.begin()+numNodesECM, nodeECMLocY.begin());

// }
    #ifdef debugModeECM
	cudaEventRecord(start8, 0);
	cudaEventSynchronize(start8);
	cudaEventElapsedTime(&elapsedTime7, start7, start8);
    #endif
/* To reduce computational cost
cout << "total Morse energy for cell-ECM is= "<< energyECM.totalMorseEnergyCellECM << endl ; 
cout << "total Morse energy for ECM-cell  is= "<< energyECM.totalMorseEnergyECMCell << endl ;
cout << "total adhesion energy for cell-ECM is= "<<  energyECM.totalAdhEnergyCellECM << endl ; 
cout << "total adhesion energy for ECM-cell  is= "<< energyECM.totalAdhEnergyECMCell << endl ; 
//assert (abs (energyECM.totalMorseEnergyCellECM-energyECM.totalMorseEnergyECMCell)<1.0) ;
//assert (abs (energyECM.totalAdhEnergyCellECM-  energyECM.totalAdhEnergyECMCell)  <1.0) ;


if (  (abs (energyECM.totalMorseEnergyCellECM-energyECM.totalMorseEnergyECMCell)>1.0) || 
	  (abs (energyECM.totalAdhEnergyCellECM-  energyECM.totalAdhEnergyECMCell)  >1.0)
   ) {

	cout << "Warning: Action and reaction forces in the ECM do not match each other" << endl ; 
}
*/ 

# ifdef debugModeECM 
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime8, start8, stop);
	std::cout << "time 1 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime1 << endl ; 
	std::cout << "time 2 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime2 << endl ; 
	std::cout << "time 3 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime3 << endl ; 
	std::cout << "time 4 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime4 << endl ; 
	std::cout << "time 5 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime5 << endl ; 
	std::cout << "time 6 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime6 << endl ; 
	std::cout << "time 7 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime7 << endl ; 
	std::cout << "time 8 spent in ECM module for moving the membrane node of cells and ECM nodes are: " << elapsedTime8 << endl ; 
#endif

//throw std::invalid_argument(" Solver called properly and I want to stop the code"); 
PrintECM(curTime); 
// std::cout<<"ERROR 9"<<std::endl;
}

void  SceECM:: PrintECM(double curTime) {
		lastPrintECM=lastPrintECM+1 ; 
        if (lastPrintECM>=freqPlotData) { 
			outputFrameECM++ ; 
			lastPrintECM=0 ;
			cout << " I am in regular print function" << endl ; 
			// First ECM output file for paraview //
			std::string vtkFileName = "./ECMFolder/ECM_" + uniqueSymbol +patch::to_string(outputFrameECM-1) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << nodeECMLocX.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
				<< 0.0 << std::endl;
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()<< " " << 3 *nodeECMLocX.size()<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << indexECM[i] << " "
				<< indexECM[i+1] << std::endl;
			}
			ECMOut << 2 << " " << indexECM[nodeECMLocX.size()-1] << " "<< indexECM[0] << std::endl; //last point to the first point


			ECMOut << "CELL_TYPES " << nodeECMLocX.size()<< endl;
			for (uint i = 0; i < nodeECMLocX.size() ; i++) {
				ECMOut << "3" << endl;
			}
			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
			ECMOut << "SCALARS Avg_Tension " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<linSpringAvgTension[i] <<endl ; 
			}
			
			ECMOut << "SCALARS Node_Type " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<peripORexcm[i] <<endl ; 
			}

			ECMOut.close();



			// second output file for curvature estimation //
			std::string txtFileName = "./ECMFolder/ECMLocationExport_" + uniqueSymbol+ patch::to_string(outputFrameECM-1) + ".txt";
			ofstream ECMLocationExport ;
			ECMLocationExport.open(txtFileName.c_str());
			//ECMExport << "ECM pouch coordinates" << std::endl;

			for (uint i = 0; i < nodeECMLocX.size(); i++) {
		//		if (peripORexcm[i]==excm) {
					ECMLocationExport<< nodeECMLocX[i] << " " << nodeECMLocY[i] << " " << 0.0 << " "<< peripORexcm[i]<<std::endl;
		//		}
			}

			//ECMExport << "ECM lumen side coordinates" << std::endl;
		//	for (uint i = 0; i < nodeECMLocX.size(); i++) {
		//		if (peripORexcm[i]==perip) {
		//			ECMLocationExport << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
		//			<< 0.0 << std::endl;
		//		}
		//	}

			ECMLocationExport.close();
		

			//Third write file for ECM
			txtFileName = "./ECMFolder/ECMTensionExport_" + uniqueSymbol+ patch::to_string(outputFrameECM-1) + ".txt";
			ofstream ECMTensionExport ;
			ECMTensionExport.open(txtFileName.c_str());

			for (uint i = 0; i < nodeECMLocX.size(); i++) {
					ECMTensionExport<< linSpringAvgTension[i]<< " " << peripORexcm[i]<< std::endl;
				}

			ECMTensionExport.close();
			///
			//Fourth write file for ECM
			energyECM.totalEnergyECMOld=energyECM.totalEnergyECM ; 
			energyECM.totalEnergyECM=       energyECM.totalMorseEnergyECMCell
			                         +      energyECM.totalAdhEnergyECMCell
								     +      energyECM.totalLinSpringEnergyECM ;
									 

			std::string cSVFileName = "./ECMFolder/EnergyExport_" + uniqueSymbol+ ".CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str(),ofstream::app);
			
			//EnergyExport <<"totalMorseEnergyCell " << "totalAdhEnergyCell "<<  "totalMorseEnergy "<<"totalAdhEnergy "<< "totalLinSpringEnergy " << std::endl;
			EnergyExport <<curTime<<","<<energyECM.totalMorseEnergyECMCell << "," << energyECM.totalAdhEnergyECMCell<< "," << energyECM.totalLinSpringEnergyECM <<"," << energyECM.totalEnergyECM <<","<<energyECM.totalEnergyPrimeECM <<std::endl;


		}

}

// This is just to create a file to be able to generate the movie with consisten frames
void  SceECM:: PrintECMRemoved(double curTime) {
		lastPrintECM=lastPrintECM+1 ; 
        if (lastPrintECM>=freqPlotData) { 
			outputFrameECM++ ; 
			lastPrintECM=0 ;

			cout << " I am in ECM removed print function" << endl ; 
			// First ECM output file for paraview //
			std::string vtkFileName = "./ECMFolder/ECM_" + uniqueSymbol +patch::to_string(outputFrameECM-1) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << nodeECMLocX.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut << -500.0  << " " << -500.0  << " "    
				<< 0.0 << std::endl; // Just out of domain
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()<< " " << 3 *nodeECMLocX.size()<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << indexECM[i] << " "
				<< indexECM[i+1] << std::endl;
			}
			ECMOut << 2 << " " << indexECM[nodeECMLocX.size()-1] << " "<< indexECM[0] << std::endl; //last point to the first point


			ECMOut << "CELL_TYPES " << nodeECMLocX.size()<< endl;
			for (uint i = 0; i < nodeECMLocX.size() ; i++) {
				ECMOut << "3" << endl;
			}
			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
			ECMOut << "SCALARS Avg_Tension " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<linSpringAvgTension[i] <<endl ; 
			}
			
			ECMOut << "SCALARS Node_Type " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<peripORexcm[i] <<endl ; 
			}

			ECMOut.close();

		}

}




AniResumeData  SceECM:: obtainResumeData() {
	AniResumeData aniResumeData ;
	thrust:: host_vector<double> hostTmpLocX; 
	thrust:: host_vector<double> hostTmpLocY; 
	thrust:: host_vector<EType>  hostTmpType; 
    
	hostTmpLocX.resize(numNodesECM) ; 
	hostTmpLocY.resize(numNodesECM) ; 
	hostTmpType.resize(numNodesECM) ; 
	cout << " I am in obtainResumeData function" << endl ; 
	thrust::copy ( 
		thrust::make_zip_iterator(
			thrust::make_tuple(nodeECMLocX.begin(),nodeECMLocY.begin(),peripORexcm.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(nodeECMLocX.begin(),nodeECMLocY.begin(),peripORexcm.begin()))+numNodesECM,
		thrust::make_zip_iterator(
			thrust::make_tuple(hostTmpLocX.begin(),hostTmpLocY.begin(),hostTmpType.begin()))); 

	cout << " I start passing to regular vector variables" << endl ; 
	CVector tmp; 
	for( int i=0  ; i<numNodesECM ; i++) {
		tmp=CVector (hostTmpLocX[i], hostTmpLocY[i], 0.0) ; 
		aniResumeData.nodePosArr.push_back(tmp) ; 
		aniResumeData.nodeECMType.push_back(hostTmpType[i]) ; 
	}
	return aniResumeData ; 
}

void SceECM::EquMotionCoef (double dt) {

   vector <double> sponLenHost(numNodesECM) ;
   vector <double> sponLenWithNext ; 
   vector <double> sponLenWithPrev ; 
   vector <double> distWithNext ; 
   vector <double> distWithPrev ;
   vector <double> dampCoefHost ; 


   sponLenWithNext.clear(); 
   sponLenWithPrev.clear(); 
   distWithNext.clear() ; 
   distWithPrev.clear() ; 
   hCoefLd.clear() ; 
   hCoefUd.clear() ;  
   hCoefD.clear()  ;
   indexNext.clear() ; 
   indexPrev.clear() ; 
   dampCoefHost.clear() ; 

   indexNext.resize(numNodesECM) ; 
   indexPrev.resize(numNodesECM) ; 
   dampCoefHost.resize(numNodesECM) ; 

   thrust::copy(sponLen.begin(),sponLen.begin()+numNodesECM, sponLenHost.begin()) ; 
   thrust::copy(dampCoef.begin(),dampCoef.begin()+numNodesECM, dampCoefHost.begin()) ; 


   double k=stiffLevel[0] ; //Assumming ECM is homogenous in mechanical properties

   for ( int i=0 ;  i< numNodesECM ; i++) {
	   indexNext.at(i)=i+1 ;
	   indexPrev.at(i)=i-1 ;
	   if (i==numNodesECM-1){
	      indexNext.at(i)=0 ; 
	   }
	   if (i==0){
	      indexPrev.at(i)=numNodesECM-1 ; 
	   }
	   sponLenWithNext.push_back( 0.5*(sponLenHost[indexNext.at(i)]+sponLenHost[i]) ); 
	   sponLenWithPrev.push_back( 0.5*(sponLenHost[indexPrev.at(i)]+sponLenHost[i]) );

	   distWithNext.push_back(sqrt( pow(tmpHostNodeECMLocX[indexNext.at(i)]-tmpHostNodeECMLocX[i],2) + 
	                                pow(tmpHostNodeECMLocY[indexNext.at(i)]-tmpHostNodeECMLocY[i],2))) ;
	   distWithPrev.push_back(sqrt( pow(tmpHostNodeECMLocX[indexPrev.at(i)]-tmpHostNodeECMLocX[i],2) + 
	                                pow(tmpHostNodeECMLocY[indexPrev.at(i)]-tmpHostNodeECMLocY[i],2)));  
   }

   for ( int i=0 ;  i< numNodesECM ; i++) {
      hCoefD.push_back (1 + k*dt/dampCoefHost.at(i)*( 2 - sponLenWithPrev.at(i)/(distWithPrev.at(i) + 0.0001 )   
	                                                    - sponLenWithNext.at(i)/(distWithNext.at(i) + 0.0001 ))) ; 
	  hCoefLd.push_back(    k*dt/dampCoefHost.at(i)*(-1 + sponLenWithPrev.at(i)/(distWithPrev.at(i) + 0.0001 ))) ; 
	  hCoefUd.push_back(    k*dt/dampCoefHost.at(i)*(-1 + sponLenWithNext.at(i)/(distWithNext.at(i) + 0.0001 ))) ; 
   }
  
#ifdef debugModeECM
   cout <<"max distance with next node is" <<*max_element ( distWithNext.begin(), distWithNext.begin()+numNodesECM) <<endl ;  
   cout <<"min distance with next node is" << *min_element ( distWithNext.begin(), distWithNext.begin()+numNodesECM) <<endl  ;

   cout <<"max distance with previous node is" <<*max_element ( distWithPrev.begin(), distWithPrev.begin()+numNodesECM) <<endl ;  
   cout <<"min distance with previous node is" << *min_element ( distWithPrev.begin(), distWithPrev.begin()+numNodesECM) <<endl  ; 

   vector < double> hCoefDAbs;
   hCoefDAbs.clear() ; 
   for ( int i=0 ;  i< numNodesECM ; i++) {
      hCoefDAbs.push_back (abs(1 + k*dt/dampCoefHost.at(i)*( 2 - sponLenWithPrev.at(i)/(distWithPrev.at(i) + 0.0001 )   
	                                                           - sponLenWithNext.at(i)/(distWithNext.at(i) + 0.0001 )))) ; 
   }
   cout <<"max main  diag. elment is "  << *max_element ( hCoefD.begin(),  hCoefD.begin() +numNodesECM) <<endl ;  
   cout <<"min main  diag. element is " << *min_element ( hCoefD.begin(),  hCoefD.begin() +numNodesECM) <<endl  ;
   cout <<"min main  Abs(diag.) element is " << *min_element ( hCoefDAbs.begin(),  hCoefDAbs.begin() +numNodesECM) <<endl  ;
   cout <<"max upper diag. element is " << *max_element ( hCoefUd.begin(), hCoefUd.begin()+numNodesECM) <<endl  ;
   cout <<"min upper diag. element is " << *min_element ( hCoefUd.begin(), hCoefUd.begin()+numNodesECM) <<endl  ;
   cout <<"max lower diag. element is " << *max_element ( hCoefLd.begin(), hCoefLd.begin()+numNodesECM) <<endl  ;
   cout <<"min lower diag. element is " << *min_element ( hCoefLd.begin(), hCoefLd.begin()+numNodesECM) <<endl  ;
   cout <<"stiffness, time step and first element of damping vector is " << endl ; 
   cout << k <<","<< dt<<"," << dampCoefHost.at(0) << endl  ; 
   cout << "constants for stiffness matrix calculated " << endl ; 
   cout << "last diagonal element is " << hCoefD.at(numNodesECM-1) << endl ;
   cout << " number of ECM nodes is "<< numNodesECM << endl ; 
# endif  
}

// Ethan's Edits new function that assigns stiffness and rest length to ECM nodes
void SceECM::AssignECMNodeStiffnessAndRestLength() {
int* cellNeighborIdPtr = thrust::raw_pointer_cast(cellNeighborId.data());
int* cellSubdomainIndxPtr = thrust::raw_pointer_cast(cellsPointerECM->getCellInfoVecs().cellSubdomainIndx.data());

thrust::transform(
    thrust::make_zip_iterator(thrust::make_tuple(peripORexcm.begin(), thrust::counting_iterator<int>(0))),
    thrust::make_zip_iterator(thrust::make_tuple(peripORexcm.begin(), thrust::counting_iterator<int>(0))) + numNodesECM,
    thrust::make_zip_iterator(
        thrust::make_tuple(stiffLevel.begin(), sponLen.begin())
    ),
	// MechProp()
    MechProp(cellNeighborIdPtr, cellSubdomainIndxPtr)
);
	// cout << " Mechanical properties after assignment is " << stiffLevel[0] << endl ; 

}

void SceECM::MoveCellNodesByECMForces(int totalNodeCountForActiveCellsECM,int currentActiveCellCount, double dt, double Damp_CoefCell, double mitoticThreshold) 
{
double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 
// bool* isActiveECM = thrust::raw_pointer_cast(
// 			&isActiveECM[0]);

EType* peripORexcmAddr= thrust::raw_pointer_cast (
			&peripORexcm[0]) ; 
// double* nodeGrowProAddr = thrust::raw_pointer_cast(
// 			&nodesPointerECM->getInfoVecs().nodeGrowPro[0]);
double* cellGrowthProgress = thrust::raw_pointer_cast(
			&cellsPointerECM->getCellInfoVecs().growthProgress[0]);

// move the nodes of epithelial cells 
//// find the closest ECM node to each each cell //

 int numCells = cellsPointerECM->getCellInfoVecs().basalLocX.size() ;

thrust::counting_iterator<int> iBegin(0) ; 
thrust::counting_iterator<int> iBegin2(0) ; 
//////////////////////////////////////////
 thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_permutation_iterator(
						cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor2(
												maxAllNodePerCell))),
					make_transform_iterator (iBegin,
							DivideFunctor2(maxAllNodePerCell)),
					make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					nodesPointerECM->getInfoVecs().nodeLocX.begin(),
					nodesPointerECM->getInfoVecs().nodeLocY.begin(), 
					nodesPointerECM->getInfoVecs().nodeIsActive.begin(),
					nodesPointerECM->getInfoVecs().memNodeType1.begin()
					)), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_permutation_iterator(
						cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
									make_transform_iterator(iBegin2,
											DivideFunctor2(
												maxAllNodePerCell))),
					make_transform_iterator (iBegin,
							DivideFunctor2(maxAllNodePerCell)),
					 make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					 nodesPointerECM->getInfoVecs().nodeLocX.begin(),
                     nodesPointerECM->getInfoVecs().nodeLocY.begin(),
					 nodesPointerECM->getInfoVecs().nodeIsActive.begin(),
					 nodesPointerECM->getInfoVecs().memNodeType1.begin()
					 ))+totalNodeCountForActiveCellsECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodesPointerECM->getInfoVecs().nodeLocX.begin(),
					nodesPointerECM->getInfoVecs().nodeLocY.begin(),
					adhPairECM_Cell.begin(),
					morseEnergyCell.begin(),
					adhEnergyCell.begin())),
				MoveNodes2_Cell(nodeECMLocXAddr,nodeECMLocYAddr,maxMembrNodePerCell,numNodesECM,dt,Damp_CoefCell,peripORexcmAddr,currentActiveCellCount, cellGrowthProgress, mitoticThreshold));//, isActiveECM));
}


void SceECM::CalLinSpringForce()
{

double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 
double* stiffLevelAddr=thrust::raw_pointer_cast (
			&stiffLevel[0]) ; 
double*  sponLenAddr =thrust::raw_pointer_cast (
			&sponLen[0]) ; 
// bool* isActiveECM = thrust::raw_pointer_cast(
// 			&isActiveECM[0]);


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					linSpringAvgTension.begin(),
					linSpringEnergy.begin())),
				LinSpringForceECM(numNodesECM,nodeECMLocXAddr,nodeECMLocYAddr,stiffLevelAddr,sponLenAddr));//, isActiveECM));

//////////////////////////////////// find the closest Cell to each ECM node ///////////



///////////////////////////////////

//cout << " I am after FindCellNeighbor functor" << endl ; 

}

void SceECM::CalBendSpringForce()
{

const double eCMBendStiff=6.0 ; // need to be an input

double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 
// bool* isActiveECM = thrust::raw_pointer_cast(
// 			&isActiveECM[0]);

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					fBendCenterX.begin(),
					fBendCenterY.begin(),
					fBendLeftX.begin(),
					fBendLeftY.begin(),
					fBendRightX.begin(),
					fBendRightY.begin())),
				CalBendECM(nodeECMLocXAddr,nodeECMLocYAddr,numNodesECM,eCMBendStiff));//, isActiveECM));

double* fBendLeftXAddr= thrust::raw_pointer_cast (
			&fBendLeftX[0]) ; 
double* fBendLeftYAddr= thrust::raw_pointer_cast (
			&fBendLeftY[0]) ; 
double* fBendRightXAddr= thrust::raw_pointer_cast (
			&fBendRightX[0]) ; 
double* fBendRightYAddr= thrust::raw_pointer_cast (
			&fBendRightY[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					fBendCenterX.begin(),
					fBendCenterY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 fBendCenterX.begin(),
                                         fBendCenterY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin())),
				SumBendForce(fBendLeftXAddr,fBendLeftYAddr,fBendRightXAddr,fBendRightYAddr,numNodesECM));


}
void SceECM::CalCellForcesOnECM(double mitoticThreshold) 
{
bool* nodeIsActiveAddr= thrust::raw_pointer_cast (
			& (nodesPointerECM->getInfoVecs().nodeIsActive[0])) ; 

int * adhPairECM_CellAddr= thrust::raw_pointer_cast (
			&adhPairECM_Cell[0]) ; 

//Old locations are chosen to make sure action-reaction balance of forces between ECM and cell nodes are fully satisfied. 
double* nodeCellLocXAddr= thrust::raw_pointer_cast (
			&nodeCellLocXOld[0]) ; 
double* nodeCellLocYAddr= thrust::raw_pointer_cast (
			&nodeCellLocYOld[0]) ;
// bool* isActiveECM = thrust::raw_pointer_cast(
// 			&isActiveECM[0]);
double* nodeCellGrowProAddr = thrust::raw_pointer_cast(
			&nodesPointerECM->getInfoVecs().nodeGrowPro[0]);
 

int numCells = cellsPointerECM->getCellInfoVecs().basalLocX.size() ;

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					cellNeighborId.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                     nodeECMLocY.begin(),
					 cellNeighborId.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin(),
					morseEnergy.begin(),
					adhEnergy.begin())),
				MorseAndAdhForceECM(numCells,maxAllNodePerCell,maxMembrNodePerCell,nodeCellLocXAddr,nodeCellLocYAddr,nodeIsActiveAddr,adhPairECM_CellAddr/*, isActiveECM*/, nodeCellGrowProAddr, mitoticThreshold));
}

void SceECM::CalSumForcesOnECM()
{
double dummy=0.0 ;

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					totalForceECMX.begin(),
					totalForceECMY.begin())),
				TotalECMForceCompute(dummy));
}

void SceECM::CalSumOnlyExplicitForcesOnECM() {
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					totalExplicitForceECMX.begin(),
					totalExplicitForceECMY.begin())),
				TotalExplicitECMForceCompute());
}


void SceECM::CalRHS(double dt)
{
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					totalExplicitForceECMX.begin(),
					totalExplicitForceECMY.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					dampCoef.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					totalExplicitForceECMX.begin(),
					totalExplicitForceECMY.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					dampCoef.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					rHSX.begin(),
					rHSY.begin())),
				RHSCompute(dt));
	
}

void  SceECM::MoveNodesBySumAllForces(double dt)
{ 
// move the nodes of ECM 

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					dampCoef.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					dampCoef.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
				MoveNodesECM(dt));
}
void SceECM::FindNeighborCandidateForCellsAndECMNodes() 
{
double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 
double * basalCellLocXAddr= thrust::raw_pointer_cast (
                            & ( cellsPointerECM->getCellInfoVecs().basalLocX[0]) ) ; 
double * basalCellLocYAddr= thrust::raw_pointer_cast ( 
                            & ( cellsPointerECM->getCellInfoVecs().basalLocY[0]) ) ;
EType* peripORexcmAddr= thrust::raw_pointer_cast (
			&peripORexcm[0]) ; 
// bool * isActiveECMAddr = thrust::raw_pointer_cast(
// 			&isActiveECM[0]);

int numCells = cellsPointerECM->getCellInfoVecs().basalLocX.size() ;


//// find the closest ECM node to each each cell //
thrust:: transform (
		thrust::make_zip_iterator (
					thrust:: make_tuple (
						cellsPointerECM->getCellInfoVecs().basalLocX.begin(),
						cellsPointerECM->getCellInfoVecs().basalLocY.begin())),
		thrust::make_zip_iterator (
					thrust:: make_tuple (
					 	cellsPointerECM->getCellInfoVecs().basalLocX.begin(),
                     	cellsPointerECM->getCellInfoVecs().basalLocY.begin()))+numCells, 
	    cellsPointerECM->getCellInfoVecs().eCMNeighborId.begin(),
		FindECMNeighborPerCell(nodeECMLocXAddr,nodeECMLocYAddr,maxTotalNodes));//, isActiveECMAddr));

	thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
					// isActiveECM.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 nodeECMLocX.begin(),
                     nodeECMLocY.begin()))+numNodesECM,
					//  isActiveECM.begin()))+numNodesECM,
	    cellNeighborId.begin(),
		FindCellNeighborPerECMNode(basalCellLocXAddr,basalCellLocYAddr, numCells));
}

void SceECM::AssignDampCoef() {

   thrust::transform ( peripORexcm.begin() ,peripORexcm.begin() +numNodesECM, dampCoef.begin(), AssignDamping(dampBasal,dampBC,dampApical) ); 
   
#ifdef debugModeECM 
   for (int i=0 ;  i<numNodesECM ; i++) {
      if (dampCoef[i] < smallNumber) {
		 cout << "damping of element " << i << " is " << dampCoef[i] << " which is wrong" <<endl ; 
	     throw::invalid_argument ( "damping coefficients in ECM is not set correctly") ;
	  }
   }
#endif
}

int SceECM::decideIfAddECMNode_M(uint numECMNodes) {
// decide if add ecm node given current active node count
	// uint maxECMNode = numNodesECM*3;
	// bool isInitPhase= nodes->isInitPhase ; 
        //  thrust::transform(
		// 	thrust::make_zip_iterator(
		// 			thrust::make_tuple(maxDistToRiVec.begin(),
		// 							   maxTenIndxTypeVec.begin()
		// 							   )),
		// 	thrust::make_zip_iterator(
		// 			thrust::make_tuple(maxDistToRiVec.begin(),
		// 							   maxTenIndxTypeVec.begin()
		// 							   ))
		// 			+ curActCellCt,
		// 			isECMAddingNode,
		// 	ECMGrowFunc(numECMNodes, maxTotalNodes, maxLengthToAddECMNodes));
	double max_Dist = -9999.9;
	int indx = -1;
	for (int i = 0; i < numECMNodes; i++){
		if (i != numECMNodes-1){
			double tmp_Dist = sqrt((nodeECMLocX[i+1] - nodeECMLocX[i])*(nodeECMLocX[i+1] - nodeECMLocX[i]) + (nodeECMLocY[i+1] - nodeECMLocY[i])*(nodeECMLocY[i+1] - nodeECMLocY[i]));
			if (tmp_Dist > max_Dist && tmp_Dist > 2.0*lknotECMPeripGPU){
				max_Dist = tmp_Dist;
				indx = i;
			}
		}
		else{
			double tmp_Dist = sqrt((nodeECMLocX[0] - nodeECMLocX[i])*(nodeECMLocX[0] - nodeECMLocX[i]) + (nodeECMLocY[0] - nodeECMLocY[i])*(nodeECMLocY[0] - nodeECMLocY[i]));
			if (tmp_Dist > max_Dist && tmp_Dist > 2.0*lknotECMPeripGPU){
				max_Dist = tmp_Dist;
				indx = i;
			}
		}
	}

	return indx;
}

void SceECM::AddECMNode(int indx, uint numECMNodes){
		std::cout<<"insertIndx = "<<indx+1<<std::endl;
		std::cout<<"numECMNodes = "<<numECMNodes<<std::endl;
		uint insertIndx = indx + 1; // newIndex 
		double insertX, insertY;
		insertX = (nodeECMLocX[indx] + nodeECMLocX[indx+1])/2.0;
		insertY = (nodeECMLocY[indx] + nodeECMLocY[indx+1])/2.0;
		uint globalIndxEnd = numECMNodes; //membrane nodes are first. End position based on newindex
		uint globalIndexInsert = insertIndx;
		// if (insertIndx<=iDApical) {  //since the current acrive membrane nodes is one more, it can not be the last ID.
		// 	iDApical=iDApical+1 ;
		// }

		// if (insertIndx<=iDBasal) {
		// 	iDBasal=iDBasal+1 ;
		// }
		for (uint i = globalIndxEnd; i >= globalIndexInsert; i--) {
			// isActiveECM[i] = isActiveECM[i - 1];
			nodeECMLocX[i] = nodeECMLocX[i - 1];
			nodeECMLocY[i] = nodeECMLocY[i - 1];
			peripORexcm[i] = peripORexcm[i-1] ;
			sponLen[i] = sponLen[i-1];
			dampCoef[i] = dampCoef[i-1];
		}
		// isActiveECM[globalIndexInsert] = true;
		nodeECMLocX[globalIndexInsert] = insertX;
		std::cout<<"nodeECMLocX[globalIndexInsert] : "<<nodeECMLocX[globalIndexInsert]<<std::endl;
		nodeECMLocY[globalIndexInsert] = insertY;
		std::cout<<"nodeECMLocY[globalIndexInsert] : "<<nodeECMLocY[globalIndexInsert]<<std::endl;
		peripORexcm[globalIndexInsert] = peripORexcm[globalIndexInsert-1]; // to have the same type of the membrane node as at least one of its neighbors
		std::cout<<"peripORexcm[globalIndexInsert] : "<<peripORexcm[globalIndexInsert]<<std::endl;
		sponLen[globalIndexInsert] = sponLen[globalIndexInsert-1];
		std::cout<<"sponLen[globalIndexInsert] : "<<sponLen[globalIndexInsert]<<std::endl;
		dampCoef[globalIndexInsert] = dampCoef[globalIndexInsert-1];
		std::cout<<"dampCoef[globalIndexInsert] : "<<dampCoef[globalIndexInsert]<<std::endl;

		// if (_memNodeType[globalIndexInsert-1] != apical1 || _memNodeType[globalIndexInsert-1] != basal1){
		// 	_memNodeType[globalIndexInsert] = _memNodeType[globalIndexInsert+1];
		// }
		//return (curActCount + 1);
		// numECMNodes += 1;
	
};

// void SceECM::calECMGrowSpeed_M() {
// // reduce_by_key, find value of max tension and their index
// 	thrust::counting_iterator<uint> iBegin(0);
	
// 	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	
// 	thrust::transform(
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							ECMDistToRi.begin(),
// 							ECMNodeType1.begin())),
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(
// 							ECMbrDistToRi.begin(),
// 							ECMNodeType1.begin()))
// 							+ numECMNodes,
// 			thrust::make_zip_iterator(
// 					thrust::make_tuple(maxDistToRiVec.begin(),
// 							maxTenIndxTypeVec.begin())),
// 			thrust::equal_to<uint>(), MaxWInfo());

	// for (int i=0 ; i<cellInfoVecs.maxDistToRiVec.size() ; i++) {
	// 	cout << "the max distance in cell" << i << " is "<<cellInfoVecs.maxDistToRiVec[i] << endl ;
	// 	cout << "At index "<<cellInfoVecs.maxTenIndxVec[i]<<std::endl; 
	// }
// }