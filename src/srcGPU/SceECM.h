#ifndef SCEECM_H
#define SCEECM_H

#include "commonData.h"
#include "SceNodes.h" 
#include "ConfigParser.h"
#include "Solver.h"
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream> 
#include <algorithm>


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<int,double,double,int> IDDI ; 
typedef thrust ::tuple<int,int,int,double,double,bool,MembraneType1> IIIDDBT ; 
typedef thrust ::tuple<double,double> DD ; 
typedef thrust ::tuple<double,double,double,double> DDDD ; 
typedef thrust ::tuple<double,double,bool> DDB ; 
typedef thrust ::tuple<double,double,int,double,double> DDIDD ; 
typedef thrust ::tuple<double,double,double,double,double,double> DDDDDD ;
typedef thrust ::tuple<double,double,double,double,double> DDDDD ; 
typedef thrust::tuple<double, EType> DT ;


struct MechPara_ECM {
	double sceInterCellCPU_ECM[5] ;
	double wLCParaCPU_ECM[4] ;
	double linSpringCPU_ECM ; 
	double linSpringRestLenCPU_ECM ; 
}; 

class SceCells ; // forward declaration
class SceECM {
	bool isECMAddingNode;
//	SceNodes* nodes;
    double dampBasal,dampBC,dampApical ; 
	vector<bool> nodeIsActive ; 
    vector<double> hCoefLd ; 
    vector<double> hCoefUd ;  
    vector<double> hCoefD  ;
	vector<int> indexPrev ; 
	vector<int> indexNext ;
	vector <double> tmpHostNodeECMLocX; 
	vector <double> tmpHostNodeECMLocY;
	thrust::device_vector<double> dampCoef ; 

	int decideIfAddECMNode_M(uint numECMNodes);
	void AddECMNode(int indx, uint numECMNodes);
	bool   eCMRemoved ; 
	bool   isECMNeighborSet ;
	bool isECMNeighborResetPostDivision;
	void   FindNeighborCandidateForCellsAndECMNodes(); 
	void MoveCellNodesByECMForces(int totalNodeCountForActiveCellsECM,int currentActiveCellCount, double dt, double Damp_CoefCell, double mitoticThreshold);  
	void CalLinSpringForce() ; 
	void CalBendSpringForce() ; 
	void CalCellForcesOnECM(double mitoticThreshold) ; 
	void CalSumForcesOnECM() ; 
	void MoveNodesBySumAllForces(double dt); 
	void CalSumOnlyExplicitForcesOnECM(); 
	void EquMotionCoef( double dt) ;
	void CalRHS(double dt);
	void AssignDampCoef() ; 
public:
	SceECM() ;

	SceNodes * nodesPointerECM ; 
	SceCells * cellsPointerECM ;
	Solver   * solverPointer ; 
	//void Initialize_SceECM(SceNodes * nodes, SceCells * cells) ; 
	void Initialize_SceECM(SceNodes * nodes, SceCells * cells, Solver *solver) {
		nodesPointerECM =nodes ; 
		cellsPointerECM= cells ; 
    	solverPointer=solver ; 
	}


	void SetIfECMIsRemoved(bool eCMRemoved) {
		this->eCMRemoved=eCMRemoved ; 
	}

	bool GetIfECMIsRemoved() 
		const { return eCMRemoved  ;}  
	AniResumeData obtainResumeData(); 

    void ApplyECMConstrain(int currentActiveCellCount, int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_CoefCell, bool cellPolar, bool subCellPolar, bool isInitPhase, double mitoticThreshold) ; 
	void Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM, int freqPlotData, string uniqueSymbol); 
		EType ConvertStringToEType (string eNodeRead) ;
	void PrintECM(double curTime);
	void PrintECMRemoved(double curTime);
	//Solver   * solverPointer ; 
double restLenECMSpring ;
double eCMLinSpringStiff ; 
double restLenECMAdhSpring ; 
double maxLenECMAdhSpring ; 
double kAdhECM ;

double stiffnessECMBasal ;
double stiffnessECMBC ;
double stiffnessECMPerip ;
double lknotECMBasal ;
double lknotECMBC ;
double lknotECMPerip ;

EnergyECMInfo energyECM ; 
int counter ; 
int outputFrameECM ;  
int lastPrintECM ;
int numNodesECM ;
int freqPlotData ; 
uint maxAllNodePerCell ; 
uint maxMembrNodePerCell ;
uint maxTotalNodes ; 

string uniqueSymbol ; 
MechPara_ECM mechPara_ECM ; 
 
thrust::device_vector<int> indexECM ;
thrust::device_vector<int> cellNeighborId ;
thrust::device_vector<double> nodeECMLocX ; 
thrust::device_vector<double> nodeECMLocY ; 
thrust::device_vector<double> nodeECMVelX ; 
thrust::device_vector<double> nodeECMVelY ;
// thrust::device_vector<bool> isActiveECM ;

thrust::device_vector<double> nodeECMTmpLocX ; 
thrust::device_vector<double> nodeECMTmpLocY ; 

//thrust::device_vector<double> nodeDeviceLocX ; 
//thrust::device_vector<double> nodeDeviceLocY ; 
thrust::device_vector<double> nodeCellLocXOld ; 
thrust::device_vector<double> nodeCellLocYOld ;
//thrust::device_vector<MembraneType1> memNodeType ;
thrust::device_vector<int>   adhPairECM_Cell ;
 
thrust::device_vector<bool> nodeIsActive_Cell ; 

thrust::device_vector<double> linSpringForceECMX ; 
thrust::device_vector<double> linSpringForceECMY ; 
thrust::device_vector<double> linSpringAvgTension  ; 
thrust::device_vector<double> linSpringEnergy  ; 
thrust::device_vector<double> morseEnergy  ; 
thrust::device_vector<double> adhEnergy  ;

thrust::device_vector<double> morseEnergyCell ; //it should be cell size 
thrust::device_vector<double> adhEnergyCell  ; // it should be cell size
 
 
thrust::device_vector<double> bendSpringForceECMX ; 
thrust::device_vector<double> bendSpringForceECMY ;

 
thrust::device_vector<double> memMorseForceECMX ; 
thrust::device_vector<double> memMorseForceECMY ;
 
thrust::device_vector<double> fBendCenterX ;
thrust::device_vector<double> fBendCenterY ;
thrust::device_vector<double> fBendLeftX ;
thrust::device_vector<double> fBendLeftY ;
thrust::device_vector<double> fBendRightX ;
thrust::device_vector<double> fBendRightY ;

thrust::device_vector<double> totalForceECMX ; 
thrust::device_vector<double> totalForceECMY ;
thrust::device_vector<double> totalExplicitForceECMX ; 
thrust::device_vector<double> totalExplicitForceECMY ;
thrust::device_vector<EType>  peripORexcm ;

thrust::device_vector<double> rHSX ; 
thrust::device_vector<double> rHSY ;
thrust::device_vector<double> stiffLevel ;
thrust::device_vector<double> sponLen ;
};

__device__
double calMorse_ECM (const double & linkLength); 

__device__
double calMorse_ECM_mitotic (const double & linkLength, double scaling); 

__device__
double calMorseEnergy_ECM (const double & linkLength); 
//__device__
//double calWLC_ECM (const double & linkLength); 


__device__
bool IsValidAdhPair (const double & dist); 

__device__
bool IsValidAdhPairForNotInitPhase (const double & dist); 

__device__
double CalAdhECM (const double & dist);

__device__
double CalAdhEnergy (const double & dist);

__device__
void DefineECMStiffnessAndLknot (EType nodeType, double & stiffness, double & sponLen) ; 

struct InitECMLoc
{
     double  _MinLocX;
     double  _MinDist;

    InitECMLoc (double MinLocX, double MinDist) : _MinLocX(MinLocX), _MinDist(MinDist) {}

   __host__  __device__

        double operator()(const int & x, const double & y)  {
        return (_MinLocX+x*_MinDist) ; 

  }
};


struct AddECMForce
{
    const double  eCMY;


    AddECMForce(double _eCMY) : eCMY(_eCMY) {}

   __host__  __device__

        double operator()( const double & x, const double & y) const {
        if (y<eCMY) {
        return (eCMY) ; 
        }
        else {
        return (y) ; 
         }
  }
};

struct MyFunction: public thrust::unary_function <double,double>{
       double _k, _bound ; 
      

       __host__ __device__ MyFunction(double k , double bound ) :
                             _k(k),_bound(bound) {
	}
       __host__ __device__
       double operator()(const double &y) {
	if (y<=_k && y>_bound) {
		return _k ; 
	}
        else {
                return y; 
	}
	}
} ; 

struct ModuloFunctor2: public thrust::unary_function <int,int>{
       int _dividend ;  
      

       __host__ __device__ ModuloFunctor2(int dividend) :
                             _dividend(dividend) {
	}
       __host__ __device__
       int operator()(const int &num) {
                return num %_dividend;
	} 
} ; 

struct DivideFunctor2: public thrust::unary_function <int,int>{
       int _dividend ;  
      

       __host__ __device__ DivideFunctor2(int dividend) :
                             _dividend(dividend) {
	}
       __host__ __device__
       int operator()(const int &num) {
                return num /_dividend;
	} 
} ; 

struct FindCellNeighborPerECMNode: public thrust::unary_function<DD,int> {

	double * _basalCellLocX;
	double * _basalCellLocY;
	int _numCells ; 

	__host__ __device__ FindCellNeighborPerECMNode(double * basalCellLocX, double * basalCellLocY, int numCells): _basalCellLocX(basalCellLocX), _basalCellLocY(basalCellLocY), _numCells (numCells)  {
	}

	__device__ int  operator() (const DD  & dD) const {
		double  eCMLocX=	thrust::get<0>(dD) ; 
		double  eCMLocY=    thrust::get<1>(dD) ; 
		// bool	isActiveECM= thrust::get<2>(dDB) ; 
		double  distMin= 1000000 ;  // large number
		int idMin=-1 ;
		double dist ; 
		// if (isActiveECM == true){
			for ( int i=0 ;  i<_numCells ; i++) { 
				dist=sqrt((eCMLocX-_basalCellLocX[i])*(eCMLocX-_basalCellLocX[i])+
						(eCMLocY-_basalCellLocY[i])*(eCMLocY-_basalCellLocY[i])) ;
				if (dist <distMin) {
					idMin=i ;
					distMin=dist ; 
				}
			}
		// }

		return idMin ; 
	}
}; 

struct FindECMNeighborPerCell: public thrust::unary_function<DD,int> {

	double * _eCMLocX;
	double * _eCMLocY;
	int _numECMNodes ; 
	// bool * _isActiveECM ;

	__host__ __device__ FindECMNeighborPerCell(double * eCMLocX, double * eCMLocY, int numECMNodes/*, bool * isActiveECM*/): _eCMLocX(eCMLocX), _eCMLocY(eCMLocY), _numECMNodes (numECMNodes)/*, _isActiveECM(isActiveECM)*/  {
	}

	__device__ int  operator() (const DD  & dD) const {
		double  basalCellLocX=    thrust::get<0>(dD) ; 
		double  basalCellLocY=    thrust::get<1>(dD) ; 
		double  distMin= 1000000 ;  // large number
		int idMin=-1 ;
		double dist ; 
		for ( int i=0 ;  i<_numECMNodes ; i++) { 
			// if (_isActiveECM[i] == true){
				dist=sqrt((basalCellLocX-_eCMLocX[i])*(basalCellLocX-_eCMLocX[i])+
						(basalCellLocY-_eCMLocY[i])*(basalCellLocY-_eCMLocY[i])) ;
				if (dist <distMin) {
					idMin=i ;
					distMin=dist ; 
				}
			// }
		}

		return idMin ; 
	}
}; 


struct MoveNodes2_Cell: public thrust::unary_function<IIIDDBT,DDIDD> {
	 double  *_locXAddr_ECM; 
         double  *_locYAddr_ECM; 
        uint _maxMembrNodePerCell ; 
	 int _numNodes_ECM ;
	 double _dt ; 
	 double _Damp_Coef ;
	 EType*  _peripORexcmAddr ;
	 int _activeCellCount ; 
	 double* _cellGrowthProgress;
	double _mitoticThreshold;
	//  bool* _isActiveECM;
	__host__ __device__ MoveNodes2_Cell (double * locXAddr_ECM, double * locYAddr_ECM, uint maxMembrNodePerCell, int numNodes_ECM, double dt, double Damp_Coef,EType * peripORexcmAddr, int activeCellCount, double* cellGrowthProgress, double mitoticThreshold/*, bool* isActiveECM*/) :
				_locXAddr_ECM(locXAddr_ECM),_locYAddr_ECM(locYAddr_ECM),_maxMembrNodePerCell(maxMembrNodePerCell),_numNodes_ECM(numNodes_ECM),_dt(dt),
			    _Damp_Coef(Damp_Coef), _peripORexcmAddr(peripORexcmAddr), _activeCellCount (activeCellCount), _cellGrowthProgress(cellGrowthProgress), _mitoticThreshold(mitoticThreshold)/*, _isActiveECM(isActiveECM)*/	{
	}
	__device__ DDIDD  operator()(const IIIDDBT & iIIDDBT) const {
	int eCMNeighborId=				thrust::get<0>(iIIDDBT) ; 
	int cellRank=					thrust::get<1>(iIIDDBT) ; 
	int nodeRankInOneCell=          thrust::get<2>(iIIDDBT) ; 
	double            locX=         thrust::get<3>(iIIDDBT) ; 
	double            locY=         thrust::get<4>(iIIDDBT) ; 
	bool              nodeIsActive= thrust::get<5>(iIIDDBT) ; 
	MembraneType1     mNodeType=     thrust::get<6>(iIIDDBT) ; 
	
	double locX_ECM, locY_ECM ; 
	double dist ;
	double fMorse ;  
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	double eMorseCell=0 ; 
	double eAdhCell= 0 ; 
	double distMin=10000 ; // large number
	double distMinX, distMinY ; 
	//double kStifMemECM=3.0 ; // need to take out 
	//double distMaxAdh=0.78125; // need to take out
	//double distSponAdh=0.0625 ;  // need to take out
	double fAdhMemECM=0.0   ; 
	double fAdhMemECMX=0.0 ; 
	double fAdhMemECMY=0.0  ; 
	int    adhPairECM=-1 ; //no adhere Pair
	int   iPair=-1 ;
	double smallNumber=0.000001;
	int eCMId ; 
	bool enteringMitotic;
	double scaling;
		
		if ( nodeIsActive && nodeRankInOneCell<_maxMembrNodePerCell ) {
			//for (int i=0 ; i<_numNodes_ECM ; i++) {
			for (int i=eCMNeighborId-150 ; i<eCMNeighborId+150 ; i++) {
				eCMId=i ; 
				// if (_isActiveECM[eCMId]==false){
				// 	continue;
				// }
				if (eCMId>_numNodes_ECM){
					eCMId=eCMId-_numNodes_ECM ;
				}
				if (eCMId<0){
					eCMId=eCMId+_numNodes_ECM ;
				}
				locX_ECM=_locXAddr_ECM[eCMId]; 
				locY_ECM=_locYAddr_ECM[eCMId];
				dist=sqrt((locX-locX_ECM)*(locX-locX_ECM)+(locY-locY_ECM)*(locY-locY_ECM)) ;
				if (_cellGrowthProgress[cellRank] > _mitoticThreshold){
					enteringMitotic = true;
					//scaling = (_nodeCellGrowProAddr[i] - _mitoticThreshold)/(1 - _mitoticThreshold);
					scaling = 1; //change made on 2/19/2022
				}
				else{
					enteringMitotic = false;
					scaling = 1;
				}
				if (enteringMitotic == true){
					fMorse=(1.0+(2.0-1.0)*scaling)*calMorse_ECM(dist);
				}
				else{
					fMorse=calMorse_ECM(dist);
				}
				eMorseCell=eMorseCell + calMorseEnergy_ECM(dist);  
				fTotalMorseX=fTotalMorseX+fMorse*(locX_ECM-locX)/dist ; 
				fTotalMorseY=fTotalMorseY+fMorse*(locY_ECM-locY)/dist ; 
				fTotalMorse=fTotalMorse+fMorse ; 
				if ( dist < distMin) {
					if( mNodeType==basal1  ) {  //adhesion only for basal nodes
						distMin=dist ; 
						distMinX=(locX_ECM-locX) ;
						distMinY=(locY_ECM-locY) ; 
						iPair=eCMId ; 
					}
				}
			}
			if (IsValidAdhPairForNotInitPhase(distMin)&& iPair!=-1) {

        		fAdhMemECM=CalAdhECM(distMin) ; 
				eAdhCell=CalAdhEnergy(distMin) ; 

				fAdhMemECMX=fAdhMemECM*distMinX/distMin ;  
				fAdhMemECMY=fAdhMemECM*distMinY/distMin ; 
				adhPairECM=iPair ; 
			}
	 }

	 return thrust::make_tuple ((locX+(fTotalMorseX+fAdhMemECMX)*_dt/_Damp_Coef),
	 							(locY+(fTotalMorseY+fAdhMemECMY)*_dt/_Damp_Coef),
								 adhPairECM,eMorseCell,eAdhCell )  ; 
}
	

	
} ;

        



struct LinSpringForceECM: public thrust::unary_function<IDD,DDDD> {
         int   _numECMNodes ; 	
         double  *_locXAddr; 
         double  *_locYAddr;
		 double *_stiffAddr; 
		 double * _sponLenAddr ; 
		//  bool * _isActiveECM;


	__host__ __device__ LinSpringForceECM (double numECMNodes, double * locXAddr, double * locYAddr, 
										   double * stiffAddr, double * sponLenAddr/*, bool * isActiveECM*/ ) :
	 _numECMNodes(numECMNodes),_locXAddr(locXAddr),_locYAddr(locYAddr),_stiffAddr (stiffAddr), _sponLenAddr (sponLenAddr)/*, _isActiveECM(isActiveECM)*/ {
	}
	 __device__ DDDD operator()(const IDD & iDD) const {
	
	int     index=    thrust::get<0>(iDD) ; 
	double  locX=     thrust::get<1>(iDD) ; 
	double  locY=     thrust::get<2>(iDD) ; 

	double locXLeft  ; 
	double locYLeft  ;

	double distLeft ;
	double forceLeft  ; 
	double forceLeftX ; 
	double forceLeftY ; 

	double locXRight ; 
	double locYRight ;
 
	double distRight ; 
	double forceRight ; 
	double forceRightX ; 
	double forceRightY ;

	bool activePairLeft ;
	bool activePairRight ;

	double energyLeft, energyRight ; 
	int indexLeft, indexRight ; 

	if (index != 0 && index != _numECMNodes) {
		locXLeft=_locXAddr[index-1] ; 
		locYLeft=_locYAddr[index-1] ;
		indexLeft=index-1 ; 
		// if (_isActiveECM[index]==true && _isActiveECM[indexLeft]==true){
		// 	activePairLeft = true;
		// }
	}
	else if (index == 0){
		locXLeft=_locXAddr[_numECMNodes-1] ;
		locYLeft=_locYAddr[_numECMNodes-1] ;
		indexLeft=_numECMNodes-1 ;
		// if (_isActiveECM[index]==true && _isActiveECM[indexLeft]==true){
		// 	activePairLeft = true;
		// }
	}
	else if (index == _numECMNodes){
		locXLeft=_locXAddr[index-1] ;
		locYLeft=_locYAddr[index-1] ;
		indexLeft=index-1 ;
		// if (_isActiveECM[index]==true && _isActiveECM[indexLeft]==true){
		// 	activePairLeft = true;
		// }
	}
	// else {
	// 	locXLeft=_locXAddr[_numNodes-1] ;
	// 	locYLeft=_locYAddr[_numNodes-1] ;
	// 	indexLeft=_numNodes-1 ; 
	// }

	// if (activePairLeft == true){
		distLeft=sqrt( ( locX-locXLeft )*(locX-locXLeft) +( locY-locYLeft )*(locY-locYLeft) ) ;
		double kStiffLeft =0.5*(  _stiffAddr[indexLeft]  +  _stiffAddr[index]  ) ; 
		double sponLenLeft=0.5*(_sponLenAddr[indexLeft]  +_sponLenAddr[index]  ) ; 
		//	forceLeft=calWLC_ECM(distLeft) ; 
			//forceLeft=_linSpringStiff*(distLeft-_restLen) ; 
			forceLeft=kStiffLeft*(distLeft-sponLenLeft) ; 
			energyLeft=0.5*kStiffLeft*(distLeft-sponLenLeft)*(distLeft-sponLenLeft) ; 
			forceLeftX=forceLeft*(locXLeft-locX)/distLeft ; 
			forceLeftY=forceLeft*(locYLeft-locY)/distLeft ; 
	// }
	
	if (index != _numECMNodes-1) {
		locXRight=_locXAddr[index+1] ; 
		locYRight=_locYAddr[index+1] ;
		indexRight=index+1 ; 
		// if (_isActiveECM[index]==true && _isActiveECM[indexRight]==true){
		// 	activePairRight = true;
		// }
	}
	if (index == _numECMNodes-1){
		locXRight=_locXAddr[0] ; 
		locYRight=_locYAddr[0] ;
		indexRight=0  ;
		// if (_isActiveECM[index]==true && _isActiveECM[indexRight]==true){
		// 	activePairRight = true;
		// }
	}
	else {
		locXRight=_locXAddr[index+1] ; 
		locYRight=_locYAddr[index+1] ;
		indexRight=index+1 ;
		// if (_isActiveECM[index]==true && _isActiveECM[indexRight]==true){
		// 	activePairRight = true;
		// }
	}
	// if (activePairRight==true){
		distRight=sqrt( ( locX-locXRight )*(locX-locXRight) +( locY-locYRight )*(locY-locYRight) ) ; 
		double kStiffRight =0.5*(_stiffAddr  [indexRight]  +  _stiffAddr[index]  ) ; 
		double sponLenRight=0.5*(_sponLenAddr[indexRight]  +_sponLenAddr[index]  ) ; 
				//forceRight=_linSpringStiff*(distRight-_restLen) ; 
		forceRight=kStiffRight*(distRight-sponLenRight) ;
		energyRight=0.5*kStiffRight*(distRight-sponLenRight)*(distRight-sponLenRight)  ;

	//  	forceRight=calWLC_ECM(distRight) ; 
		forceRightX=forceRight*(locXRight-locX)/distRight ; 
		forceRightY=forceRight*(locYRight-locY)/distRight ; 
	//  }
	return thrust::make_tuple(forceLeftX+forceRightX,forceLeftY+forceRightY,0.5*(forceLeft+forceRight), energyLeft+energyRight) ;
//	}
        

}

	
} ;

 
 struct MorseAndAdhForceECM: public thrust::unary_function<IDDI,DDDD> {
         int  _numCells ; 	
         uint  _maxNodePerCell ; 	
         uint  _maxMembrNodePerCell ; 	
         double  *_locXAddr_Cell; 
         double  *_locYAddr_Cell; 
	 bool    *_nodeIsActive_Cell ;  
	 int     *_adhPairECM_Cell ; 
	//  bool * _isActiveECM;
	 double * _nodeCellGrowProAddr;
	 double _mitoticThreshold;
	__host__ __device__ MorseAndAdhForceECM (int numCells, uint maxNodePerCell, uint maxMembrNodePerCell, double * locXAddr_Cell, double * locYAddr_Cell, bool * nodeIsActive_Cell, int * adhPairECM_Cell/*, bool* isActiveECM*/, double* nodeCellGrowProAddr, double mitoticThreshold):
	_numCells(numCells), _maxNodePerCell(maxNodePerCell), _maxMembrNodePerCell(maxMembrNodePerCell),_locXAddr_Cell(locXAddr_Cell),_locYAddr_Cell(locYAddr_Cell),_nodeIsActive_Cell(nodeIsActive_Cell),_adhPairECM_Cell(adhPairECM_Cell), /*_isActiveECM(isActiveECM),*/ _nodeCellGrowProAddr(nodeCellGrowProAddr), _mitoticThreshold(mitoticThreshold) {
	}
	 __device__ 
	DDDD operator()(const IDDI & iDDI) const {
	
	int     index=    thrust::get<0>(iDDI) ; 
	double  locX=     thrust::get<1>(iDDI) ; 
	double  locY=     thrust::get<2>(iDDI) ; 
	int   cellId=     thrust::get<3>(iDDI) ; 

	double fMorse ; 
	double locX_C, locY_C ; 
	double dist ;
	double fAdhX=0 ; 
	double fAdhY=0 ; 
	double eAdh=0 ; 
	double eMorse=0 ; 
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	//double kAdhECM=3 ;  //need to take out 
	//double distAdhSpon=0.0625 ; // need to take out 
	//double distAdhMax=0.78125 ; // need to take out
	double fAdh ; 
	// we are already in active cells. Two more conditions: 1-it is membrane 2-it is active node
	bool enteringMitotic;
	double scaling;

	int cellIdAfter=cellId +1 ; 
	if (cellIdAfter>_numCells-1) {
		cellIdAfter=0 ; 
	}
	int cellIdBefore=cellId -1 ;  
	if (cellIdBefore<0) {
		cellIdBefore=_numCells-1 ; 
	}

    for (int i=cellIdBefore*_maxNodePerCell ; i<(cellIdBefore+1)*_maxNodePerCell ; i++) {
    //for (int i=0 ; i<_numCells*_maxNodePerCell ; i++) {
		if (_nodeCellGrowProAddr[i] > _mitoticThreshold){
			enteringMitotic = true;
			//scaling = (_nodeCellGrowProAddr[i] - _mitoticThreshold)/(1 - _mitoticThreshold);
			scaling = 1; //change made on 2/2/2022
		}
		else{
			enteringMitotic = false;
			scaling = 1;
		}
		if (_nodeIsActive_Cell[i] && (i%_maxNodePerCell)<_maxMembrNodePerCell){
			locX_C=_locXAddr_Cell[i]; 
			locY_C=_locYAddr_Cell[i];
			dist=sqrt((locX-locX_C)*(locX-locX_C)+(locY-locY_C)*(locY-locY_C)) ;
			if (enteringMitotic == true){
				fMorse=(1.0+(2.0-1.0)*scaling)*calMorse_ECM(dist);
			}
			else{
				fMorse=calMorse_ECM(dist);
			}  
			eMorse=eMorse + calMorseEnergy_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ;
			if (_adhPairECM_Cell[i]==index) {
				fAdh=CalAdhECM(dist) ; 
				eAdh=eAdh+CalAdhEnergy(dist) ; 
				fAdhX=fAdhX+fAdh*(locX_C-locX)/dist ; 
				fAdhY=fAdhY+fAdh*(locY_C-locY)/dist ; 
			}			 
		}
	}
	
	for (int i=cellId*_maxNodePerCell ; i<(cellId+1)*_maxNodePerCell ; i++) {
		if (_nodeCellGrowProAddr[i] > _mitoticThreshold){
			enteringMitotic = true;
			// scaling = (_nodeCellGrowProAddr[i] - _mitoticThreshold)/(1 - _mitoticThreshold);
			scaling = 1; // change made on 2/2/2022
		}
		else{
			enteringMitotic = false;
			scaling = 1;
		}
		if (_nodeIsActive_Cell[i] && (i%_maxNodePerCell)<_maxMembrNodePerCell){
			locX_C=_locXAddr_Cell[i]; 
			locY_C=_locYAddr_Cell[i];
			dist=sqrt((locX-locX_C)*(locX-locX_C)+(locY-locY_C)*(locY-locY_C)) ;
			if (enteringMitotic == true){
				fMorse=(1.0+(2.0-1.0)*scaling)*calMorse_ECM(dist);
			}
			else{
				fMorse=calMorse_ECM(dist);
			}
			// fMorse=calMorse_ECM(dist);  
			eMorse=eMorse + calMorseEnergy_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ;
			if (_adhPairECM_Cell[i]==index) {
				fAdh=CalAdhECM(dist) ; 
				eAdh=eAdh+CalAdhEnergy(dist) ; 
				fAdhX=fAdhX+fAdh*(locX_C-locX)/dist ; 
				fAdhY=fAdhY+fAdh*(locY_C-locY)/dist ; 
			}			 
		}
	}

	for (int i=cellIdAfter*_maxNodePerCell ; i<(cellIdAfter+1)*_maxNodePerCell ; i++) {
		if (_nodeCellGrowProAddr[i] > _mitoticThreshold){
			enteringMitotic = true;
			// scaling = (_nodeCellGrowProAddr[i] - _mitoticThreshold)/(1 - _mitoticThreshold);
			scaling = 1; // change made on 2/2/2022
		}
		else{
			enteringMitotic = false;
			scaling = 1;
		}
		if (_nodeIsActive_Cell[i] && (i%_maxNodePerCell)<_maxMembrNodePerCell){
			locX_C=_locXAddr_Cell[i]; 
			locY_C=_locYAddr_Cell[i];
			dist=sqrt((locX-locX_C)*(locX-locX_C)+(locY-locY_C)*(locY-locY_C)) ;
			if (enteringMitotic == true){
				fMorse=(1.0+(2.0-1.0)*scaling)*calMorse_ECM(dist);
			}
			else{
				fMorse=calMorse_ECM(dist);
			}
			// fMorse=calMorse_ECM(dist);  
			eMorse=eMorse + calMorseEnergy_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ;
			if (_adhPairECM_Cell[i]==index) {
				fAdh=CalAdhECM(dist) ; 
				eAdh=eAdh+CalAdhEnergy(dist) ; 
				fAdhX=fAdhX+fAdh*(locX_C-locX)/dist ; 
				fAdhY=fAdhY+fAdh*(locY_C-locY)/dist ; 

			}			 
		}
	}
	



	return thrust::make_tuple(fTotalMorseX+fAdhX,fTotalMorseY+fAdhY,eMorse,eAdh) ;

	}

}; 




struct TotalECMForceCompute: public thrust::unary_function<DDDDDD,DD> {

	double _dummy ; 

	__host__ __device__ TotalECMForceCompute(double dummy):_dummy(dummy) {
	}

	__host__ __device__ DD operator() (const DDDDDD & dDDDDD) const {

	double fLinSpringX=  thrust:: get<0>(dDDDDD); 
	double fLinSpringY=  thrust:: get<1>(dDDDDD); 
	double fBendSpringX= thrust:: get<2>(dDDDDD); 
	double fBendSpringY= thrust:: get<3>(dDDDDD); 
	double fMembX       = thrust:: get<4>(dDDDDD); 
	double fMembY       = thrust:: get<5>(dDDDDD); 


	return thrust::make_tuple(fLinSpringX+fBendSpringX+fMembX,fLinSpringY+fBendSpringY+fMembY); 
//	return thrust::make_tuple(fLinSpringX+fMembX,fLinSpringY+fMembY); 
	}
}; 

struct RHSCompute: public thrust::unary_function<DDDDD,DD> {

	double _dt ;

	__host__ __device__ 
	          RHSCompute(double dt): _dt(dt)
	{
	}
	__host__ __device__ DD operator() (const DDDDD & dDDDD) const {
    
		double fExplicitX= thrust:: get<0>(dDDDD); 
		double fExplicitY= thrust:: get<1>(dDDDD); 
		double locX      = thrust:: get<2>(dDDDD); 
		double locY      = thrust:: get<3>(dDDDD);
		double dampCoef  = thrust:: get<4>(dDDDD);
		
		return thrust::make_tuple(fExplicitX*_dt/dampCoef + locX , 
		    	                  fExplicitY*_dt/dampCoef + locY); 

	}
	

}; 


struct TotalExplicitECMForceCompute: public thrust::unary_function<DDDD,DD> {


	__host__ __device__ TotalExplicitECMForceCompute(){
	}

	__host__ __device__ DD operator() (const DDDD & dDDD) const {

	double fBendSpringX= thrust:: get<0>(dDDD); 
	double fBendSpringY= thrust:: get<1>(dDDD); 
	double fMembX       = thrust:: get<2>(dDDD); 
	double fMembY       = thrust:: get<3>(dDDD); 


	return thrust::make_tuple(fBendSpringX+fMembX,fBendSpringY+fMembY); 
	}
}; 




struct MechProp: public thrust::unary_function<EType,DD> {


	__host__ __device__ MechProp() {
	}

	__device__ DD operator() (const EType  & nodeType) const {

		double stiffness ;
		double sponLen   ;    

		DefineECMStiffnessAndLknot (nodeType, stiffness, sponLen) ;  

		return thrust::make_tuple(stiffness,sponLen); 
	}
}; 



struct MoveNodesECM: public thrust::unary_function<DDDDD,DD> {

	double _dt ; 
	__host__ __device__ MoveNodesECM (double dt): _dt(dt) {
	}
	__host__ __device__ DD operator() (const DDDDD & dDDDD) const  {
		double 			  locXOld= thrust:: get <0> (dDDDD) ;
		double 			  locYOld= thrust:: get <1> (dDDDD) ;
		double 			  fx= 	   thrust:: get <2> (dDDDD) ;
		double 			  fy= 	   thrust:: get <3> (dDDDD) ;
		double            dampCoef=thrust:: get <4> (dDDDD) ; 
		return thrust::make_tuple (locXOld+fx*_dt/dampCoef, locYOld+fy*_dt/dampCoef) ;
	}
}; 


struct CalBendECM: public thrust::unary_function<IDD, DDDDDD> {
	double* _locXAddr;
	double* _locYAddr;
	int  _numECMNodes ;
	double _eCMBendStiff ;  
	// bool* _isActiveECM ;

	__host__ __device__ CalBendECM (double* locXAddr, double* locYAddr, int numECMNodes, double eCMBendStiff/*, bool* isActiveECM*/) :
				_locXAddr(locXAddr), _locYAddr(locYAddr),_numECMNodes(numECMNodes),_eCMBendStiff(eCMBendStiff)/*, _isActiveECM(isActiveECM)*/{
	}
	
	__host__ __device__ DDDDDD operator()(const IDD &  iDD) const {

		int   nodeRank = thrust::get<0>(iDD);
		double locX = thrust::get<1>(iDD);
		double locY = thrust::get<2>(iDD);

		double kb=_eCMBendStiff ; 
		double leftPosX,leftPosY;
		double lenLeft;

		double rightPosX,rightPosY;
		double lenRight;
		double cosTheta, bendLeftX, bendRightX, bendLeftY, bendRightY, bendCenterX, bendCenterY;
		bool activeTrio;
			
			int index_left = nodeRank - 1;
			if (index_left == -1) {
				index_left = _numECMNodes - 1;
			}
		
				leftPosX = _locXAddr[index_left];
				leftPosY = _locYAddr[index_left];
				lenLeft = sqrt((leftPosX - locX) * (leftPosX - locX) + (leftPosY - locY) * (leftPosY - locY) );


			int index_right = nodeRank + 1;
			if (index_right ==  _numECMNodes) {
				index_right = 0;
			}
				rightPosX = _locXAddr[index_right];
				rightPosY = _locYAddr[index_right];
				lenRight = sqrt((rightPosX - locX) * (rightPosX - locX) + (rightPosY - locY) * (rightPosY - locY) );

			// if (activeTrio == true){
				cosTheta=( (leftPosX-locX)*(rightPosX-locX)+(leftPosY-locY)*(rightPosY-locY) )/(lenRight*lenLeft) ; 

				
				bendLeftX= -kb*(rightPosX-locX)/(lenRight*lenLeft)+
										kb*cosTheta/(lenLeft*lenLeft)*(leftPosX-locX) ; 

				bendRightX=-kb*(leftPosX-locX)/(lenRight*lenLeft)+
										kb*cosTheta/(lenRight*lenRight)*(rightPosX-locX) ; 

				bendLeftY =-kb*(rightPosY-locY)/(lenRight*lenLeft)+
										kb* cosTheta/(lenLeft*lenLeft)*(leftPosY-locY) ; 

				bendRightY=-kb*(leftPosY-locY)/(lenRight*lenLeft)+
										kb* cosTheta/(lenRight*lenRight)*(rightPosY-locY) ;

				bendCenterX=-(bendLeftX+bendRightX) ; 
				bendCenterY=-(bendLeftY+bendRightY) ; 
			// }	
			return thrust::make_tuple(bendCenterX,bendCenterY,
					bendLeftX, bendLeftY, bendRightX, bendRightY);

	}
}; 


struct SumBendForce: public thrust::unary_function<IDD,DD> {

	double* _fBendLeftXAddr;
	double* _fBendLeftYAddr;
	double* _fBendRightXAddr;
	double* _fBendRightYAddr;
	int  	_numECMNodes ;

	__host__ __device__ SumBendForce (double* fBendLeftXAddr, double* fBendLeftYAddr, double *fBendRightXAddr, double *fBendRightYAddr, int numECMNodes) :
				_fBendLeftXAddr(fBendLeftXAddr), _fBendLeftYAddr(fBendLeftYAddr),_fBendRightXAddr(fBendRightXAddr),_fBendRightYAddr(fBendRightYAddr),_numECMNodes(numECMNodes){
	}
	
	__host__ __device__ DD operator() (const IDD & iDD) const {

		int   nodeRank = thrust::get<0>(iDD);
		double fBendCenterX = thrust::get<1>(iDD);
		double fBendCenterY = thrust::get<2>(iDD);

		int index_left = nodeRank - 1;
		if (index_left == -1) {
			index_left = _numECMNodes - 1;
		}
		
		int index_right = nodeRank + 1;
		if (index_right ==  _numECMNodes) {
			index_right = 0;
		}

	return thrust::make_tuple(fBendCenterX+_fBendLeftXAddr[index_right]+_fBendRightXAddr[index_left],
				  fBendCenterY+_fBendLeftYAddr[index_right]+_fBendRightYAddr[index_left]); 
	}
}; 




struct AssignDamping : public thrust::unary_function<EType,double> {


   double _dampBasal ;
   double _dampBC ; 
   double _dampApical ; 

  __host__ __device__  AssignDamping  (double dampBasal, double dampBC, double dampApical ):
   _dampBasal (dampBasal) , _dampBC (dampBC), _dampApical (dampApical) {
					 }

	__host__ __device__  double operator()(const EType & eCMNodeType) const {

		if ( eCMNodeType==excm ){

			return (_dampBasal) ; 
		}
			else if (eCMNodeType==bc2 ) {

				return (_dampBC) ; 
			}
				else if ( eCMNodeType==perip ){

					return (_dampApical) ; 
				}
					else {
				 		return (0.0) ; // this is an indication of error
					}
	}
	
} ; 


struct ECMGrowFunc: public thrust::unary_function<DT, bool> {
	uint _currECMCount;
	uint _bound;
	double _maxLengthToAddECMNodes;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight__host__ __device__
	__host__ __device__ ECMGrowFunc(uint currECMCount, uint bound, double maxLengthToAddECMNodes) :
	_currECMCount(currECMCount), _bound(bound), _maxLengthToAddECMNodes(maxLengthToAddECMNodes) {
	}
	//Ali __host__ __device__ BoolD operator()(const DUi& dui) {
	__host__ __device__ bool operator()(const DT& dt) {
		//Ali double progress = thrust::get<0>(dui);
        double LengthMax=thrust::get<0>(dt); //Ali
		EType newNodeType = thrust::get<1>(dt);
		if (_currECMCount < _bound   && LengthMax>_maxLengthToAddECMNodes){
			return true;
		}
		else {
			return false;
		}
		
	}
};

#endif
