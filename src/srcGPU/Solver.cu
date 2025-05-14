#include "Solver.h"

using namespace std ; 

 vector<double>  Solver::SOR3DiagPeriodic(const vector <bool> & nodeIsActive, const vector <double> & lDiag, const vector <double> & diag, const vector <double> & uDiag,
	                                      const vector <double> & rHS, 
										  const vector <int> & prevIndex,
										  const vector <int> & nextIndex,
										  vector<double>  & firstGuess) {
				
	const int maxIteration=2500 ; 
	const double beta=1.2 ; 
	
	vector <double> ans =firstGuess ;
	vector<double> ansOld=firstGuess ; 
	double maxError=10000 ;  // Just a big number to go inside while loop for the fist time. 
    const int N     =diag.size();        // --- Size of the linear system
	
	//Simple iterative SOR solver
	int numIterator=0 ; 
	while (maxError>1.0E-6 && numIterator<maxIteration) {
		maxError=0 ; 
		numIterator ++ ; 
		for (int i=0; i<N ; i++) {
			if (!nodeIsActive.at(i)) {
				continue ; 
			}
			ansOld.at(i)=ans.at(i) ; 
			ans.at(i)=beta*(rHS.at(i)-lDiag.at(i)*ans[prevIndex.at(i)]-uDiag.at(i)*ans[nextIndex.at(i)])/diag.at(i)+ (1-beta)*ansOld.at(i); 
			if ( abs (ansOld.at(i)-ans.at(i)) >maxError) {
				maxError=abs (ansOld.at(i)-ans.at(i));  
			}
		}
	}
	cout << "In SOR solver after " << numIterator <<" iteration, maximum difference in two successuve iterations is "<< maxError << endl ;  ; 
	return ans; 
 }
