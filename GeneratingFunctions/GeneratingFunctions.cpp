
/*
* 	(C) Peter Mann (2018) [University of St Andrews]
* 	Computes fraction of infected nodes as a function of the disease
* 	transmissibility for the SIR process.
*
*	$ g++ -std=c++11 auxillary.cpp SIR.cpp -o out
*	$ ./out
*
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "auxillary.h"   // include auxillary functions declaration

#define N 10000          // Number of steps
#define C 3              // Average degree
#define Accuracy 1e-10   // Accuracy required

class GF {
	/* Experiment class to compute generating functions for an SIR process. 
	   Class takes the degree distribution as its argument as a vector
	   (see auxillary file for further info).
	*/

	std::vector<double> _Pk;
	double u, T;
	int i;

	public:
		// constructor for class, takes degree distribution
		GF (std::vector<double> Pk): _Pk(Pk) {
		}
 
 		// function to compute average degree
		double averageDegree (const std::vector<double> &_Pk);

		// function to compute G_0
		double G_0 ( const std::vector<double> &_Pk, double ave_k, double z );

		// function to compute G_1
		double G_1 ( const std::vector<double> &_Pk, double ave_k, double z );

		// function to compute infected fraction
		double evaluate( const std::vector<double> &_Pk, double ave_k, double T );

		// function runs the experiment
		void _do () {
			// open an output file
			std::ofstream output_file ("fraction_infected.txt");
			if (output_file.is_open()){
  				// compute the average degree
				double ave_k = averageDegree(_Pk);

				// iterate through updates
				for (i=1; i<=N; i++){

					// compute transmissibility
					T = (double) i/N;

					// compute infected fraction & print to file
					output_file << T << " " << evaluate(_Pk, ave_k, T) << "\n";
				}
				// close output file after use
				output_file.close();
			}
			// if can't open output file ... report it!
			else std::cout << "Unable to open output file\n ";
		}

};

double GF::evaluate ( const std::vector<double> &Pk, double ave_k, double T ) {
/*
 * Class: GF
 * Function:  evaluate
 * --------------------
 * Computes the number of infected nodes following an epidemic:
 * 
 *    				S = 1 - G_0(1-T+Tu)
 *  
 * after first solving the self-consistent equation:
 * 
 * 				  u = G_1(1-T+Tu)
 * 
 * The parameter `u` is the probability that a neighbour remains
 * susceptible during the disease spreading. If there is an 
 * outbreak of the disease, it will have a unique solution in the 
 * unit interval. 
 *
 * arg Pk: degree distribution (vector)
 * arg ave_k: average degree (double)
 * arg T: transmissibility (double)
 * 
 * returns: fraction of infected nodes (double)
 */
	double u_old; // refresh variable for old u value
	u = 0.5;      // initialise u

	// solves self-consistent equation
	do {
		// update old value
		u_old = u; 

		// iterate new value
    		u = G_1(Pk, ave_k, 1 - T + T*u); 

    		// while difference is above accuracy
  	} while (fabs(u-u_old)>Accuracy);

  	// compute fraction of infected nodes
  	return 1 - G_0(Pk, ave_k, 1 - T + T*u);
}



double GF::G_0 ( const std::vector<double> &Pk, double ave_k, double z ) {
/*
 * Class: GF
 * Function:  G_0
 * --------------------
 * Computes the G_0 generating function given an 
 * average degree and the degree distribution.
 *
 * arg Pk: degree distribution (vector)
 * arg ave_k: average degree (double)
 * arg z: the generating function argument (double)
 *
 * returns: \sum_k p_k * z**k
 */
 	double sum = 0.0;
 	for (auto iter = Pk.begin(); iter != Pk.end(); ++iter) {
    		int k = std::distance(Pk.begin(), iter);
    		sum += Pk.at(k) * pow(z, k);
	}
	return sum;
}

double GF::G_1 ( const std::vector<double> &Pk, double ave_k, double z ) {
/*
 * Class: GF
 * Function:  G_1
 * --------------------
 * Computes the G_1 generating function given an 
 * average degree and the degree distribution.
 *
 * arg Pk: degree distribution (vector)
 * arg ave_k: average degree (double)
 * arg z: the generating function argument (double)
 *
 * returns: \frac{1}{\langle k \rangle}\sum_{k=0}^\infty (k+1)p_{k+1}(x)^k
 */
 	double sum = 0.0;
 	for (auto iter = Pk.begin(); iter != Pk.end(); ++iter) {
    		int k = std::distance(Pk.begin(), iter);
    		sum += k * Pk.at(k) * pow(z, k-1);
	}
	return (double) sum / ave_k;
}

double GF::averageDegree ( const std::vector<double> &Pk ) {
/* 
 * Class: GF
 * Function:  averageDegree
 * --------------------
 * Computes the average degree by reading the
 * degree distribution from Pk (vector) at index
 * k.  
 *
 * arg Pk: the degree distribution (vector)
 *
 * returns: average degree (double)
 */	
 	double sum = 0.0;
 	for (auto iter = Pk.begin(); iter != Pk.end(); ++iter) {
    		int k = std::distance(Pk.begin(), iter);
    		sum += k * Pk.at(k);
	}
	return sum;
}


int main(int argc, char **argv){

	// compute Pk
  	std::vector<double> Pk = degreeDist();

	// construct class
	GF experiment(Pk);

	// run class member function
	experiment._do();

}



