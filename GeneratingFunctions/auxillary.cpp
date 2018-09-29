
/*
* (C) Peter Mann (2018) [University of St Andrews]
* Helper function(s) in auxillary file.
*/

#include <iostream>
#include <fstream>
#include <vector>

#include "auxillary.h"

std::vector<double> degreeDist(){
/*
 * Function:  degreeDist
 * --------------------
 * Computes the degree distribution by reading it from 
 * a file named input_data.txt.
 *
 * returns: degree distribution (array, double)
 */	
 	// create a vector
 	std::vector<double> Pk;

 	// open input file
 	std::ifstream input_file("input_data.txt", std::ios::in);
 	if (input_file.is_open()){
 		// while values in file, grab values and populate Pk
 		double num;
 		while (input_file >> num){
 			Pk.push_back(num);
 		}

 		// close the file stream
		input_file.close();

		// return degree distribution vector
		return Pk;
 	}
 	else{
 		std::cout << "error opening input data file!\n ";
 		exit(1);
 	}	
 }