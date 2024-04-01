//
// Created by Sahand Adibnia on 3/30/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Fock_Matrix.h"
#include "Data_Abstraction.h"
#include "Diag_Conv_CNDO2.h"

#ifndef HW5_POSITION_GRADIENT_H
#define HW5_POSITION_GRADIENT_H

arma::mat S_deriv(const iteration_data &conv_it_data);

arma::mat gamma_deriv(const iteration_data &conv_it_data);

arma::mat V_nuc_deriv(const iteration_data &conv_it_data);

// Calculation of energy gradient
arma::mat x_mat(const iteration_data &conv_it_data);
arma::mat y_mat(const iteration_data &conv_it_data);

double electron_gradient_x(const iteration_data &conv_it_data, int deriv_atom_index, int direction);
arma::mat electron_gradient_mat(const iteration_data &conv_it_data);

arma::mat total_gradient_mat(const iteration_data &conv_it_data);


void write_gradient_to_file(const iteration_data &conv_it_data, std::string output_location);

#endif //HW5_POSITION_GRADIENT_H
