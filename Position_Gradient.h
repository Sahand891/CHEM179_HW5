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

#endif //HW5_POSITION_GRADIENT_H
