//
// Created by Sahand Adibnia on 3/31/24.
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
#include "Position_Gradient.h"

#ifndef HW5_FINITE_DIFFERENCES_H
#define HW5_FINITE_DIFFERENCES_H

// Code does not work at the moment, give a systematically incorrect gradient

double central_differences(int deriv_atom_index, const std::vector<Atom> &atoms, double h, int direction);
arma::vec central_diff_grad(int deriv_atom_index, const std::vector<Atom> &atoms, double h);
arma::mat central_diff_mat(const std::vector<Atom> &atoms, double h);

#endif //HW5_FINITE_DIFFERENCES_H
