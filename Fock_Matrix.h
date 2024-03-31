//
// Created by Sahand Adibnia on 3/30/24.
//

#ifndef HW5_FOCK_MATRIX_H
#define HW5_FOCK_MATRIX_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Data_Abstraction.h"

// Return matrix elements and matrix for gamma
double gamma_AB(Atom A, Atom B);
arma::mat gamma_matrix(std::vector<Atom> atoms);

// Compute overlap matrix from a vector of atoms
arma::mat S(const std::vector<Atom> &Atoms);

// Compute a core, 1-electron Hamiltonian
arma::mat h(const std::vector<Atom> &atoms);


// Function that takes in MO coefficient matrix for alpha, a number of alpha electrons, returns the rectangular occupied coefficient matrix
arma::mat C_occ_alpha(arma::mat c_alpha, int p);
// Calculates density matrix for alpha electrons from occupied coefficient matrix for alpha electrons
arma::mat P_alpha(arma::mat c_occ_alpha);

// Returns vector of total electron density on each atom
arma::vec P_AA(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta);

// Fock matrix for alpha electrons
arma::mat F_alpha(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta);

#endif //HW5_FOCK_MATRIX_H
