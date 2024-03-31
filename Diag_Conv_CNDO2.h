//
// Created by Sahand Adibnia on 3/30/24.
//

#ifndef HW5_DIAG_CONV_CNDO2_H
#define HW5_DIAG_CONV_CNDO2_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>
#include "Fock_Matrix.h"
#include "Data_Abstraction.h"


struct iteration_data {

    std::vector<Atom> atoms;

    int p;
    int q;

    arma::mat P_alpha_old;
    arma::mat P_beta_old;

    arma::mat fock_alpha;
    arma::mat fock_beta;

    arma::mat C_alpha_new;
    arma::mat C_beta_new;

    arma::mat P_alpha_new;
    arma::mat P_beta_new;

    int iteration_count;

    // check for convergence just based on alpha density matrix
    bool converged = arma::approx_equal(P_alpha_new,P_alpha_old,"absdiff",1e-6);

    // Other information in case we need it
    std::vector<AO> AOs = atoms_to_AOs(atoms);
    arma::vec P_total_new = P_AA(atoms, P_alpha_new, P_beta_new);

    arma::mat gamma = gamma_matrix(atoms);
    arma::mat overlap = S(atoms);
    arma::mat H_core = h(atoms);



};


arma::mat find_MO_coefs(arma::mat F);
arma::vec get_energy_eigs(arma::mat F);


iteration_data start_CNDO2(const std::vector<Atom> &atoms, bool do_print=false, bool auto_electrons=true, int p_input=0, int q_input=0);
// Stores a vector of iteration data, which can then be put in an output file nicely!!
std::vector<iteration_data> converge_CNDO2(const iteration_data &it_data, std::vector<iteration_data> &it_data_vec, bool do_print=false);
iteration_data obtain_converged_data(std::vector<iteration_data> &converged_it_data);

double nuc_repl_energy(const std::vector<Atom> &atoms);
double electron_energy(const arma::mat &F_alpha, const arma::mat &F_beta, const arma::mat &P_alpha, const arma::mat &P_beta, const arma::mat &H_core);
double compute_total_energy(const iteration_data &conv_it_data);

// Function that outputs converged CNDO/2 data
iteration_data perform_CNDO2(const std::vector<Atom> &atoms);


#endif //HW5_DIAG_CONV_CNDO2_H
