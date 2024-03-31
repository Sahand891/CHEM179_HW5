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


// Overlap Matrix
double S_uv(contr_Gaussian contrG1, contr_Gaussian contrG2) {

    double final_overlap = 0;
    for (int i = 0; i < 3; i++) {
        double norm_i = contrG1.norm_constants[i];
        double contr_coef_i = contrG1.contr_coefs[i];
        for (int j = 0; j < 3; j++) {
            double overlap = S_ab(contrG1.primGs[i], contrG2.primGs[j]);
            double norm_j = contrG2.norm_constants[j];
            double contr_coef_j = contrG2.contr_coefs[j];
            final_overlap += contr_coef_i * contr_coef_j * norm_i * norm_j * overlap;
        }
    }
    return final_overlap;

}


arma::mat S(const std::vector<Atom> &Atoms) {

    std::vector<AO> AOs = atoms_to_AOs(Atoms);
    int matrix_size = AOs.size();

    // Initialize a matrix of the appropriate size, just filled with zeroes
    arma::mat S(matrix_size, matrix_size, arma::fill::zeros);

    // Iterate through indices in matrix, computing appropriate overlap integral!
    for (size_t i = 0; i < matrix_size; i++) {
        //std::cout<< i<< std::endl;
        for (size_t j = 0; j < matrix_size; j++) {
            double Sij = S_uv(AOs[i].cGTO, AOs[j].cGTO);
            S(i, j) = Sij;
            //std::cout<< i <<j<< std::endl;
            //S.print();
        }
    }

    return S;
}

// One-electron Hamiltonian
arma::mat h(const std::vector<Atom> &atoms) {

    std::vector<AO> AOs = atoms_to_AOs(atoms);

    // initializing appropriately sized matrix
    arma::mat final_h(AOs.size(), AOs.size(), arma::fill::zeros);

    // Diagonal matrix elements
    double term1, term2, term3;
    int i=0; // to keep the index
    for (auto& AO : AOs) {
        // Reset values of all terms on each iteration
        term1=0, term2=0, term3=0;

        Atom atom = AO.origin_atom; // atom corresponding to the AO

        // term2, term3 the same—only dependent on atom, not also the orbital
        term2 = (atom.Z_ - 0.5)*gamma_AB(atom, atom);
        for (auto& atom2 : atoms) {
            term3 += atom2.Z_*gamma_AB(atom,atom2);
        }
        term3 -= atom.Z_*gamma_AB(atom,atom); // get rid of the unnecessary term from for loop above

        // term1 is a bit more complicated, we need to care about the orbital angular momentum too (s or p)
        if (AO.L == 0) { // the AO we're on is an s AO
            term1 = atom.IA_s;
        } else if (AO.L == 1) { // the AO we're on is a p AO (L=1)
            term1 = atom.IA_p;
        }

//        std::cout << term1 << std::endl;
//        std::cout << term2 << std::endl;
//        std::cout << term3 << std::endl;


        final_h(i,i) = -term1-term2-term3;

        i += 1;
    }


    // Now the off-diagonal matrix elements
    for (int u=0; u < AOs.size(); u++) {
        AO AO_1 = AOs[u];

        for (int v=0; v < AOs.size(); v++) {
            AO AO_2 = AOs[v];
            if (u != v) { // off-diagonals only
                double S_ele = S_uv(AO_1.cGTO,AO_2.cGTO);
                final_h(u,v) = 0.5*(AO_1.origin_atom.beta + AO_2.origin_atom.beta)*S_ele;
            }
        }
    }

    return final_h;

}











// Density Matrix

// Let's just try calculating the density the matrix more simply from the OCCUPIED coefficient matrices:
// Getting the occupied coefficient matrix from a total coefficient matrix
arma::mat C_occ_alpha(arma::mat c_alpha, int p) {
    // p is the number of alpha electrons here

    int num_rows = sqrt(c_alpha.size()); // since the size of a matrix is the total # of elements, and c_alpha is always square

    // Initialize appropriately sized matrix
    arma::mat final_mat(num_rows, p, arma::fill::zeros);

    // slice the c_alpha square matrix into a rectangular matrix, just selecting the first p columns
    for (size_t j=0; j < p; j++) {
        for (size_t i=0; i < num_rows; i++) {
            final_mat(i,j) = c_alpha(i,j);
        }
    }
    return final_mat;
}

arma::mat P_alpha(arma::mat c_occ_alpha) {
    // Simply the occupied coefficient matrix times its transpose!
    return c_occ_alpha*c_occ_alpha.t();
}


// Total electron density on one atom—the sum of diagonal total matrix elements
arma::vec P_AA(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta){

    // First get the total density matrix
    arma::mat P_tot = P_alpha + P_beta;

    // Initialize an arma::vec of the appropriate size
    arma::vec final_vec(atoms.size(), arma::fill::zeros);

    double sum=0;
    // Assumes order of "atoms" vector corresponds to density matrix AO order
    int i =0;
    while (i < atoms.size()) {
        if (atoms[i].A == 1) { // if you're dealing with a Hydrogen atom as the atom you're at
            final_vec(i) = P_tot(i,i);
            i += 1;
        } else { // if it's not a hydrogen atom, it's gonna be a 2nd row element atom (C, N, O, F)
            // Now we take the sum of the next FOUR elements
            final_vec(i) = P_tot(i,i) + P_tot(i+1,i+1) + P_tot(i+2,i+2) + P_tot(i+3,i+3);
            i += 4;
        }
    }

    return final_vec;
}






// Gamma Calculation

double zero_zero_integral(Atom A, Atom B, int k, int k_prime, int l, int l_prime) {

    double final_val = 0;

    // Only care about valence s orbital, so always take the FIRST cGTO
    double alpha_k_A = A.cGTOs[0].primGs[k].alpha;
    double alpha_k_prime_A = A.cGTOs[0].primGs[k_prime].alpha;
    double alpha_k_B = B.cGTOs[0].primGs[l].alpha;
    double alpha_k_prime_B = B.cGTOs[0].primGs[l_prime].alpha;

    double sigma_A = 1 / (alpha_k_A + alpha_k_prime_A);
    double U_A = pow((M_PI * sigma_A), 1.5);

    double sigma_B = 1 / (alpha_k_B + alpha_k_prime_B);
    double U_B = pow((M_PI * sigma_B), 1.5);

    // Calculating d = (R_A - R_B)**2
    arma::vec d_vec = {(A.X - B.X), (A.Y - B.Y), (A.Z - B.Z)};
    double d = pow(arma::norm(d_vec),2);

    double V_squared = 1 / (sigma_A + sigma_B);
    double T = V_squared * d;


    if (T == 0) {
        final_val = U_A * U_B * sqrt(2 * V_squared) * sqrt(2.0 / M_PI);
    } else {
        final_val = U_A * U_B * sqrt(1.0 / d) * erf(sqrt(T));
    }

    return 27.2114079527*final_val; // convert gamma elements to eV from atomic units

}

double gamma_AB(Atom A, Atom B) {

    // Note: only looks at valence s orbitals of A and B
    // This means we're looking at the FIRST cGTO ALWAYS
    // but at each index i, take the ith primitive Gaussian's contraction coefficient and normalization constant!

    double sum = 0;

    int count = 0;
    for (int k = 0; k < 3; k++) {
        double d_k_A = A.cGTOs[0].contr_coefs[k] * A.cGTOs[0].norm_constants[k];
        for (int k_prime = 0; k_prime < 3; k_prime++) {
            double d_kprime_A = A.cGTOs[0].contr_coefs[k_prime] * A.cGTOs[0].norm_constants[k_prime];
            for (int l = 0; l < 3; l++) {
                double d_l_B = B.cGTOs[0].contr_coefs[l] * B.cGTOs[0].norm_constants[l];
                for (int l_prime = 0; l_prime < 3; l_prime++) {
                    double d_lprime_B = B.cGTOs[0].contr_coefs[l_prime] * B.cGTOs[0].norm_constants[l_prime];

                    double integral = zero_zero_integral(A, B, k, k_prime, l, l_prime);

                    sum += d_k_A * d_kprime_A * d_l_B * d_lprime_B * integral;
                    count += 1;

                }
            }
        }
    }

    //std::cout << count << std::endl;

    return sum;
}

arma::mat gamma_matrix(std::vector<Atom> atoms) {

    int matrix_size = atoms.size();
    arma::mat final_mat(matrix_size, matrix_size, arma::fill::zeros);

    // Iterate through indices in matrix, computing appropriate value of gamma!
    for (size_t i = 0; i < matrix_size; i++) {
        for (size_t j = 0; j < matrix_size; j++) {
            double mat_ij = gamma_AB(atoms[i], atoms[j]);
            final_mat(i, j) = mat_ij;
        }
    }

    return final_mat;

}







// Fock matrix construction

arma::mat F_alpha(const std::vector<Atom> &atoms, arma::mat P_alpha, arma::mat P_beta) {

    std::vector<AO> AOs = atoms_to_AOs(atoms);

    // initializing appropriately sized matrix
    arma::mat final_F(AOs.size(), AOs.size(), arma::fill::zeros);

    // Diagonal matrix elements
    double term1, term2, term3;
    int i=0; // to keep the index of the AOs (which is same as row/column index of diagonal elements of F matrix)
    int j=0; // keep track of atom index
    for (auto& AO : AOs) {
        // Reset values of all terms on each iteration
        term1=0, term2=0, term3=0;

        Atom atom = AO.origin_atom; // atom corresponding to the AO

        arma::mat P_total = P_AA(atoms, P_alpha, P_beta);

        // term2, term3 only depend on the AO and the atom, not the L of the AO

        // P_AA is a VECTOR of electron density per ATOM ... so the question is, which NUMBER atom are we on??
        term2 = ((P_total[j] - atom.Z_) - (P_alpha(i,i) - 0.5))*gamma_AB(atom,atom); // indexing problem is with P_AA - which number atom are we on???


        // For calculating term 3, need this extra k index, which tells me which atom in the atoms vector that atom2 is (atom2's index), which tells me what value from P_AA to extract!
        for (int k = 0; k < atoms.size(); k++) {
            Atom atom2 = atoms[k];
            term3 += (P_total[k] - atom2.Z_)*gamma_AB(atom,atom2);
        }
        term3 -= (P_total[j] - atom.Z_)*gamma_AB(atom,atom);; // get rid of the unnecessary term from for loop above


        // term1 is a bit more complicated, we need to care about the orbital angular momentum too (s or p)
        if (AO.L == 0) { // the AO we're on is an s AO
            term1 = atom.IA_s;
        } else if (AO.L == 1) { // the AO we're on is a p AO (L=1)
            term1 = atom.IA_p;
        }



        final_F(i,i) = -term1+term2+term3;

        // j index tracks the atom we're on
        // if we just looked at H atom, then we go to the NEXT atom (b/c H will only have 1 AO)
        // if we just looked at C,N,O,F, then we go to the next AO, but stay at the same atom
        if (atom.A == 1) {
            j += 1;
        }

        // i index tracks the AO we're on
        i += 1;
    }


    // Now the off-diagonal matrix elements
    for (int u=0; u < AOs.size(); u++) {
        AO AO_1 = AOs[u];
        Atom atom_1 = AO_1.origin_atom;

        for (int v=0; v < AOs.size(); v++) {
            AO AO_2 = AOs[v];
            Atom atom_2 = AO_2.origin_atom;

            if (u != v) { // off-diagonals only
                double S_ele = S_uv(AO_1.cGTO,AO_2.cGTO);
                final_F(u,v) = 0.5*(AO_1.origin_atom.beta + AO_2.origin_atom.beta)*S_ele - P_alpha(u,v)*gamma_AB(atom_1, atom_2);
            }
        }
    }

    return final_F;

}
