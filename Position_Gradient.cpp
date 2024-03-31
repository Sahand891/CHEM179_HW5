//
// Created by Sahand Adibnia on 3/30/24.
//

#include "Position_Gradient.h"

// Overlap Derivative

double analytical_S_1D_deriv(double X_a, double alpha_a, int l_a, double X_b, double alpha_b, double l_b) {

    double term1 = -l_a*analytical_overlap_int_1D(X_a, alpha_a, l_a-1, X_b, alpha_b, l_b);
    double term2 = 2*alpha_a*analytical_overlap_int_1D(X_a, alpha_a, l_a+1, X_b, alpha_b, l_b);

    return term1+term2;
};

double S_deriv_one_element(AO &AO_1, AO &AO_2, int direction) {

    int k=0; // first AO primG index
    double sum=0;

    // If the two AOs exist on the same atom, the derivative of the overlap integral when moving that atom (or any other atom) is 0!
    if (AO_1.origin_atom.X == AO_2.origin_atom.X and AO_1.origin_atom.Y == AO_2.origin_atom.Y and AO_1.origin_atom.Z == AO_2.origin_atom.Z and AO_1.origin_atom.A == AO_2.origin_atom.A) {
        return 0;
    }

    for (auto& primG_1 : AO_1.cGTO.primGs) {
        int l=0; // second AO primG index
        for (auto& primG_2 : AO_2.cGTO.primGs) {
            double coef = AO_1.cGTO.d_primes[k] * AO_2.cGTO.d_primes[l];
            double I_x = analytical_overlap_int_1D(AO_1.origin_atom.X, primG_1.alpha, primG_1.l, AO_2.origin_atom.X, primG_2.alpha, primG_2.l);
            double I_y = analytical_overlap_int_1D(AO_1.origin_atom.Y, primG_1.alpha, primG_1.m, AO_2.origin_atom.Y, primG_2.alpha, primG_2.m);
            double I_z = analytical_overlap_int_1D(AO_1.origin_atom.Z, primG_1.alpha, primG_1.n, AO_2.origin_atom.Z, primG_2.alpha, primG_2.n);
            if (direction == 0) {
                double I_deriv = analytical_S_1D_deriv(AO_1.origin_atom.X, primG_1.alpha, primG_1.l, AO_2.origin_atom.X, primG_2.alpha, primG_2.l);
                sum += coef*I_y*I_z*I_deriv;
            } else if (direction == 1) {
                double I_deriv = analytical_S_1D_deriv(AO_1.origin_atom.Y, primG_1.alpha, primG_1.m, AO_2.origin_atom.Y, primG_2.alpha, primG_2.m);
                sum += coef*I_x*I_z*I_deriv;
            } else {
                double I_deriv = analytical_S_1D_deriv(AO_1.origin_atom.Z, primG_1.alpha, primG_1.n, AO_2.origin_atom.Z, primG_2.alpha, primG_2.n);
                sum += coef*I_x*I_y*I_deriv;
            }
            l += 1;
        }
        k += 1;
    }


    return sum;

}

arma::mat S_deriv(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;
    std::vector<AO> AOs = conv_it_data.AOs;

    // Initialize appropriately sized matrix - each row is a dimension for a certain atom, and each column is an AO
    arma::mat final_mat(3, pow(AOs.size(),2), arma::fill::zeros);

    for (size_t i = 0; i < final_mat.n_rows; i++) { // iterate through directions x, y, z
        int col_index=0;
        for (auto& ao_1 : AOs) {
            // Fix one AO as AO_1, then iterate through rest of AOs as AO_2
            for (auto& ao_2: AOs) {
                final_mat(i,col_index) = S_deriv_one_element(ao_1, ao_2, i);
                col_index += 1;
            }
        }
    }

    return final_mat;
}



// Gamma Derivative - based on s orbitals only

double zero_zero_deriv(Atom A, Atom B, int k, int k_prime, int l, int l_prime, int direction) {

    // Calculating d = (R_A - R_B)**2
    arma::vec d_vec = A.pos_vec - B.pos_vec;
    double d = pow(arma::norm(d_vec),2);

    if (d == 0) { // if we're looking at the derivative of gammaAA, essentially
        return 0;
    }

    // Only care about valence s orbital, so always take the FIRST cGTO
    double alpha_k_A = A.cGTOs[0].primGs[k].alpha;
    double alpha_k_prime_A = A.cGTOs[0].primGs[k_prime].alpha;
    double alpha_k_B = B.cGTOs[0].primGs[l].alpha;
    double alpha_k_prime_B = B.cGTOs[0].primGs[l_prime].alpha;

    double sigma_A = 1 / (alpha_k_A + alpha_k_prime_A);
    double U_A = pow((M_PI * sigma_A), 1.5);

    double sigma_B = 1 / (alpha_k_B + alpha_k_prime_B);
    double U_B = pow((M_PI * sigma_B), 1.5);

    double V_squared = 1 / (sigma_A + sigma_B);
    double T = V_squared * d;


    // Compiling individual components of the derivative

    double term1 = U_A*U_B*(A.X-B.X) / d;
    if (direction == 1) {
        term1 = U_A*U_B*(A.Y-B.Y) / d;
    } else if (direction == 2) {
        term1 = U_A*U_B*(A.Z-B.Z) / d;
    }
    double term2 = -erf(sqrt(T)) / sqrt(d);
    double term3 = 2*sqrt(V_squared)*exp(-T) / sqrt(M_PI);

    return term1*(term2+term3) * 27.2114079527; // multiply by 27.211 to convert Hartree to eV (1 Hartree = 27.2114079527 eV)

}

double gamma_deriv_one_element(Atom &A, Atom &B, int direction) {

    double sum=0;

    for (int k = 0; k < 3; k++) {
        double d_k_A = A.cGTOs[0].d_primes[k];
        for (int k_prime = 0; k_prime < 3; k_prime++) {
            double d_kprime_A = A.cGTOs[0].d_primes[k_prime];
            for (int l = 0; l < 3; l++) {
                double d_l_B = B.cGTOs[0].d_primes[l];
                for (int l_prime = 0; l_prime < 3; l_prime++) {
                    double d_lprime_B = B.cGTOs[0].d_primes[l_prime];

                    double deriv = zero_zero_deriv(A, B, k, k_prime, l, l_prime, direction);

                    sum += d_k_A * d_kprime_A * d_l_B * d_lprime_B * deriv;

                }
            }
        }
    }


    return sum;

}

arma::mat gamma_deriv(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    // Initialize appropriately sized matrix - each row is a dimension for a certain atom, and each column is an AO
    arma::mat final_mat(3, pow(atoms.size(),2), arma::fill::zeros);

    for (size_t i = 0; i < final_mat.n_rows; i++) { // iterate through directions x, y, z
        int col_index=0;
        for (auto& atom1 : atoms) {
            // Fix one atom as atom1, then iterate through rest of atoms as atom2
            for (auto& atom2: atoms) {
                final_mat(i,col_index) = gamma_deriv_one_element(atom1, atom2, i);
                col_index += 1;
            }
        }
    }

    return final_mat;

}











// Nuclear Repulsion Derivative

double calc_V_nuc_deriv() {


    return 0;
}


arma::mat V_nuc_deriv(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    // Initialize appropriately sized matrix - each row is a dimension, each column is an atom
    arma::mat final_mat(3, atoms.size(), arma::fill::zeros);


    return final_mat;
}
