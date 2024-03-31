//
// Created by Sahand Adibnia on 3/30/24.
//

#include "Position_Gradient.h"

double analytical_S_1D_deriv(double X_a, double alpha_a, int l_a, double X_b, double alpha_b, double l_b) {

    double term1 = -l_a*analytical_overlap_int_1D(X_a, alpha_a, l_a-1, X_b, alpha_b, l_b);
    double term2 = 2*alpha_a*analytical_overlap_int_1D(X_a, alpha_a, l_a+1, X_b, alpha_b, l_b);

    return term1+term2;
};

double S_deriv_one_value(AO &AO_1, AO &AO_2, int direction) {

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
    arma::mat final_mat(3*(atoms.size()-1), pow(AOs.size(),2), arma::fill::zeros);

    for (size_t i = 0; i < final_mat.n_rows; i++) { // iterate through directions x, y, z
        int col_index=0;
        for (auto& ao_1 : AOs) {
            // Fix one AO as AO_1, then iterate through rest of AOs as AO_2
            for (auto& ao_2: AOs) {
                final_mat(i,col_index) = S_deriv_one_value(ao_1, ao_2, i);
                col_index += 1;
            }
        }
    }

    return final_mat;
}



double calc_V_nuc_deriv() {


    return 0;
}


arma::mat V_nuc_deriv(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    // Initialize appropriately sized matrix - each row is a dimension, each column is an atom
    arma::mat final_mat(3, atoms.size(), arma::fill::zeros);


    return final_mat;
}
