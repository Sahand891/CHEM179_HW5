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

    // Initialize appropriately sized matrix - each row is a dimension, and each column is an AO-AO pair (e.g. 1sA-1sA 1sA-1sB 1sB-1sA 1sB-1sB for H2)
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

    // Initialize appropriately sized matrix - each row is a dimension, each column an atom-atom pair (e.g. AA, AB, BA, BB, for H2 case)
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

// Derivative with respect to atom A
double V_nuc_deriv_one_element(Atom &A, std::vector<Atom> &atoms, int direction) {

    double sum=0;
    for (auto& atom : atoms) {

        double term1 = A.Z_*atom.Z_;
        double term2 = (A.X - atom.X);
        if (direction == 1) {
            term2 = (A.Y - atom.Y);
        } else if (direction == 2) {
            term2 = (A.Z - atom.Z);
        }

        // if atom=A, ignore!
        if (arma::norm((A.pos_vec - atom.pos_vec)) == 0) {
            sum += 0;
        } else {
            sum += -term1*term2 / pow(arma::norm((A.pos_vec - atom.pos_vec)), 3);
        }
    }

    return sum * 27.2114079527; // Hartree to eV conversion
}


arma::mat V_nuc_deriv(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    // Initialize appropriately sized matrix - each row is a dimension, each column is an atom
    arma::mat final_mat(3, atoms.size(), arma::fill::zeros);

    for (size_t i = 0; i < 3; i++) {
        int col_index=0;
        for (auto& atom : atoms) {
            final_mat(i,col_index) = V_nuc_deriv_one_element(atom, atoms, i);
            col_index += 1;
        }
    }


    return final_mat;
}




// Creating energy gradient


double x_uv(int u, int v, const iteration_data &conv_it_data) {
    // u, v indicate AO indices

    double beta_A = conv_it_data.AOs[u].origin_atom.beta;
    double beta_B = conv_it_data.AOs[v].origin_atom.beta;
    double P_uv_tot = conv_it_data.P_alpha_new(u,v) + conv_it_data.P_beta_new(u,v);

    return (beta_A+beta_B)*P_uv_tot;


}

arma::mat x_mat(const iteration_data &conv_it_data) {

    std::vector<AO> AOs = conv_it_data.AOs;

    arma::mat final_mat(AOs.size(),AOs.size(),arma::fill::zeros);

    for (int u=0; u < AOs.size(); u++) {
        for (int v=0; v < AOs.size(); v++) {
            final_mat(u,v) = x_uv(u,v,conv_it_data);
        }
    }

    return final_mat;
}


double y_ab(int a, int b, const iteration_data &conv_it_data) {
    // a, b indicate atom indices

    std::vector<Atom> atoms = conv_it_data.atoms;

    Atom atomA = atoms[a];
    Atom atomB = atoms[b];

    double P_AA_tot = conv_it_data.P_total_new[a]; // P_total is an armadillo vector
    double P_BB_tot = conv_it_data.P_total_new[b];

    // Calculating the tricky summation term
    arma::mat P_alpha = conv_it_data.P_alpha_new;
    arma::mat P_beta = conv_it_data.P_beta_new;

    // Determining the AO indices associated with each atom
    std::vector<int> a_AO_indices = obtain_AO_indices_of_atom(a, atoms);
    std::vector<int> b_AO_indices = obtain_AO_indices_of_atom(b, atoms);

    double summation_term = 0;
    for (auto& u : a_AO_indices) {
        for (auto& v : b_AO_indices) {
            summation_term += P_alpha(u,v)*P_alpha(u,v) + P_beta(u,v)*P_beta(u,v);
        }
    }

    return P_AA_tot*P_BB_tot - atomB.Z_*P_AA_tot - atomA.Z_*P_BB_tot - summation_term;

}

arma::mat y_mat(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    arma::mat final_mat(atoms.size(),atoms.size(),arma::fill::zeros);

    for (int a=0; a < atoms.size(); a++) {
        for (int b=0; b < atoms.size(); b++) {
            final_mat(a,b) = y_ab(a,b,conv_it_data);
        }
    }

    return final_mat;
}


// We take the gradient with respect to a certain atom's position! That atom is the atom at the deriv_atom_index
double electron_gradient_x(const iteration_data &conv_it_data, int deriv_atom_index, int direction) {

    arma::mat x = x_mat(conv_it_data);
    arma::mat y = y_mat(conv_it_data);

    //x.print();
    //y.print();

    std::vector<Atom> atoms = conv_it_data.atoms;
    std::vector<AO> AOs = conv_it_data.AOs;

    Atom deriv_atom = atoms[deriv_atom_index];
    std::vector<int> AOs_indices_deriv_atom = obtain_AO_indices_of_atom(deriv_atom_index, atoms);

    double overlap_sum = 0;
    // Had to implement these nested for loops in special way to avoid double counting (adding both 01 and 10, which cancel!)
    for (auto& u : AOs_indices_deriv_atom) {
        for (int v=0; v < AOs.size(); v++) {
            // Check if the AOs are the same—if so, don't calculate anything
            //std::cout << u << v << std::endl;
            if (u == v) {
                overlap_sum += 0;
            } else {
                overlap_sum += x(u,v)*S_deriv_one_element(AOs[u], AOs[v], direction);
                //std::cout << x(u,v)*S_deriv_one_element(AOs[u], AOs[v], direction) << std::endl;
            }
            //std::cout << overlap_sum << std::endl;
        }
    }


    double gamma_sum=0;
    for (int b=0; b < atoms.size(); b++) {
        // Check if the atom we're on is the same as the atom we're deriving with respect to—if so, don't calculate anything
        if (deriv_atom_index == b) {
            gamma_sum += 0;
        } else {
            gamma_sum += y(deriv_atom_index,b)*gamma_deriv_one_element(deriv_atom, atoms[b], direction);
        }
    }


    return overlap_sum+gamma_sum;

}

arma::mat electron_gradient_mat(const iteration_data &conv_it_data) {

    std::vector<Atom> atoms = conv_it_data.atoms;

    // Initialize appropriately sized matrix - rows are x,y,z, 1 column per atom (each column gives gradient of atom position)
    arma::mat final_mat(3, atoms.size(), arma::fill::zeros);

    for (int j=0; j < final_mat.n_cols; j++) {
        for (int i=0; i < final_mat.n_rows; i++) { // go through directions first
            final_mat(i,j) = electron_gradient_x(conv_it_data, j, i);
        }
    }

    return final_mat;

}


arma::mat total_gradient_mat(const iteration_data &conv_it_data) {
    return electron_gradient_mat(conv_it_data) + V_nuc_deriv(conv_it_data);
}