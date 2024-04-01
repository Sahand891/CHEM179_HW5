//
// Created by Sahand Adibnia on 3/31/24.
//

#include "Finite_Differences.h"

// Code does not work at the moment, give a systematically incorrect gradient



// We take the derivative with respect to movement of a certain atom
double central_differences(int deriv_atom_index, const std::vector<Atom> &atoms, double h, int direction) {

    // Adding h to specified direction for specified atom
    std::vector<Atom> atoms_pos_h = atoms;
    if (direction == 0) {
        atoms_pos_h[deriv_atom_index].X += h;
    } else if (direction == 1) {
        atoms_pos_h[deriv_atom_index].Y += h;
    } else {
        atoms_pos_h[deriv_atom_index].Z += h;
    }

    // Subtracting h from specified direction for specified atom
    std::vector<Atom> atoms_neg_h = atoms;
    if (direction == 0) {
        atoms_neg_h[deriv_atom_index].X -= h;
    } else if (direction == 1) {
        atoms_neg_h[deriv_atom_index].Y -= h;
    } else {
        atoms_neg_h[deriv_atom_index].Z -= h;
    }

    iteration_data conv_itdata_pos_h = perform_CNDO2(atoms_pos_h);
    iteration_data conv_itdata_neg_h = perform_CNDO2(atoms_neg_h);

    double term1 = compute_total_energy(conv_itdata_pos_h);
    double term2 = compute_total_energy(conv_itdata_neg_h);
    return -(term1 - term2) / (2*h);

}


// Gradient vector for a specific atom in the molecule
arma::vec central_diff_grad(int deriv_atom_index, const std::vector<Atom> &atoms, double h) {

    // Initialize appropriately sized vector
    arma::vec final_vec(3, arma::fill::zeros);

    // Iterate over direction
    for (int direc=0; direc < 3; direc++) {
        final_vec[direc] = central_differences(deriv_atom_index, atoms, h, direc);
    }

    return final_vec;
}

// Matrix of gradient vectors for all atoms in molecule
arma::mat central_diff_mat(const std::vector<Atom> &atoms, double h) {

    // Initialize appropriately sized matrix
    arma::mat final_mat(3, atoms.size(), arma::fill::zeros);

    // Iterate over each atom and direction
    for (int atom_index=0; atom_index < atoms.size(); atom_index++) {
        for (int direc=0; direc < 3; direc++) {
            final_mat(direc, atom_index) = central_differences(atom_index, atoms, h, direc);
        }
    }

    return final_mat;

}
