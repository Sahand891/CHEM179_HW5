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
#include "Position_Gradient.h"

int main() {

    std::string path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/H2.txt";
    std::vector<Atom> H2_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);

    // Now performing the CNDO/2 energy optimization
    iteration_data final_H2_itdata = perform_CNDO2(H2_Atoms);
    std::cout << compute_total_energy(final_H2_itdata) << std::endl;

    iteration_data final_HF_itdata = perform_CNDO2(HF_Atoms);
    std::cout << compute_total_energy(final_HF_itdata) << std::endl;

    iteration_data final_HO_itdata = perform_CNDO2(HO_Atoms);
    std::cout << compute_total_energy(final_HO_itdata) << std::endl;



    // For H2
    arma::mat H2_S_deriv = S_deriv(final_H2_itdata);
    //H2_S_deriv.print();

    // For HF
    arma::mat HF_S_deriv = S_deriv(final_HF_itdata);
    //HF_S_deriv.print();


    // For HO
    arma::mat HO_S_deriv = S_deriv(final_HO_itdata);
    HO_S_deriv.print();

    return 0;
}