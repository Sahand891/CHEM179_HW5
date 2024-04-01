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
#include "Finite_Differences.h"

int main() {

    // Reading in the atoms and performing CNDO/2 calculations to get the energy at the specified atomic distances
    std::string path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/H2.txt";
    std::vector<Atom> H2_Atoms = read_atoms_from_file(path);
    iteration_data final_H2_itdata = perform_CNDO2(H2_Atoms);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/HF.txt";
    std::vector<Atom> HF_Atoms = read_atoms_from_file(path);
    iteration_data final_HF_itdata = perform_CNDO2(HF_Atoms);

    path = "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/Xiao/HW5/sample_input/HO.txt";
    std::vector<Atom> HO_Atoms = read_atoms_from_file(path);
    iteration_data final_HO_itdata = perform_CNDO2(HO_Atoms);

    // Computing the gradient (as a matrix, each column being gradient with respect to movement of one atom) and outputting into a nice txt file
    write_gradient_to_file(final_H2_itdata, "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/my_outputs/H2.txt");
    write_gradient_to_file(final_HF_itdata, "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/my_outputs/HF.txt");
    write_gradient_to_file(final_HO_itdata, "/Users/sahandadibnia/homeworks/CHEM 179/Homework/HW5/my_outputs/HO.txt");

    return 0;
}