//
// Created by Sahand Adibnia on 3/30/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>
#include <armadillo>

#ifndef HW5_DATA_ABSTRACTION_H
#define HW5_DATA_ABSTRACTION_H

// All the important data structures first

struct prim_Gaussian {

    /* Every unnormalized primitive Gaussian (3D) has:
     * 1. A set of coordinates: X, Y, Z
     * 2. An exponent (alpha)
     * 3. Dimensional angular momentum (exponents for X, Y, Z coords)
     * 4. Total angular momentum (sum of dimensional ang momentums)
     */

    double X, Y, Z, alpha;
    int l=0, m=0, n=0;
    int L = l+m+n;

};

struct contr_Gaussian {
    std::vector<prim_Gaussian> primGs;
    std::vector<double> contr_coefs;
    std::vector<double> norm_constants;

    arma::vec d_primes = arma::vec(contr_coefs) % arma::vec(norm_constants); // element-wise multiplication
};

struct STO3G {
    std::vector<double> exps;
    std::vector<double> con_coef_s;
    std::vector<double> con_coef_p;
};


struct Atom {
    // General atom data type

    double X,Y,Z;
    std::vector<contr_Gaussian> cGTOs;
    int A; // atomic number

    // Semi-emiprical CNDO/2 parameters
    double IA_s;
    double IA_p;
    double beta;

    // valence atomic number
    int Z_ = std::max(1, A-2);

    // position vector
    arma::vec pos_vec = {X,Y,Z};
};


// Atomic orbital struct
struct AO {
    int A, L; // each AO contains information of which atom is came from!! Best to iterate through AOs!!!
    Atom origin_atom;
    contr_Gaussian cGTO;
};



// The relevant functions to construct and mess around with the data structures
std::vector<Atom> read_atoms_from_file(const std::string& path);
Atom construct_Atom(double X, double Y, double Z, std::string atom_symbol);
double analytical_overlap_int_1D(double X_a, double alpha_a, int l_a, double X_b, double alpha_b, double l_b);
double S_ab(prim_Gaussian primG1, prim_Gaussian primG2);

std::vector<AO> atoms_to_AOs(const std::vector<Atom> &Atoms);
int count_total_electrons(const std::vector<Atom> &Atoms);
int count_alpha_electrons(const std::vector<Atom> &Atoms);
int count_beta_electrons(const std::vector<Atom> &Atoms);



#endif //HW5_DATA_ABSTRACTION_H
