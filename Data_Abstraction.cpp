//
// Created by Sahand Adibnia on 3/30/24.
//

#include "Data_Abstraction.h"
#include "Utils.h"

std::vector<AO> atoms_to_AOs(const std::vector<Atom> &Atoms) {

    std::vector<AO> AOs;

    for (auto &atom: Atoms) {
        int atomic_num = atom.A;
        for (auto &cGTO: atom.cGTOs) {
            // just need to take angular momentum of ONE of the primitive Gaussians in a cGTO
            int L = cGTO.primGs[0].L; // will tell me if I have an s or a p orbital
            AO AO_to_append = {atomic_num, L, atom, cGTO};
            AOs.push_back(AO_to_append);
        }
    }

    return AOs;

}


// Counting number of electrons in a molecule
int count_total_electrons(const std::vector<Atom> &Atoms) {

    // Here, n is the total number of electrons provided by atomic orbitals of the atoms
    int n = 0;
    for (auto &atom: Atoms) {
        if (atom.A == 1) { // H atom
            n += 1;
        } else { // C, N, O, or F atom
            n += atom.Z_;
        }
    }
    return n;
}

int count_alpha_electrons(const std::vector<Atom> &Atoms) {
    int n = count_total_electrons(Atoms);
    return ceil(n/2.0);
}

int count_beta_electrons(const std::vector<Atom> &Atoms) {
    int n = count_total_electrons(Atoms);
    return floor(n/2.0);
}


// Overlap integral between two primitive Gaussians, based on HW2 implementation
double analytical_overlap_int_1D(double X_a, double alpha_a, int l_a, double X_b, double alpha_b, double l_b) {

    double coef = exp(-alpha_a*alpha_b*pow((X_a-X_b), 2)/(alpha_a+alpha_b)) * sqrt(M_PI / (alpha_a+alpha_b));

    double sum = 0;
    for (int i=0; i <= l_a; i++) {
        for (int j=0; j <= l_b; j++) {
            if ((i+j) % 2 == 0) { // only care about terms with i+j is even
                double term1 = nChooseR(l_a,i);
                double term2 = nChooseR(l_b,j);
                double Xp = (alpha_a*X_a + alpha_b*X_b) / (alpha_b + alpha_a);
                double term3 = double_factorial((i+j-1)) * pow((Xp - X_a),l_a-i) * pow((Xp-X_b), l_b-j) / pow((2*(alpha_a+alpha_b)),(i+j)/2);
                sum += term1*term2*term3;
            } else {
                sum += 0;
            }
        }
    }

    return coef*sum;

}

double S_ab(prim_Gaussian primG1, prim_Gaussian primG2) {

    // Overlap integral between two 3D primitive Gaussians

    // All integrals in x,y,z directions are all calculated the same, then multiplied to get total overlap integral!

    // x direction
    int l_a = primG1.l;
    int l_b = primG2.l;
    double x_int = analytical_overlap_int_1D(primG1.X, primG1.alpha, l_a, primG2.X, primG2.alpha,l_b);

    // y-direction
    int m_a = primG1.m;
    int m_b = primG2.m;
    double y_int = analytical_overlap_int_1D(primG1.Y, primG1.alpha, m_a, primG2.Y, primG2.alpha,m_b);

    // z-direction
    int n_a = primG1.n;
    int n_b = primG2.n;
    double z_int = analytical_overlap_int_1D(primG1.Z, primG1.alpha, n_a, primG2.Z, primG2.alpha,n_b);

    return x_int*y_int*z_int;

}

// For calculating overlap integral between two Gaussians
double normalizer(prim_Gaussian primG) {
    double overlap = S_ab(primG, primG);
    return sqrt(1/overlap);
}


Atom compile_2nd_row_atom(double X, double Y, double Z, STO3G sto3g, int A, double IA_s, double IA_p, double beta) {

    // 2s - 3 primitive Gaussians
    prim_Gaussian primG1_2s = {X, Y, Z, sto3g.exps[0]};
    prim_Gaussian primG2_2s = {X, Y, Z, sto3g.exps[1]};
    prim_Gaussian primG3_2s = {X, Y, Z, sto3g.exps[2]};
    std::vector<prim_Gaussian> primGs_2s = {primG1_2s, primG2_2s, primG3_2s};
    // 2px - 3 primitive Gaussians
    prim_Gaussian primG1_2px = {X, Y, Z, sto3g.exps[0], 1,0,0};
    prim_Gaussian primG2_2px = {X, Y, Z, sto3g.exps[1], 1,0,0};
    prim_Gaussian primG3_2px = {X, Y, Z, sto3g.exps[2], 1, 0,0};
    std::vector<prim_Gaussian> primGs_2px = {primG1_2px, primG2_2px, primG3_2px};
    // 2py - 3 primitive Gaussians
    prim_Gaussian primG1_2py = {X, Y, Z, sto3g.exps[0], 0,1,0};
    prim_Gaussian primG2_2py = {X, Y, Z, sto3g.exps[1], 0,1,0};
    prim_Gaussian primG3_2py = {X, Y, Z, sto3g.exps[2], 0, 1,0};
    std::vector<prim_Gaussian> primGs_2py = {primG1_2py, primG2_2py, primG3_2py};
    // 2pz - 3 primitive Gaussians
    prim_Gaussian primG1_2pz = {X, Y, Z, sto3g.exps[0], 0,0,1};
    prim_Gaussian primG2_2pz = {X, Y, Z, sto3g.exps[1], 0,0,1};
    prim_Gaussian primG3_2pz = {X, Y, Z, sto3g.exps[2], 0, 0,1};
    std::vector<prim_Gaussian> primGs_2pz = {primG1_2pz, primG2_2pz, primG3_2pz};

    // Vector of contraction coefficients - 2s
    std::vector<double> contr_coefs_2s = sto3g.con_coef_s;
    // Vector of contraction coefficients - 2p (same for all 2p orbitals)
    std::vector<double> contr_coefs_2p = sto3g.con_coef_p;

    // Vector of normalization constants - 2s
    std::vector<double> norm_constants_2s = {normalizer(primG1_2s), normalizer(primG2_2s), normalizer(primG3_2s)};
    // Vector of normalization constants - 2p (same for all 2p orbitals)
    std::vector<double> norm_constants_2p = {normalizer(primG1_2px), normalizer(primG2_2px), normalizer(primG3_2px)};

    // Make the vector of contracted Gaussians (there will be 4)
    // 2s contracted Gaussian
    contr_Gaussian cGTO_2s = {primGs_2s, contr_coefs_2s, norm_constants_2s};
    // 2px contracted Gaussian
    contr_Gaussian cGTO_2px = {primGs_2px, contr_coefs_2p, norm_constants_2p};
    // 2px contracted Gaussian
    contr_Gaussian cGTO_2py = {primGs_2py, contr_coefs_2p, norm_constants_2p};
    // 2px contracted Gaussian
    contr_Gaussian cGTO_2pz = {primGs_2pz, contr_coefs_2p, norm_constants_2p};
    std::vector<contr_Gaussian> cGTOs = {cGTO_2s, cGTO_2px, cGTO_2py,cGTO_2pz};

    Atom final_atom = {X,Y,Z,cGTOs,A,IA_s,IA_p,beta};

    return final_atom;

}

Atom construct_H_Atom(double X, double Y, double Z) {

    STO3G H_STO3G = {{3.42525091, 0.62391373, 0.16885540},{0.15432897, 0.53532814, 0.44463454}};
    STO3G sto3g = H_STO3G;
    // 3 primitive Gaussian per orbital, only 1 orbital in H atom
    prim_Gaussian primG1 = {X, Y, Z, sto3g.exps[0]};
    prim_Gaussian primG2 = {X, Y, Z, sto3g.exps[1]};
    prim_Gaussian primG3 = {X, Y, Z, sto3g.exps[2]};
    std::vector<prim_Gaussian> primGs = {primG1, primG2, primG3};

    // Vector of contraction coefficients
    std::vector<double> contr_coefs = sto3g.con_coef_s;

    // Vector of normalization constants
    std::vector<double> norm_constant = {normalizer(primG1), normalizer(primG2), normalizer(primG3)};

    // Make the contracted Gaussian, put it into a 1-element vector (for data type consistency)
    contr_Gaussian cGTO = {primGs, contr_coefs, norm_constant};
    std::vector<contr_Gaussian> cGTOs = {cGTO};

    Atom final_atom = {X,Y,Z,cGTOs,1,7.176,0,-9};

    return final_atom;

}

Atom construct_Atom(double X, double Y, double Z, std::string atom_symbol) {

    if (atom_symbol == "H") {
        return construct_H_Atom(X,Y,Z);
    }

    // just default to C
    STO3G C_STO3G = {{2.94124940, 0.68348310, 0.22228990}, {-0.09996723, 0.39951283, 0.70011547}, {0.15591627, 0.60768372, 0.39195739}};
    STO3G sto3g = C_STO3G;
    int A = 6;
    double IA_s = 14.051;
    double IA_p = 5.572;
    double beta = -21;

    if (atom_symbol == "C") {
        sto3g = C_STO3G;
        A = 6;
        IA_s = 14.051;
        IA_p = 5.572;
        beta = -21;
    } else if (atom_symbol == "N") {
        STO3G N_STO3G = {{0.3780455879E+01, 0.8784966449E+00, 0.2857143744E+00},{-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00},{0.1559162750E+00, 0.6076837186E+00, 0.3919573931E+00}};
        sto3g = N_STO3G;
        A = 7;
        IA_s = 19.316;
        IA_p = 7.275;
        beta = -25;
    } else if (atom_symbol == "O") {
        STO3G O_STO3G = {{5.033151319, 1.169596125, 0.3803889600},{-0.09996722919, 0.3995128261, 0.7001154689},{0.1559162750, 0.6076837186, 0.3919573931}};
        sto3g = O_STO3G;
        A = 8;
        IA_s = 25.390;
        IA_p = 9.111;
        beta = -31;
    } else if (atom_symbol == "F") {
        STO3G F_STO3G = {{0.6464803249E+01, 0.1502281245E+01, 0.4885884864E+00},{-0.9996722919E-01, 0.3995128261E+00, 0.7001154689E+00},{0.15591627, 0.60768372, 0.39195739}};
        sto3g = F_STO3G;
        A = 9;
        IA_s = 32.272;
        IA_p = 11.080;
        beta = -39;
    } else {
        throw std::invalid_argument("There is some problem with the input Atom format.");
    }

    return compile_2nd_row_atom(X, Y, Z, sto3g, A, IA_s, IA_p, beta);

}



// Reads in atoms from file, returns a vector of atoms
std::vector<Atom> read_atoms_from_file(const std::string &path) {

    std::ifstream inputFile(path);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
    }

    std::vector<Atom> Atoms;

    // Read values from the .txt file, skipping the first line (which contains information about the number of atoms at the charge)
    std::string line;
    getline(inputFile, line); // get the number of atoms
    size_t num_atoms = stoull(line);
    int A;
    double X, Y, Z;
    while (inputFile >> A >> X >> Y >> Z) {
        //std::cout << A << X << Y << Z << std::endl;
        if (A==1) {
            //std::cout << abs(A-1) << std::endl;
            Atoms.push_back(construct_H_Atom(X,Y,Z));
        } else if (A==6) {
            Atoms.push_back(construct_Atom(X,Y,Z,"C"));
        } else if (A==7) {
            Atoms.push_back(construct_Atom(X,Y,Z,"N"));
        } else if (A==8) {
            Atoms.push_back(construct_Atom(X,Y,Z,"O"));
        } else if (A==9) {
            Atoms.push_back(construct_Atom(X,Y,Z,"F"));
        }
    }

    inputFile.close();

    return Atoms;

}




