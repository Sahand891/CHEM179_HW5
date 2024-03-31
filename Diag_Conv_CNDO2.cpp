//
// Created by Sahand Adibnia on 3/30/24.
//

#include "Diag_Conv_CNDO2.h"

arma::mat find_MO_coefs(arma::mat F) {
    arma::vec energies;
    arma::mat C;
    arma::eig_sym(energies, C, F); // can use eig_sym because F is guranteed to be symmetric

    return C;
}


// Calculate energy eigenvalues of a converged Fock matrix
arma::vec get_energy_eigs(arma::mat F) {
    arma::vec energies;
    arma::mat C;
    arma::eig_sym(energies, C, F); // can use eig_sym because F is guaranteed to be symmetric

    return energies;
}


// Calculate distance between two atoms
double atom_distance(const Atom &A, const Atom &B) {
    return sqrt(pow((A.X - B.X),2) + pow((A.Y - B.Y),2) + pow((A.Z - B.Z),2));
}

// Calculate nuclear repulsion energy for set of atoms
double nuc_repl_energy(const std::vector<Atom> &atoms) {
    double sum=0;
    for (auto& A : atoms) {
        for (auto& B : atoms) {
            if (A.X == B.X and A.Y == B.Y and A.Z == B.Z) { // check if they're the same atom
                sum += 0;
            } else {
                double R_AB = atom_distance(A,B);
                sum += A.Z_*B.Z_ / R_AB;
            }
        }
    }
    return sum * 27.211 / 2; // multiply by eV/a.u. ratio, divide by 2 to eliminate double counting
}


// Calculate electron energy from converged density and Fock matrices
double electron_energy(const arma::mat &F_alpha, const arma::mat &F_beta, const arma::mat &P_alpha, const arma::mat &P_beta, const arma::mat &H_core) {

    int matrix_length = sqrt(F_alpha.size());

    // Iterating through indices of the matrices
    double alpha_term=0, beta_term=0;
    for (int u=0; u < matrix_length; u++) {
        for (int v=0; v < matrix_length; v++) {
            alpha_term += P_alpha(u,v)*(H_core(u,v)+F_alpha(u,v));
            beta_term += P_beta(u,v)*(H_core(u,v)+F_beta(u,v));
        }
    }

    return (alpha_term+beta_term) / 2;
}


// Calculate total energy of a converged system (sum of nuclear repulsion and electron energies)
double compute_total_energy(const iteration_data &conv_it_data) {

    double nucl_rep_energy = nuc_repl_energy(conv_it_data.atoms);
    double electron_e = electron_energy(conv_it_data.fock_alpha, conv_it_data.fock_beta, conv_it_data.P_alpha_new, conv_it_data.P_beta_new, conv_it_data.H_core);

    return nucl_rep_energy+electron_e;
}


// This function starts CNDO2 density matrix optimization with an initial guess of the core Hamiltonian matrix, and also prints out the 0th iteration
iteration_data start_CNDO2(const std::vector<Atom> &atoms, bool do_print, bool auto_electrons, int p_input, int q_input) {

    arma::mat gamma = gamma_matrix(atoms);

    arma::mat overlap = S(atoms);

    arma::mat H_core = h(atoms);

    // Compute appropriate matrices to input into iteration_data object

    int p = p_input;
    int q = q_input;

    if (auto_electrons) { // automatically calculate # of alpha and beta electrons if so specified
        p = count_alpha_electrons(atoms);
        q = count_beta_electrons(atoms);
    }

    // C_alpha = C_beta = C for the zeroth iteration since we're using core Hamiltonian to start

    // 0th iteration: Diagonalizing H_core to get the MO coefficients
    arma::mat C = find_MO_coefs(H_core);

    // Calculating "new" density matrices for the 0th iteration (i.e., starting density matrices)
    arma::mat C_occ_a = C_occ_alpha(C, p);
    arma::mat C_occ_b = C_occ_alpha(C, q);
    arma::mat P_a = P_alpha(C_occ_a);
    arma::mat P_b = P_alpha(C_occ_b);

    // For the first iteration, old matrix just set to be a 1x1 matrix of 0's, and specifically set converged bool to be false (SUPER IMPORTANT!)
    // We won't even use P_alpha_old and P_beta_old in the 0th iteration anyway

    iteration_data it_data = {atoms,
                              p,
                              q,
                              arma::mat(1,1,arma::fill::zeros),
                              arma::mat(1,1,arma::fill::zeros),
                              H_core,
                              H_core,
                              C,
                              C,
                              P_a,
                              P_b,
                              0,
                              false};


    if (do_print) {
        // Now print out everything nicely so we can see it!!
        std::cout << "gamma" << std::endl;
        gamma.print();
        std::cout << "Overlap" << std::endl;
        overlap.print();
        std::cout << "p = " << p << " q = " << q << std::endl;
        std::cout << "H_core" << std::endl;
        H_core.print();
        // Printing out zeroth iteration
        std::cout << "Iteration: 0" << std::endl;
        std::cout << "Fa" << std::endl;
        H_core.print();
        std::cout << "Fb" << std::endl;
        H_core.print();
        std::cout << "after solving the eigenvalue equation: " << std::endl;
        std::cout << "Ca" << std::endl;
        C.print();
        std::cout << "Cb" << std::endl;
        C.print();
        std::cout << "p = " << it_data.p << " q = " << it_data.q << std::endl;
        std::cout << "Pa_new" << std::endl;
        P_a.print();
        std::cout << "Pb_new" << std::endl;
        P_b.print();
        std::cout << "P_total" << std::endl;
        P_AA(atoms, P_a, P_b).print();
    }

    return it_data;

}


std::vector<iteration_data> converge_CNDO2(const iteration_data &it_data, std::vector<iteration_data> &it_data_vec, bool do_print) {

    // Base case
    if (it_data.converged) {
        return it_data_vec;
    }

    // Extracting constant information for each iteration
    std::vector<Atom> atoms = it_data.atoms;
    int p = it_data.p;
    int q = it_data.q;


    // Old density matrices for this iteration at the new ones from the previous iteration
    arma::mat P_alpha_old = it_data.P_alpha_new;
    arma::mat P_beta_old = it_data.P_beta_new;


    // Make the new fock matrix based on the density matrices from the previous iteration
    arma::mat Fock_alpha = F_alpha(atoms, P_alpha_old, P_beta_old);
    arma::mat Fock_beta = F_alpha(atoms, P_beta_old, P_alpha_old);


    // Diagonalize the fock matrices and extra MO coefficients
    arma::mat C_alpha_new = find_MO_coefs(Fock_alpha);
    arma::mat C_beta_new = find_MO_coefs(Fock_beta);


    // Calculate the new density matrices based on the new MO coefficient matrices
    // Compare new vs old density matrices
    arma::mat C_occ_alpha_new = C_occ_alpha(C_alpha_new, p);
    arma::mat P_alpha_new = P_alpha(C_occ_alpha_new);

    arma::mat C_occ_beta_new = C_occ_alpha(C_beta_new, q);
    arma::mat P_beta_new = P_alpha(C_occ_beta_new);


    // Compile all the appropriate info into a new iteration_data object
    iteration_data new_it_data = {atoms,
                                  p,
                                  q,
                                  P_alpha_old,
                                  P_beta_old,
                                  Fock_alpha,
                                  Fock_beta,
                                  C_alpha_new,
                                  C_beta_new,
                                  P_alpha_new,
                                  P_beta_new,
                                  it_data.iteration_count+1};

    // Adding that new set of iteration data to the original vector of iteration data
    it_data_vec.push_back(new_it_data);

    if (do_print) {
        std::cout << "Iteration: " << it_data.iteration_count << std::endl;
        std::cout << "Fa" << std::endl;
        Fock_alpha.print();
        std::cout << "Fb" << std::endl;
        Fock_beta.print();
        std::cout << "after solving the eigenvalue equation: " << std::endl;
        std::cout << "Ca" << std::endl;
        C_alpha_new.print();
        std::cout << "Cb" << std::endl;
        C_beta_new.print();
        std::cout << "p = " << p << " q = " << q<< std::endl;
        std::cout << "Pa_new" << std::endl;
        P_alpha_new.print();
        std::cout << "Pb_new" << std::endl;
        P_beta_new.print();
        std::cout << "P_total" << std::endl;
        new_it_data.P_total_new.print();
    }

    return converge_CNDO2(new_it_data, it_data_vec, do_print);

}

iteration_data obtain_converged_data(std::vector<iteration_data> &converged_it_data) {
    return converged_it_data.back();
}



void write_CNDO2_to_file(std::string output_location, iteration_data &starting_it_data, std::vector<iteration_data> &it_data_vec) {

    std::ofstream outFile(output_location);

    // Check if the file is successfully opened
    if (!outFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
    }

    // Start with the data from starting_it_data (i.e., the H_core stuff)
    outFile << "Gamma Matrix" << std::endl;
    starting_it_data.gamma.print(outFile);
    outFile << "Overlap Matrix" << std::endl;
    starting_it_data.overlap.print(outFile);
    outFile << "p = " << starting_it_data.p << " q = " << starting_it_data.q << std::endl;
    outFile << "H_core" << std::endl;
    starting_it_data.fock_alpha.print(outFile);

    // Now print out all the iteration data!
    for (auto& it_data : it_data_vec) {
        outFile << "Iteration: " << it_data.iteration_count << std::endl;
        outFile << "Fa" << std::endl;
        it_data.fock_alpha.print(outFile);
        outFile << "Fb" << std::endl;
        it_data.fock_beta.print(outFile);
        outFile << "After solving the eigenvalue equation: " << std::endl;
        outFile << "Ca" << std::endl;
        it_data.C_alpha_new.print(outFile);
        outFile << "Cb" << std::endl;
        it_data.C_beta_new.print(outFile);
        outFile << "p = " << it_data.p << " q = " << it_data.q << std::endl;
        outFile << "Pa_new" << std::endl;
        it_data.P_alpha_new.print(outFile);
        outFile << "Pb_new" << std::endl;
        it_data.P_beta_new.print(outFile);
        outFile << "P_total" << std::endl;
        it_data.P_total_new.print(outFile);
    }

    // Now print out the final, converged data
    iteration_data conv_data = obtain_converged_data(it_data_vec);

    arma::mat Ea = get_energy_eigs(conv_data.fock_alpha);
    arma::mat Eb = get_energy_eigs(conv_data.fock_beta);

    double nuc_repl_e = nuc_repl_energy(conv_data.atoms);
    double electron_e = electron_energy(conv_data.fock_alpha, conv_data.fock_beta, conv_data.P_alpha_new, conv_data.P_beta_new, conv_data.H_core);

    outFile << "\nFinal Converged Molecule:" << std::endl;
    outFile << "Ea (Energy eigenvalue matrix for alpha electrons)" << std::endl;
    Ea.print(outFile);
    outFile << "Eb (Energy eigenvalue matrix for beta electrons)" << std::endl;
    Eb.print(outFile);
    outFile << "Ca (MO coefficient matrix for alpha electrons)" << std::endl;
    conv_data.C_alpha_new.print(outFile);
    outFile << "Cb (MO coefficient matrix for beta electrons)" << std::endl;
    conv_data.C_beta_new.print(outFile);
    outFile << "Nuclear Repulsion Energy: " << nuc_repl_e << " eV" << std::endl;
    outFile << "Electron Energy: " << electron_e << " eV" << std::endl;
    outFile << "This molecule has energy: " << electron_e+nuc_repl_e << " eV" << std::endl;

    outFile.close();

}



iteration_data perform_CNDO2(const std::vector<Atom> &atoms) {

    iteration_data initial_itdata = start_CNDO2(atoms);
    std::vector<iteration_data> start_vec = {initial_itdata};
    std::vector<iteration_data> conv_itdatas = converge_CNDO2(initial_itdata, start_vec);
    iteration_data final_itdata = obtain_converged_data(conv_itdatas);

    return final_itdata;
}
