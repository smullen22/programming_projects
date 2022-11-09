#include "header.h"

#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
typedef Eigen::Tensor<double, 4> Tensor4D;
typedef Eigen::Tensor<double, 3> Tensor3D;

 
int main() {
 
CCSD ccsd(7,5);

ccsd.Hcore_form("/home/smullen22/IntroProjects/ProgrammingProjects/Project#05/input/h2o/STO-3G/t.dat", "/home/smullen22/IntroProjects/ProgrammingProjects/Project#05/input/h2o/STO-3G/v.dat");
ccsd.coeff_initial("/home/smullen22/IntroProjects/ProgrammingProjects/Project#03/MO_Project3_data.txt"); //my data from proj3 but off by less than .0000001, probably due to eigen solver used?

ccsd.AO_to_MO();
ccsd.MO_to_spin();

ccsd.spin_Fock();
ccsd.spin_energies();

ccsd.Emp2();

double Ecc_old;
Ecc_old = ccsd.emp2;
ccsd.t1_initial();
ccsd.t2_initial();

ccsd.Tau_tilda(ccsd.t1_in,ccsd.t2_in);
ccsd.Tau(ccsd.t2_in, ccsd.t1_in);

ccsd.F_ae(ccsd.t1_in, ccsd.tau_tilda);
ccsd.F_mi(ccsd.t1_in, ccsd.tau_tilda);
ccsd.F_me(ccsd.t1_in);
ccsd.W_mnij(ccsd.t1_in, ccsd.tau);
ccsd.W_abef(ccsd.t1_in,ccsd.tau);
ccsd.W_mbej(ccsd.t1_in, ccsd.t2_in);

Matrix T1i(ccsd.nsm_occ,ccsd.nsm_vir);
Tensor4D T2i(ccsd.nsm_occ, ccsd.nsm_occ, ccsd.nsm_vir, ccsd.nsm_vir);
Matrix T1i1(ccsd.nsm_occ,ccsd.nsm_vir);
Tensor4D T2i1(ccsd.nsm_occ, ccsd.nsm_occ, ccsd.nsm_vir, ccsd.nsm_vir);

T1i = ccsd.update_T1(ccsd.t1_in, ccsd.t2_in);
T2i = ccsd.update_T2(ccsd.t1_in, ccsd.t2_in);

double Ecc_new;
Ecc_new = ccsd.Ecc(ccsd.updated_T1, ccsd.updated_T2);
int iteration = 1;
cout << "iteration: " << iteration << '\t' << Ecc_new << endl;
iteration =2;
//Tensor3D error(iteration, ccsd.nsm_occ*ccsd.nsm_vir+1,ccsd.nsm_occ*ccsd.nsm_vir);

    ccsd.Tau_tilda(ccsd.updated_T1,ccsd.updated_T2);
    ccsd.Tau(ccsd.updated_T2, ccsd.updated_T1);
    ccsd.F_ae(ccsd.updated_T1, ccsd.tau_tilda);
    ccsd.F_mi(ccsd.updated_T1,ccsd.tau_tilda);
    ccsd.F_me(ccsd.updated_T1);
    ccsd.W_mnij(ccsd.updated_T1, ccsd.tau);
    ccsd.W_abef(ccsd.updated_T1, ccsd.tau);
    ccsd.W_mbej(ccsd.updated_T1, ccsd.updated_T2);
    T1i = ccsd.update_T1(ccsd.updated_T1, ccsd.updated_T2);
    T2i = ccsd.update_T2(ccsd.updated_T1, ccsd.updated_T2);
    Ecc_new = ccsd.Ecc(ccsd.updated_T1, ccsd.updated_T2);
    cout << "iteration: " << iteration << '\t' << Ecc_new << endl;
    

while (abs(Ecc_new-Ecc_old) >= 1E-12) {
     iteration += 1;
    T1i1.setZero();
    T2i1.setZero();
    Ecc_old=Ecc_new;
    ccsd.Tau_tilda(ccsd.updated_T1,ccsd.updated_T2);
    ccsd.Tau(ccsd.updated_T2, ccsd.updated_T1);
    ccsd.F_ae(ccsd.updated_T1, ccsd.tau_tilda);
    ccsd.F_mi(ccsd.updated_T1,ccsd.tau_tilda);
    ccsd.F_me(ccsd.updated_T1);
    ccsd.W_mnij(ccsd.updated_T1, ccsd.tau);
    ccsd.W_abef(ccsd.updated_T1, ccsd.tau);
    ccsd.W_mbej(ccsd.updated_T1, ccsd.updated_T2);
    T1i1 = ccsd.update_T1(ccsd.updated_T1, ccsd.updated_T2);
    T2i1 = ccsd.update_T2(ccsd.updated_T1, ccsd.updated_T2);
    ccsd.Error(iteration, T1i, T1i1, T2i, T2i1);
    //cout << ccsd.B << endl;
    ccsd.Sol_Diis(iteration);
    Ecc_new = ccsd.DisEcc(iteration);;
    cout << "iteration: " << iteration << '\t' << Ecc_new << endl;
    T1i = T1i1;
    T2i = T2i1;
}

// Etot = Ecc_new + Emp2 + Escf(proj 3)


    return 0;
}