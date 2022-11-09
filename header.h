#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <cassert>
#include "/home/smullen22/IntroProjects/ProgrammingProjects/Project#04/eigen/Eigen/Dense"
#include "/home/smullen22/IntroProjects/ProgrammingProjects/Project#04/eigen/Eigen/Eigenvalues"
#include "/home/smullen22/IntroProjects/ProgrammingProjects/Project#04/eigen/Eigen/Core"
#include "/home/smullen22/IntroProjects/ProgrammingProjects/Project#04/eigen/Eigen/LU"
#include "/home/smullen22/IntroProjects/ProgrammingProjects/Project#04/eigen/unsupported/Eigen/CXX11/Tensor"


#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
typedef Eigen::Tensor<double, 4> Tensor4D;
typedef Eigen::Tensor<double, 3> Tensor3D;

using namespace std;
 
class CCSD
{ 
    public: 
    int nao;
    int nocc=5;
    int nvir=2;
    int nsm=14;
    int nsm_occ=10;
    int nsm_vir=4;
    double *TEI;
    double *TEI_MO;
    Matrix Hcore;
    Matrix coeff_in;
    Tensor4D TEI_spin;
    Matrix Fock_spin;
    Matrix energies;
    double emp2;
    Matrix t1_in;
    Tensor4D t2_in;
    Tensor4D tau_tilda;
    Tensor4D tau;
    Tensor4D Wmnij;
    Tensor4D Wabef;
    Tensor4D Wmbej;
    Matrix Fae;
    Matrix Fmi;
    Matrix Fme;
    Matrix updated_T1;
    Tensor4D updated_T2;
    Matrix err;
    int numv;
    Tensor3D disTen;
    Tensor3D Err;
    double diisEcc;
    Matrix B;
     
    Matrix Hcore_form(const char *file1, const char *file2);
    Matrix coeff_initial(const char *file);
    double AO_to_MO();
    Tensor4D MO_to_spin();
    Matrix spin_Fock();
    Matrix spin_energies();
    double Emp2();
    Matrix t1_initial();
    Tensor4D t2_initial();
    Tensor4D Tau_tilda(Matrix t1, Tensor4D t2);
    Matrix F_ae(Matrix t1, Tensor4D tau_tilda);
    Matrix F_mi(Matrix t1, Tensor4D tau_tilda);
    Matrix F_me(Matrix t1);
    Tensor4D Tau(Tensor4D t2, Matrix t1);
    Tensor4D W_mnij(Matrix t1, Tensor4D tau);
    Tensor4D W_abef(Matrix t1, Tensor4D tau);
    Tensor4D W_mbej(Matrix t1, Tensor4D t2);
    Matrix update_T1(Matrix t1_old, Tensor4D t2);
    Tensor4D update_T2(Matrix t1_old, Tensor4D t2_old);
    double Ecc(Matrix updated_t1, Tensor4D t2_updated);
    Tensor3D Error(int it, Matrix T1i, Matrix T1i1, Tensor4D T2i, Tensor4D T2i1);
    Tensor3D Sol_Diis(int it);
    double DisEcc(int it);

    CCSD(int nao, int nocc);
    ~CCSD();
};
