#include "header.h"

#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) : (((j)*((j)+1)/2)+(i)))

using namespace std;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
typedef Eigen::Tensor<double, 4> Tensor4D;
typedef Eigen::Tensor<double, 3> Tensor3D;

// reads and stores kinetic energy and nuclear attraction integrals
// calculates Hcore
Matrix CCSD::Hcore_form(const char *file1, const char *file2) {
  ifstream ke(file1);
  assert(ke.good());
  Matrix T(nao,nao);
  double b;
  int j,k; 
  for (int i=0; i<nao*nao; i++) {
    ke >> j >> k >> T(j-1, k-1);
    b = T(j-1, k-1);
    T(k-1, j-1) = b;
  } 

  ifstream nuc_at(file2);
  assert(nuc_at.good());
  Matrix V(nao,nao);
  double a;
  for (int i=0; i<49; i++) {
    nuc_at >> j >> k >> V(j-1, k-1);
    a = V(j-1, k-1);
    V(k-1, j-1) = a;
  }

  Hcore.resize(nao, nao); 
  for (int i=0; i<nao; i++) {
    for (int j=0; j<nao; j++) {
      Hcore(i,j) = T(i,j) + V(i,j);
    }
  }
  return Hcore;
}

//create Matrix of MO coefficients using last coefficient matrix from proj3
Matrix CCSD::coeff_initial(const char *file) {
  int num1,num2;
  coeff_in.resize(nao,nao);
  ifstream MO_data(file);
  assert(MO_data.good());
  for (int i=0; i<nao; i++) {
    for (int j=0; j<nao;j++) {
      MO_data >> coeff_in(i,j);
    }
  }
  return coeff_in;
}

//Transform the Two-Electron Integrals to the MO Basis (proj4)
double CCSD::AO_to_MO() {
 int i, ijkl;
  int l=0;
  int s=0;
  int p, q, r, pq, rs, pqrs;
  TEI_MO = new double[5000];

  for(i=0,ijkl=0; i < nao; i++) {
    for(int j=0; j <= i; j++) {
      for(int k=0; k <= i; k++) {
        for(l=0; l <= (i==k ? j : k); l++,ijkl++) {
          for(p=0; p < nao; p++) {
            for(q=0; q < nao; q++) {
              pq = INDEX(p,q);
              for(r=0; r < nao; r++) {
                for(s=0; s < nao; s++) {
                  rs = INDEX(r,s);
                  pqrs = INDEX(pq,rs);
                  TEI_MO[ijkl] += coeff_in(p,i) * coeff_in(q,j) * coeff_in(r,k) * coeff_in(s,l) * TEI[pqrs];
                }
              }
            }
          }
        }
      }
    }
  } 
  return *TEI_MO;
}

// Translate integrals to spin-orbital basis 
Tensor4D CCSD::MO_to_spin() {
  double value1, value2;
  TEI_spin.resize(nsm, nsm, nsm, nsm);
  int pr,qs,ps,qr,prqs,psqr;
  for(int p=0; p < nsm; p++) {
    for(int q=0; q < nsm; q++) {
      for(int r=0; r < nsm; r++) {
        for(int s=0; s < nsm; s++) {
          pr = INDEX(p/2,r/2); // to say its in mo 0 or mo 1
          qs = INDEX(q/2,s/2);
          ps = INDEX(p/2,s/2);
          qr = INDEX(q/2,r/2);
          prqs = INDEX(pr,qs);
          psqr = INDEX(ps,qr);
          value1 = TEI_MO[prqs] * (p%2 == r%2) * (q%2 == s%2); // saying the contribution of spin ie if doin up and spin down equals 0
          value2 = TEI_MO[psqr] * (p%2 == s%2) * (q%2 == r%2);
          TEI_spin(p,q,r,s) = value1-value2;

        } 
      }
    }
  }
  return TEI_spin;
}

// creating spin orbital fock matrix
Matrix CCSD::spin_Fock() { 
  Matrix Hcore_spin(nsm,nsm);
  for (int p=0;p<nsm; p++) {
    for (int q=0;q<nsm;q++) {
      Hcore_spin(p,q) = 0;
    }
  }
  Matrix tr(nao,nao);
  tr = coeff_in.transpose();
  for (int i=0; i<nao;i++) {
    for (int j=0;j<nao;j++) {
      for (int p=0;p<nsm;p++) {
        for (int q=0;q<nsm;q++) {
          Hcore_spin(p,q) += tr(p/2,i)*Hcore(i,j)*coeff_in(j,q/2);
        }
      }
    }
  }
  //cout << Hcore_spin << endl;
  Fock_spin.resize(nsm,nsm);
  double sum=0;
  for (int p=0;p<nsm;p++) {
    for (int q=0;q<nsm;q++) {
      Fock_spin(p,q) = Hcore_spin(p,q)* (p%2 == q%2);
      for (int m=0;m<nsm_occ;m++) {
        Fock_spin(p,q) += TEI_spin(p,m,q,m);
      }
    }
  } 
  return Fock_spin;
}

//creates matrix of eigenvalues of spin fock matrix which is equivalent to spin orbital energies
Matrix CCSD::spin_energies() {
  Eigen::SelfAdjointEigenSolver<Matrix> solve(Fock_spin);
  Matrix evec = solve.eigenvectors();
  Matrix eval = solve.eigenvalues();
  energies.resize(nsm,nsm);
  for (int i=0; i<nsm; i++) {
    for (int j=0;j<nsm;j++) {
      if (i==j) {

        energies(i,j) = eval(i);
      }
      else {
        energies(i,j)=0;
      }
    }
  }
  return energies;
}

//calculates Emp2 energy
double CCSD::Emp2() {
  int tia=0;
  emp2=0.0;
  for (int i=0; i<nsm_occ;i++) {
      for (int j=0;j<nsm_occ; j++) {
          for (int a=nsm_occ;a<nsm;a++) {
        for (int b=nsm_occ;b<nsm;b++) {
          double tiajb = TEI_spin(i,j,a,b)/(energies(i,i)+energies(j,j)-energies(a,a)-energies(b,b));
          emp2 += TEI_spin(i,j,a,b)*tiajb;
        }
      }
    }
  }
  emp2 = .25*emp2;
  return emp2;
}

//creates a matrix of the initial t1 values (all 0)
Matrix CCSD::t1_initial() {
  t1_in.resize(nsm_occ, nsm_vir);
  t1_in.setZero(nsm_occ, nsm_vir);
  return t1_in;
}

//creates 4D tensor of the initial t2 values
Tensor4D CCSD::t2_initial() {
  t2_in.resize(nsm_occ, nsm_occ, nsm_vir, nsm_vir);
  double den;
  for ( int i = 0; i < nsm_occ; i++) {
    for ( int j=0; j < nsm_occ; j++) {
      for (int a=0; a < nsm_vir; a++) {
        for (int b=0; b < nsm_vir; b++) {
          den = Fock_spin(i,i) + Fock_spin(j,j) - Fock_spin(a+nsm_occ,a+nsm_occ) - Fock_spin(b+nsm_occ,b+nsm_occ);
          t2_in(i,j,a,b) = TEI_spin(i,j,a+nsm_occ,b+nsm_occ)/den;
          if (abs(t2_in(i,j,a,b)) < 3.3E-16) {
            t2_in(i,j,a,b)=0;
          }
          else {
            t2_in(i,j,a,b) = t2_in(i,j,a,b);
          }
        }
      }
    }
  }
  return t2_in;
}

Tensor4D CCSD::Tau_tilda(Matrix t1, Tensor4D t2) {
  tau_tilda.resize(nsm_occ, nsm_occ, nsm_vir, nsm_vir);
  tau_tilda = t2;
  for (int i=0; i < nsm_occ; i++) {
    for (int j=0; j < nsm_occ; j++) {
      for (int a=0; a < nsm_vir; a++) {
        for (int b=0; b < nsm_vir; b++) {
          tau_tilda(i,j,a,b) += 0.5 *t1(i,a)*t1(j,b);
					tau_tilda(i,j,a,b) -= 0.5 *t1(i,b)*t1(j,a);
          if (abs(tau_tilda(i,j,a,b))<1E-15) {
            tau_tilda(i,j,a,b)=0;
          }
          else {
            tau_tilda(i,j,a,b)=tau_tilda(i,j,a,b);
          }
        }
      }
    }
  }
  return tau_tilda;
}

Matrix CCSD::F_ae(Matrix t1, Tensor4D tau_tilda) {
  Fae.resize(nsm_vir,nsm_vir);
  Fae.setZero(nsm_vir,nsm_vir);
  int Kronecker=0;
  double value1=0, value2=0, value3=0, value4=0;
  for (int a=0;a<nsm_vir;a++) {
    for (int e=0;e<nsm_vir;e++) {
      int Kronecker=0;
      double value1=0, value2=0, value3=0, value4=0;
      value1 = Fock_spin(a+nsm_occ,e+nsm_occ);
  if (a==e) {
    Kronecker = 1;
  }
  else {
    Kronecker=0;
  }
  value1 = (1-Kronecker)*value1;
  for (int m=0; m<nsm_occ;m++) {
    value2 += Fock_spin(m,e+nsm_occ)*t1(m,a);
    for (int f=0;f<nsm_vir;f++) {
      value3 += t1(m,f)*TEI_spin(m,a+nsm_occ,f+nsm_occ,e+nsm_occ); 
      for (int n=0;n<nsm_occ;n++) {
        value4 += tau_tilda(m,n,a,f)*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
      }
    }
  }
  Fae(a,e) = value1 - .5*value2 + value3 - .5*value4;
    }
  }
  return Fae;
}

Matrix CCSD::F_mi(Matrix t1, Tensor4D tau_tilda) {
  Fmi.resize(nsm_occ,nsm_occ);
  Fmi.setZero(nsm_occ,nsm_occ);
  for (int m=0;m<nsm_occ;m++) {
    for (int i=0;i<nsm_occ;i++) {
      int Kronecker=0;
      double value1=0, value2=0, value3=0, value4=0;
      if (i==m) {
        Kronecker=1;
      }
      else {
        Kronecker=0;
      }
      value1 = (1-Kronecker)*Fock_spin(m,i);
      for (int e=0; e<nsm_vir;e++) {
        value2 += t1(i,e)*Fock_spin(m,e+nsm_occ);
        for (int n=0; n<nsm_occ; n++) {
          value3 += t1(n,e)*TEI_spin(m,n,i,e+nsm_occ);
          for (int f=0;f<nsm_vir;f++) {
            value4 += tau_tilda(i,n,e,f)*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
          }
        }
      }
  Fmi(m,i) = value1 + .5*value2 + value3 + .5*value4;
    }
  }
  return Fmi;
}


Matrix CCSD::F_me(Matrix t1) { 
  Fme.resize(nsm_occ,nsm_vir);
  Fme.setZero(nsm_occ,nsm_vir);
  for (int m=0;m<nsm_occ;m++) {
    for (int e=0;e<nsm_vir;e++) {
      double value1=0;
      for (int n=0; n<nsm_occ;n++) {
        for (int f=0; f<nsm_vir;f++) {
          value1 += t1(n,f)*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
        }
      }
      Fme(m,e) = Fock_spin(m,e+nsm_occ)+value1;
    }
  }
  return Fme;
}

Tensor4D CCSD::Tau(Tensor4D t2, Matrix t1) {
  tau.resize(nsm_occ,nsm_occ,nsm_vir,nsm_vir);
  tau = t2;
  for (int i=0; i < nsm_occ; i++) {
    for (int j=0; j < nsm_occ; j++) {
      for (int a=0; a< nsm_vir; a++) {
        for (int b=0; b<nsm_vir ; b++) {
          tau(i,j,a,b) += t1(i,a)*t1(j,b); 
					tau(i,j,a,b) -= t1(i,b)*t1(j,a);
          if (abs(tau(i,j,a,b))<1E-15) {
            tau(i,j,a,b)=0;
          }
          else {
            tau(i,j,a,b)=tau(i,j,a,b);
          }
        }
      }
    }
  }

  return tau;
}

Tensor4D CCSD::W_mnij(Matrix t1, Tensor4D tau) {
  Wmnij.resize(nsm_occ,nsm_occ,nsm_occ,nsm_occ);
  for (int m=0;m<nsm_occ;m++) {
    for (int n=0;n<nsm_occ;n++) {
      for (int i=0;i<nsm_occ;i++) {
        for (int j=0;j<nsm_occ;j++) {
          double value1=0, value2=0, value3=0, value4=0;
          value1 = TEI_spin(m,n,i,j);
          for (int e=0; e<nsm_vir; e++) {
            value2 += t1(j,e)*TEI_spin(m,n,i,e+nsm_occ);
            value3 +=t1(i,e)*TEI_spin(m,n,j,e+nsm_occ);
            for (int f=0; f<nsm_vir; f++) {
              value4 += tau(i,j,e,f)*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
            }
          }
          Wmnij(m,n,i,j) = value1 + value2 - value3 + .25*value4;
        }
      }
    }
  }
  return Wmnij;
}

Tensor4D CCSD::W_abef(Matrix t1, Tensor4D tau) {
  Wabef.resize(nsm_vir,nsm_vir,nsm_vir,nsm_vir);
  for (int a=0;a<nsm_vir;a++) {
    for (int b=0;b<nsm_vir;b++) {
        for (int e=0;e<nsm_vir;e++) {
          for (int f=0;f<nsm_vir;f++) {
            double value1=0, value2=0, value3=0, value4=0;
            value1 = TEI_spin(a+nsm_occ,b+nsm_occ,e+nsm_occ,f+nsm_occ);
            for (int m=0;m<nsm_occ;m++) {
              value2 += t1(m,b)*TEI_spin(a+nsm_occ,m,e+nsm_occ,f+nsm_occ);
              value3 += t1(m,a)*TEI_spin(b+nsm_occ,m,e+nsm_occ,f+nsm_occ);
              for (int n=0;n<nsm_occ;n++) {
                value4 += tau(m,n,a,b)*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
              }
            }
                Wabef(a,b,e,f) = value1-value2+value3+.25*value4;
          }
        }
      }
    }
  return Wabef; 
}

Tensor4D CCSD::W_mbej(Matrix t1, Tensor4D t2) {
  Wmbej.resize(nsm_occ,nsm_vir,nsm_vir,nsm_occ);
  double value1, value2, value3, value4;
  for (int m=0;m<nsm_occ;m++) {
    for (int b=0;b<nsm_vir;b++) {
      for (int e=0;e<nsm_vir;e++) {
        for (int j=0;j<nsm_occ;j++) {
          value1=0,value2=0,value3=0,value4=0;
          value1 = TEI_spin(m,b+nsm_occ,e+nsm_occ,j);
          for (int f=0;f<nsm_vir;f++) {
            value2 += t1(j,f)*TEI_spin(m,b+nsm_occ,e+nsm_occ,f+nsm_occ);
          }
          for (int n=0;n<nsm_occ;n++) {
            value3 += t1(n,b)*TEI_spin(m,n,e+nsm_occ,j);
            for (int f=0;f<nsm_vir;f++) {
              value4 += (.5*t2(j,n,f,b)+t1(j,f)*t1(n,b))*TEI_spin(m,n,e+nsm_occ,f+nsm_occ);
            }
          }
          Wmbej(m,b,e,j) = value1 + value2 - value3 - value4;
        }
      }
    }
  }
  return Wmbej;
}
//functions above(ie tau_tilda, Fae,Wmnij etc.) written to make code more efficient so you arent calculating each time for updated amplitudes
Matrix CCSD::update_T1(Matrix t1_old, Tensor4D t2) {
  updated_T1.resize(nsm_occ, nsm_vir);
  updated_T1.setZero(nsm_occ,nsm_vir);
  for (int i=0;i<nsm_occ;i++) {
    for (int a=0; a<nsm_vir; a++) {
      double value1=0,value2=0,value3=0,value4=0,value5=0,value6=0;
      for (int e=0;e<nsm_vir;e++) {
        value1 += t1_old(i,e)*Fae(a,e);
      }
        for (int m=0; m<nsm_occ; m++) {
          value2 += t1_old(m,a)*Fmi(m,i);
          for (int e=0; e<nsm_vir; e++) {
            value3 += t2(i,m,a,e)*Fme(m,e);
            for (int f=0; f<nsm_vir; f++) {
              value5 += t2(i,m,e,f)*TEI_spin(m,a+nsm_occ,e+nsm_occ,f+nsm_occ);
            }
            for (int n=0;n<nsm_occ;n++) {
              value6 += t2(m,n,a,e)*TEI_spin(n,m,e+nsm_occ,i);
            }
          }
        }
        for (int n=0; n<nsm_occ;n++) {
          for (int f=0; f<nsm_vir;f++) {
            value4 += t1_old(n,f)*TEI_spin(n,a+nsm_occ,i,f+nsm_occ);
          }
        }
      
      updated_T1(i,a) = Fock_spin(i,a+nsm_occ) + value1 - value2 + value3 - value4 - .5*value5 - .5*value6;
      double Ds = Fock_spin(i,i) - Fock_spin(a+nsm_occ,a+nsm_occ);
      updated_T1(i,a) /= Ds;
    }
  }
  return updated_T1;
}

Tensor4D CCSD::update_T2( Matrix t1_old, Tensor4D t2_old) {
	double value1, value2, value3, value4, value5, value6, value7, value8;
  double value1P, value;
  updated_T2.resize(nsm_occ,nsm_occ,nsm_vir,nsm_vir);
  for (int i=0; i< nsm_occ; i++) {
    for (int j=0; j< nsm_occ; j++) {
      for (int a=0; a<nsm_vir; a++) {
        for (int b=0; b<nsm_vir; b++){
					value=0;
					value2 = 0;
					value5 = 0;
					value7 = 0;
					for (int e=0; e< nsm_vir;e++) {
						value1 = 0;
						value1P = 0;
						for (int m=0; m < nsm_occ; m++) {
							value1 += t1_old(m,b)*Fme(m, e);
							value1P += t1_old(m,a)*Fme(m, e);
						}
						value2 += t2_old(i,j,a,e)*(Fae(b, e) - 0.5*value1)- t2_old(i,j,b,e)*(Fae(a, e) - 0.5*value1P);
						for (int f=0; f< nsm_vir; f++) {
							value5 += tau(i,j,e,f)*Wabef(a,b,e,f);
						}	
						value7 += t1_old(i,e)*TEI_spin(a+nsm_occ,b+nsm_occ,e+nsm_occ,j) - t1_old(j,e)*TEI_spin(a+nsm_occ,b+nsm_occ,e+nsm_occ,i);
					}
					
					value3 = 0;
					value4 = 0;
					value6 = 0;
					value8 = 0;	
					for (int m=0; m< nsm_occ; m++) {
						double v3a = 0;
						double v3b = 0;	
						for ( int e=0; e< nsm_vir; e++) {
							v3a += t1_old(j,e)*Fme(m, e);
							v3b += t1_old(i,e)*Fme(m, e);
						}
						value3 += t2_old(i,m,a,b)*(Fmi(m, j) + 0.5*v3a)- t2_old(j,m,a,b)*(Fmi(m,i)) + 0.5*v3b;
						for (int n=0; n < nsm_occ; n++) {
							value4 += tau(m,n,a,b)*Wmnij(m, n, i, j);
						}
						for (int e=0; e< nsm_vir; e++) {
							value6 += t2_old(i,m,a,e)*Wmbej(m, b, e, j) - t1_old(i,e)*t1_old(m,a)*TEI_spin(m,b+nsm_occ,e+nsm_occ,j);
							value6 -= t2_old(i,m,b,e)*Wmbej(m, a, e, j) - t1_old(i,e)*t1_old(m,b)*TEI_spin(m,a+nsm_occ,e+nsm_occ,j);
							value6 -= t2_old(j,m,a,e)*Wmbej(m, b, e, i) - t1_old(j,e)*t1_old(m,a)*TEI_spin(m,b+nsm_occ,e+nsm_occ,i);
							value6 += t2_old(j,m,b,e)*Wmbej(m, a, e, i) - t1_old(j,e)*t1_old(m,b)*TEI_spin(m,a+nsm_occ,e+nsm_occ,i);
						}
				    value8 += t1_old(m,a)*TEI_spin(m,b+nsm_occ,i,j) - t1_old(m,b)*TEI_spin(m,a+nsm_occ,i,j);
					}
					double Dijab = Fock_spin(i,i) + Fock_spin(j,j) - Fock_spin(a+nsm_occ,a+nsm_occ) - Fock_spin(b+nsm_occ,b+nsm_occ);
					value = TEI_spin(i,j,a+nsm_occ, b+nsm_occ);
					updated_T2(i,j,a,b) = value + value2 - value3 + value4/2 + value5/2 + value6 + value7 - value8;
					updated_T2(i,j,a,b) /= Dijab;
				}
      }
    }
  } 
		//cout << update_taibj << endl;
	return updated_T2;
}
//calculates CCSD energy
double CCSD::Ecc(Matrix updated_t1, Tensor4D t2_updated) {
  double Ecc=0;
  double value1=0,value2=0,value3=0;
  for (int i=0; i<nsm_occ;i++) {
    for (int a=0; a<nsm_vir;a++) {
      for (int j=0;j<nsm_occ;j++) {
        for (int b=0;b<nsm_vir;b++) {
          value2 += TEI_spin(i,j,a+nsm_occ,b+nsm_occ)*t2_updated(i,j,a,b);
          value3 += TEI_spin(i,j,a+nsm_occ,b+nsm_occ)*updated_t1(i,a)*updated_t1(j,b);
        }
      }
      value1 += Fock_spin(i,a+nsm_occ)*updated_t1(i,a);
    }
  }
  Ecc = value1 + .25*value2 + .5*value3;
  return Ecc; 
}

Tensor3D CCSD::Error(int it, Matrix T1i, Matrix T1i1, Tensor4D T2i, Tensor4D T2i1) {
  T1i.resize(1,nsm_occ*nsm_vir);
  T1i1.resize(1,nsm_occ*nsm_vir);
  Err.resize(2*nsm,nsm_occ*nsm_vir+1,nsm_occ*nsm_vir);
  int sum=0;
  if (it-2 <= 6) {
   for (int i = 0; i < nsm_occ; i++) {
     for (int a = 0; a < nsm_vir; a++) {
      for (int j = 0; j < nsm_occ; j++) {
          for (int b = 0;b < nsm_vir; b++) {
            int ia = INDEX(i,a);
            int jb = INDEX(j,b);
            Err(it-2,ia,jb) = T2i1(i,j,a,b) - T2i(i,j,a,b);
            if (i==0 && a==0) {
            Err(it-2,nsm_occ*nsm_vir+1,jb) = T1i1(i,a)-T1i(i,a); 
            }
          }
        }
      }
    }
  }
  else {
      for (int i = 0; i < nsm_occ; i++) {
     for (int a = 0; a < nsm_vir; a++) {
      for (int j = 0; j < nsm_occ; j++) {
          for (int b = 0;b < nsm_vir; b++) {
            int ia = INDEX(i,a);
            int jb = INDEX(j,b);
            Err(it-9,ia,jb) = T2i1(i,j,a,b) - T2i(i,j,a,b);
            if (i==0 && a==0) {
            Err(it-9,nsm_occ*nsm_vir+1,jb) = T1i1(i,a)-T1i(i,a);
            }
          }
        }
      }
    }
  }
  return Err;
}


Tensor3D CCSD::Sol_Diis(int it) {
  int val = it-1;
  B.resize(it,it);
  B.setZero(it,it);
  for (int x=0;x<it;x++) {
  for (int y=0;y<it;y++) {
  for (int i=0; i< nsm_occ; i++) {
    for (int a=0;a<nsm_vir;a++) {
      for (int j=0; j< nsm_occ; j++) {
        for (int b=0;b<nsm_vir;b++) {
          int ia = INDEX(i,a);
          int jb = INDEX(j,b);
          B(x,y) += Err(it,ia,jb)*Err(val,ia,jb);
          if (i==0 && a==0) {
            B(x,y)+=Err(it,nsm_occ*nsm_vir+1,jb)*Err(val,nsm_occ*nsm_vir+1,jb);
          }
          B(x,val) = -1;
          B(y,val) = -1;
          B(val,x) = -1;
          B(val,y) = -1;
        }
      }
    }
  }
  }
  }
    B(val,val)=0;
    Eigen::VectorXd C;
		C.setZero(it);
		C(it-1) = -1;
		//cout << "C" << endl;
		//cout << C << endl; 
    Eigen::HouseholderQR<Matrix> mat(B);
		Eigen::VectorXd ans;
        ans.setZero(it);
        ans = mat.solve(C);
        //cout << "ans" << ans << endl;
    disTen.resize(it,nsm_occ*nsm_vir,nsm_occ*nsm_vir+1);
	for (int i=0; i < it; i++) {
    for (int j=0; j < nsm_occ; j++) {
          for (int a=0;a<nsm_vir;a++) {
            int ja = INDEX(j,a);
            for (int k=0; k <nsm_occ; k++) {
              for (int b=0;b<nsm_vir;b++) {
          int kb = INDEX(k,b);
				disTen(it-2,ja,kb) += ans(i)*Err(it-2,ja,kb);
        if (j==0 && a==0) {
        disTen(i,nsm_occ*nsm_vir+1,kb) += Err(i, nsm_occ*nsm_vir+1,kb);
        }
              }
            }
    }
    }
  }

return disTen;
}
double CCSD::DisEcc(int it) {
  double diisEcc=0;
  double value1=0,value2=0,value3=0;
  for (int i=0; i<nsm_occ;i++) {
    for (int a=0; a<nsm_vir;a++) {
      for (int j=0;j<nsm_occ;j++) {
        for (int b=0;b<nsm_vir;b++) {
          int ia = INDEX(i,a);
          int jb = INDEX(j,b);
          value2 += TEI_spin(i,j,a+nsm_occ,b+nsm_occ)*disTen(it,ia,jb);
          value3 += TEI_spin(i,j,a+nsm_occ,b+nsm_occ);
          value1 += Fock_spin(i,a+nsm_occ);
          if (i==0 && a==0) {
          value3 += disTen(it,nsm_occ*nsm_vir+1,jb);
          value1 += disTen(it,nsm_occ*nsm_vir+1,jb)*disTen(it,nsm_occ*nsm_vir+1,jb);
          }
        }
      }
    }
  }
  diisEcc = value1 + .25*value2 + .5*value3;
  return diisEcc; 
}


CCSD::CCSD(int number, int number2) {
  nao = number;
  nocc = number2;
    
  TEI = new double[5000];
  FILE *el_in;
  el_in = fopen("/home/smullen22/IntroProjects/ProgrammingProjects/Project#05/input/h2o/STO-3G/eri.dat", "r");


  int u, v, l, s, uv, ls, uvls;
  double value;
    
  while (fscanf(el_in, "%d %d %d %d %lf", &u, &v, &l, &s, &value) != EOF) {
    u=u-1;
    v=v-1;
    l=l-1;
    s=s-1;
    uv = INDEX(u,v);
    ls = INDEX(l,s);
    uvls = INDEX(uv, ls);
    TEI[uvls] = value;
  }
  fclose(el_in); 
}


CCSD::~CCSD() {}



