#ifndef NON_OBJ_PROPAGATOR_H
#define NON_OBJ_PROPAGATOR_H

#include <vector>
#include <complex>
#include <cmath>

// my headers
#include "Constants.hpp"

// my preprocessor defs
#define delta_alpha 13.2 // 1.14*aupolarizability in O2 // 13.2 for CO2 // (6.98374727360818*aupolarizability) //1.96 cubic angstroms fo rCO_2, // .696 for n2  // ~7 in atomic units for o2
// delta alpha 13.2 a.u. for co2 from ref Yukikazu Itikawa, JPC Ref. Data v31 p749 (2002)
#define REAL(z,i) ((z)(2*i))
#define IMAG(z,i) ((z)(2*i+1))

typedef std::vector< std::complex<double> > cvec_t;

using namespace std;

typedef std::vector<double> record_t;
typedef std::vector< record_t > data2d_t;

template <typename T>
inline T PI(void){return boost::math::constants::pi<T>();}
template <typename T>
inline T two_PI(void){return boost::math::constants::two_pi<T>();}
template <typename T>
inline T root_PI(void) { return boost::math::constants::root_pi<T>(); }
template <typename T>
inline T half_PI(void) { return boost::math::constants::half_pi<T>(); }


void hc2ampphase(std::vector< std::vector< double * > > &data, const unsigned ntsteps);
void hc2ampphase_print(std::vector< std::vector< double * > > &data,const unsigned nsamples,std::ofstream &outstream);
void hc2ampphase_print(std::vector< std::vector< double * > > &data,const unsigned nsamples,std::string & basename);
void unwrap(const unsigned samples,double * array);
void removeslope(const unsigned samples,double * array);

struct PARAMS { // this is also crazy, redo this by implementing classes
  double *strengths, *omegas, *t0s, *Ctaus, *phases;
  int npulses;
  int m;
  int dim;
  int *jPtr, *vibsPtr; 
  double *ejPtr;
  double *pjPtr;
  double *pjnewPtr;
  int maxj;
  int sizej, nvibs, newvibsize, currentv;
  double *aajmPtr,*bbjmPtr,*ccjmPtr;
  int jstart;
  int dimref;
  int ntsteps;
  double kTinau;
  double *vibensPtr, *pvPtr;
  double tstepsize;
  data2d_t * legendresPtr;
  record_t * thetasPtr;
  data2d_t * distroPtr;
};

// proto-typing functions //  This is crazy, fix this with some classes
void freepropstep(double * t,const double * dt, double y[],void *paraPtrvoid);
bool nearpulse(const double t,PARAMS * paraPtr);
bool inpulse(const double t,PARAMS * paraPtr);
bool inpulse(const double t,PARAMS * paraPtr,double *FF);
bool inpulse(const double t,PARAMS * paraPtr,double *FF,double *dFFdt);

int func (double t,const double y[],double f[],void *paraPtrvoid);
int jac(double t,const double y[],double *dfdy,double dfdt[],void *paraPtrvoid);
void sqrnormalizey(double *y,PARAMS *paraPtr);
void addtosignaldistro(double *y,double *signal,double *imsignal,int tind,PARAMS *paraPtr);
void addtosignal(double *y,double *signal,double *imsignal,int tind,PARAMS *paraPtr);
void passtosignal(double *signal,double *imsignal,PARAMS *paraPtr);
void passtosignaldistro(double *signal,double *imsignal,PARAMS *paraPtr);
void addtopjnew(const double *y,PARAMS *paraPtr);
void passtopjnew(PARAMS *paraPtr);
void setrealj(int *realj,const int *i,PARAMS *paraPtr);
void samplecoupling(double *y,PARAMS *paraPtr);
void samplefield(double t,double FF);

void printLegendres(std::string & file,record_t & th,data2d_t & legs);
void computeLegendres(record_t & th, data2d_t & legs,PARAMS *paraPtr);

void setjejvecs(PARAMS *paraPtr);
void setvibsvibens(PARAMS *paraPtr);
void setelemsjm(const int m,const int sizej,double *aajmPtr,double *bbjmPtr,double *ccjmPtr);
void inity(double *y,PARAMS *paraPtr);
void setpjvec(PARAMS *paraPtr);
void setpvvec(PARAMS *paraPtr);

void rotationalenergies_oco(PARAMS *paraPtr);
void rotationalenergies_nno(PARAMS *paraPtr);
void rotationalenergies_oo(PARAMS *paraPtr);
void rotationalenergies_ii(PARAMS *paraPtr);
void rotationalenergies_nn(PARAMS *paraPtr);
void rotationalenergies_nn_nodistortion(PARAMS *paraPtr);

void vibrationalenergies_oco(PARAMS *paraPtr);
void vibrationalenergies_nno(PARAMS *paraPtr);
void vibrationalenergies_oo(PARAMS *paraPtr);
void vibrationalenergies_ii(PARAMS *paraPtr);
void vibrationalenergies_nn(PARAMS *paraPtr);

int getnewsize(const int oldsize,const double *vec);
void sumnormvec(const int size,double *vec);
void scalevec(const int size,double *vec,double scale);
double sumvec(const int size,const double *vec);

void print2col(const int sizej, const double *vec1, const double *vec2);
void print2col(const int sizej, const int *intvec,  const double *doublevec);

int getmaxdim(PARAMS *paraPtr);

int writepops(PARAMS *paraPtr,string pjstartname,string pjnewname);

#endif
