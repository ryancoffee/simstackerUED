/* compile with: (to get gsl linked on slac unix)
   g++ -Wall -I/usr/local/include -c -DHAVE_INLINE non_obj_propagator.cpp
   g++ -L/usr/local/lib non_obj_propagator.o -lgsl -lgslcblas -lm
 */

// standard includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <ctime>
#include <string>
//#include <initializer_list> // this requires -std=c++11 in the CFLAGS
#include <vector>
#include <assert.h>
#include <algorithm>
#include <iterator>


// gsl includes
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_rng.h>

// fftw3 includes
#include <fftw3.h>


// Boost imports // 
#include <boost/math/constants/constants.hpp>

#include "./non_obj_propagator.hpp"
//#include "./data2d.hpp"

// OpenCV includes //
//#include "opencv2/core/core.hpp"

// kmeans includes //
//#include "./kmeans.hpp"

// kmlocal include //
//#include "kmlocal/KMlocal.h"

using namespace std;

#define MIN_POWSPEC 1e-3 // this is really dirty...

int main(int argc, char *argv[]) {

	double PI = boost::math::constants::pi<double>();
	std::cout << "testing boost::constants\tpi = \t" << PI << std::endl;

	// getting cmd line args
	if (argc < 9) {
		cerr << "Syntax is $ main maxj T(K) twinstart(ps) tend(ps) tstep (ps) numpulses strength(10^16W/cm^2) decay(>0) width(fs) delay(fs)" << endl;
		cerr << "Writes files pj_nn.dat, cossq_nn.dat, field_nn.dat" << endl;
		cerr << "./non_obj_propagator 60 300 11 12 .010 1 .01 1 30 12000" << std::endl;
		return 1;
	}

	// generate filenames //
	string filename,filetail="_",pjstarttail="_";
	std::string datadir("data_n2o/");
	std::string molstring("n2o");
	filetail += argv[2];
	filetail += "K_";
	filetail += argv[6];
	filetail += "pls_";
	filetail += argv[7];
	filetail += "inten_";
	filetail += argv[9];
	filetail += "fswidth_";
	filetail += argv[10];
	filetail += "fsdelay.dat";

	pjstarttail += argv[2];
	pjstarttail += "K.dat";

	// Start the clock //
	time_t start,end,rawtime;
	double diff;
	time(&start);
	filename = datadir + molstring + ".kicklog";
	filename += filetail;
	ofstream logout(filename.c_str(),ios::app);
	if (!logout){
		cerr << "Couldn't append to " << filename << endl;
		return 1;
	}

	time(&rawtime);
	logout << ctime(&rawtime) << " $ " ;
	for (int i=0;i<argc;i++){
		logout << argv[i] << " ";
	}
	logout << endl;


	// setting up local constants
	const double abstol=1e-10;
	const double reltol=1e-8;        // relative tolerance reltol = 1e-3 in matlab default
	const double kTinau = static_cast<double>(atof(argv[2]))*kb/Eh;
	const int maxj = static_cast<int>(atoi(argv[1]));
	const int sizej = maxj+1;
	const int nvibs = 4; // set this to reasonable number after debugging
	const int npulses = static_cast<int>(atoi(argv[6]));

	// setting up local variables //
	int j[sizej],vibs[nvibs],newvibsize;
	newvibsize = nvibs;
	double vibens[nvibs],pv[nvibs];
	double ej[sizej],pj[sizej],pjnew[sizej];//,pjmnew[sizej][2*sizej-1];
	//  double sqrt2jp1[sizej+2];
	double aajm[sizej],bbjm[sizej],ccjm[sizej];
	double strengths[npulses],Ctaus[npulses],t0s[npulses],omegas[npulses],phases[npulses];

	for (int i=0;i<npulses;i++){
		strengths[i]= delta_alpha*auenergy/Eh * gsl_pow_2(aufor10PW) * static_cast<double>(atof(argv[7]));
		//strengths[i] *= gsl_pow_int(static_cast<double>(atof(argv[8])), i) ;
		Ctaus[i]=M_SQRTPI * static_cast<double>(atof(argv[9]))/fsPau/2.0; // this is 1/e half-width, input is assumed to be 1/e full-width of a gaussian.  The sqrt(pi) gives a proper cos^2 pulse length that integrates to the same area as the gaussian.
		t0s[i]= i * static_cast<double>(atof(argv[10]))/fsPau;
		omegas[i]=0.0;// hc/800/Eh; // wavelength assumed to be 800nm, or 0 if non-resonant
		phases[i]=0.0 * pi;
	}

	double *aajmPtr=aajm;
	double *bbjmPtr=bbjm;
	double *ccjmPtr=ccjm;

	// setting up pointers //
	int *jPtr=j,*vibsPtr=vibs;
	double *pvPtr=pv,*ejPtr=ej,*vibensPtr=vibens,*pjPtr=pj,*pjnewPtr=pjnew;

	PARAMS params;
	void *paraPtrvoid=&params; 

	params.strengths = strengths;
	params.Ctaus = Ctaus;
	params.t0s = t0s;
	params.phases = phases;
	params.omegas = omegas;
	params.npulses = npulses;


	params.kTinau=kTinau;
	params.jPtr = jPtr;
	params.ejPtr = ejPtr;
	params.pjPtr = pjPtr;
	params.pjnewPtr = pjnewPtr;
	params.aajmPtr=aajmPtr;
	params.bbjmPtr=bbjmPtr;
	params.ccjmPtr=ccjmPtr;

	params.vibsPtr=vibsPtr;
	params.vibensPtr=vibensPtr;
	params.pvPtr=pvPtr;
	params.nvibs=nvibs;
	params.newvibsize=newvibsize;

	// Init vibs //
	setvibsvibens(&params);
	setpvvec(&params);

	params.maxj=maxj;
	params.sizej=sizej;

	const double twinstart = (double)(atof(argv[3]))*1e3/fsPau;
	const double tend = (double)(atof(argv[4]))*1e3/fsPau;
	if (tend<twinstart){
		cerr << "Whoops, twinstart is greater than tend... avoid a segmentation fault, fix this!" << endl;
		return 1;
	};

	const double tstepsize = (double)(atof(argv[5]))*1e3/fsPau;
	double tstart = twinstart;
	for (int i=0;i<npulses;i++){
		tstart = GSL_MIN_DBL(params.t0s[i] - params.Ctaus[i],tstart);
	}

	// ------------ setting up the time steps and ntsteps to accomodate the wavelet transform need for 2^j length ------- //
	double tmp  =  (( tend - twinstart )/ tstepsize) +1;
	unsigned ntsteps;
	unsigned ntpower;
	ntpower = (unsigned)((log(tmp)/log(2.)) + 1);
	ntsteps = pow(2, ntpower);
	
	params.ntsteps = ntsteps;
	params.tstepsize = tstepsize;
	clog << twinstart*fsPau << " fs: " << tstepsize*fsPau << " fs: " << (twinstart+(tstepsize*ntsteps))*fsPau << " fs in " << ntsteps << " steps." << endl;

	// ------------- setting up time and signal vectors ---------------- //

	unsigned tstep=0;
	double times[2*ntsteps],signal[2*ntsteps],imsignal[2*ntsteps];//,cossq[2]={0.0,0.0};
	for(tstep=0;tstep<ntsteps;tstep++){
		signal[tstep]=0.0;
		imsignal[tstep]=0.0;
		times[tstep]=twinstart+tstep*tstepsize;
	}
	tstep=0;


	// ---------- setting up the ODE system ---------- //

	// initializing evolution system
	double passt=tstart;                     // - delay -  (1/e half-widths) of pulse1  
	//  clog << "t starts at " << passt*fsPau << "fs" << endl;


	// --------- Init pjnew vec --------//
	for (int i=0;i<sizej;i++){
		pjnewPtr[i] = 0.0;
	}

	//  params.newvibsize=1; // use to evaluate only v=0
	for (int v=0;v<params.newvibsize;v++) {    // ultimately loop on v
		params.currentv=v;
		setjejvecs(&params);
		setpjvec(&params);


		double pop=1.0; //for testing the preservation of square magnitude
		const double poptol=0.5; // 	used in if (pop > 1+poptol || pop < 1-poptol) ... 

		for(int m=0;m<sizej-2;m++) {                                               // ultimately loop on m
			//    clog << "working on m = " << m << "  "; 

			params.m=m;

			setelemsjm(m,sizej,aajmPtr,bbjmPtr,ccjmPtr);   // ,coeffaamPtr,coeffbbmPtr,coeffccmPtr,aaj0Ptr,bbj0Ptr,ccj0Ptr,sqrt2jp1Ptr);

			// ------- now I want to make a function to set dim ----------- //
			int dim=getmaxdim(&params); //(sizej-m)*2;// (jwin*2+1)*2;  // here's where we remove jwin, it seems dim is twice too big since only even(odd) Js couple

			params.dim=dim;
			params.ejPtr=ejPtr;
			const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;  // fast step
			//    const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;  // smooth step
			gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,dim);
			gsl_odeiv_control *c = gsl_odeiv_control_y_new(abstol,reltol);
			gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(dim);
			gsl_odeiv_system sys = {func,jac,(size_t)dim,paraPtrvoid};

			double h = .1/fsPau;                                   // initial step size
			double y[dim];        

			for(int jstart = 0;jstart<sizej-m;jstart++) {           // ultimately loop on j
				params.jstart = jstart;
				if (jstart+m<10 || pjPtr[jstart+m]>5e-4) {

					double t=tstart;
					inity(y,&params);



					int status=GSL_SUCCESS;
					unsigned tind;
					tind=0;
					for(tind=0;tind<ntsteps;tind++){
						if (!nearpulse(t,&params)){             // somehow do a freepropstep in here

							freepropstep(&t,&tstepsize,y,&params);

						} else {                              // solve ODEs for coeffs

							while (t<times[tind] && status==GSL_SUCCESS){
								status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,times[tind],&h,y);
								sqrnormalizey(y,&params);
								//	    fieldout << t*fsPau/1e3 << "\t" << params.field << endl;
								if (status !=GSL_SUCCESS){
									cerr << "evolve status: " << gsl_strerror(status) << endl;
									break;
								}
							}
						}
						addtosignal(y,signal,imsignal,tind,&params);

						passt = t;
					}
					//	clog << "exited while loop, status: " << status << " = " << gsl_strerror(status) << endl;

					addtopjnew(y,&params);
					if (pop> 1+poptol || pop < 1-poptol){
						cerr << "not preserving square magnitude for j" << params.jstart+params.m << ", m" << params.m << endl;
						cerr << "square magnitude = " << pop << endl;
						return 1;
					}

					if (m==0 && jstart==m && v==0 ){ // Test coupling for j,m=0,0 and j=26,m=13
						samplecoupling(y,&params);
					}
				} else {
					passtosignal(signal,imsignal,&params);
					passtopjnew(&params);
				}
			} // close j loop 

			// ---------- free the ODE memory ---------- //

			// needs to be in the m-loop since m determines the dimension... maybe later see if always keeping track of sizej equations is better.
			gsl_odeiv_evolve_free(e);
			gsl_odeiv_control_free(c);
			gsl_odeiv_step_free(s);

			//    clog << "done." << endl;    

		} // close m loop

		time(&end);
		diff=difftime(end,start);
		clog << "just finished v = " << params.currentv << " at " << diff << "s." << endl;
		//  clog << "t ends at " << passt*fsPau << "fs" << endl;
	} // close v loop





	// stop the clock
	time(&end);
	diff=difftime(end,start);
	clog << ntsteps << " timesteps, " << diff << " s runtime.\a" << endl;
	logout << ntsteps << " timesteps, " << diff << " s runtime." << endl;

	string pjstartname,pjnewname;
	pjstartname = datadir + molstring + ".pjstart";
	pjstartname += pjstarttail;
	pjnewname = datadir + molstring + ".pj";
	pjnewname += filetail;

	int wrote = writepops(&params,pjstartname,pjnewname);
	if (wrote !=0)
		return 1;


	//open output file
	filename = datadir + molstring + ".cossq";
	filename += filetail;
	ofstream cossqout(filename.c_str(),ios::out);
	std::string collectfilename(datadir);
	collectfilename += molstring + ".collect";

	ofstream cossqcollect(collectfilename.c_str(),std::ios::app);
	if (!cossqout){
		cerr << filename << " could not be opened" << endl;
		return 1;
	}
	for(unsigned i=0;i<ntsteps;i++){
		double psdelay;
		psdelay = times[i]*fsPau/1e3;
		cossqout << psdelay << " " << signal[i] << " " << imsignal[i] << endl;
		cossqcollect << signal[i] << "\t";
	}
	cossqcollect << "\n";
	cossqout.close();
	cossqcollect.close();

	return 0;
}

void unwrap(const unsigned samples,double * array){
	for (unsigned i=1;i<samples;i++){
		double difference;
		difference = array[i] - array[i-1];
		if (difference > M_PI){
			for (unsigned j=i;j<samples;j++)
				array[j] -= 2.*M_PI;
		} else if (difference < -M_PI) {
			for (unsigned j=i;j<samples;j++)
				array[j] += 2.*M_PI;
		}
	}
	double ref = array[0];
	for (unsigned i=0;i<samples;i++)
		array[i] -= ref;
}
void removeslope(const unsigned samples,double * array){
	array[samples/2] -= M_PI*samples/2;
	for (unsigned i=1;i<samples/2;i++){
		array[i] -= M_PI*i;
		array[samples-i] += M_PI*i;
	}

}

void hc2ampphase(std::vector< std::vector< double * > > &data,const unsigned nsamples){
	const unsigned ens = data.size();
	const unsigned chans = data[0].size();
	//std::cout << ens << "\t" << chans << std::endl;
	assert(chans==3);
	double re,im;
	const unsigned srcind = 0, ampind = 1, phaseind = 2;
	for(unsigned e=0;e<ens;e++){
		re = data[e][srcind][0];
		data[e][ampind][0] = fabs(re)/sqrt(nsamples);
		data[e][phaseind][0] = atan2(0.0,re);
		re = data[e][srcind][nsamples/2];
		data[e][ampind][nsamples/2] = fabs(re)/sqrt(nsamples);
		data[e][phaseind][nsamples/2] = atan2(0.0,re);
		for (unsigned k=1;k<nsamples/2;k++){
			// the FFTW does not normalize //
			re = data[e][srcind][k];
			im = data[e][srcind][nsamples-k];
			data[e][ampind][k] = data[e][ampind][nsamples-k] = sqrt( pow((double)re,(int)2) + pow((double)im,(int)2) / nsamples);
			data[e][phaseind][k] = data[e][phaseind][nsamples-k] = atan2(im,re);
		}
	}
	for (unsigned e=0;e<ens;e++){
		if (data[e][phaseind][0])
		unwrap(nsamples,data[e][phaseind]);
		removeslope(nsamples,data[e][phaseind]);
	}

}
void hc2ampphase_print(std::vector< std::vector< double * > > &data,const unsigned nsamples,std::ofstream &outstream){
	for(unsigned e=0;e<data.size();e++){
		for (unsigned k=0;k<nsamples;k++){
			outstream << e << "\t" << k << "\t" << data[e][0][k] << "\t" << data[e][1][k] << "\t" << (data[e][2][k]/M_PI) << "\n";
		}
		outstream << "\n";
	}

}
void hc2ampphase_print(std::vector< std::vector< double * > > &data,const unsigned nsamples,std::string & basename){
	static std::ofstream outstream;
	std::string filename;
	filename = basename + ".amp";
	outstream.open(filename.c_str(),std::ios::out);
	for(unsigned e=0;e<data.size();e++){
		for (unsigned k=0;k<nsamples;k++){
			outstream << data[e][1][k] << "\t";
		}
		outstream << "\n";
	}
	outstream.close();
	filename = basename + ".phase";
	outstream.open(filename.c_str(),std::ios::out);
	for(unsigned e=0;e<data.size();e++){
		for (unsigned k=0;k<nsamples;k++){
			outstream << data[e][2][k] << "\t";
		}
		outstream << "\n";
	}
	outstream.close();

}



	// ---------- member functions ---------- //



	void print2col(const int size, const int *intvec,  const double *doublevec){
		for(int i=0;i<size;i++){
			cout << *intvec << " " << *doublevec << endl;
			intvec++;
			doublevec++;
		}
		intvec-=size;
		doublevec-=size;
	}

	void print2col(const int size, const double *vec1, const double *vec2){
		for(int i=0;i<size;i++){
			cout << *vec1 << " " << *vec2 << endl;
			vec1++;
			vec2++;
		}
		vec1-=size;
		vec2-=size;
	}

	double sumvec(const int size,const double *vec){
		double sum=0;
		for(int i=0;i<size;i++){
			sum+=*vec;
			vec++;
		}
		vec -=size;
		return sum;
	}

	void scalevec(const int size,double *vec,double scale){
		for(int i=0;i<size;i++){
			*vec*=scale;
			vec++;
		}
		vec-=size;
	}

	void sumnormvec(const int size,double *vec){
		double sum=0;
		for(int i=0;i<size;i++){
			sum+=*vec;
			vec++;
		}
		vec -=size;
		scalevec(size,vec,1/sum);
	}

	// Initialize vib and rot vectors

	void setvibsvibens(PARAMS *paraPtr){
		// first set vibs via vibsPtr
		for(int vv=0;vv<paraPtr->nvibs;vv++){
			paraPtr->vibsPtr[vv]=vv;
		}
		// now set vibens via vibensPtr
		//vibrationalenergies_oco(paraPtr);
		vibrationalenergies_nno(paraPtr);
		//vibrationalenergies_oo(paraPtr);
	}

	void setjejvecs(PARAMS *paraPtr){
		// first set j via jPtr
		for(int i=0;i<paraPtr->sizej;i++){
			paraPtr->jPtr[i]=i;
		}
		// now set ej via ejPtr
		//rotationalenergies_oco(paraPtr);
		rotationalenergies_nno(paraPtr);
		//rotationalenergies_oo(paraPtr);
	}

	void setpvvec(PARAMS *paraPtr){//mark
		double *pvPtr = paraPtr->pvPtr;
		double frac;
		for(int i=0;i<paraPtr->nvibs;i++){
			frac=(paraPtr->vibensPtr[i]) / paraPtr->kTinau;
			pvPtr[i]=exp(-frac);
		}
		sumnormvec(paraPtr->nvibs,pvPtr);
		paraPtr->newvibsize=getnewsize(paraPtr->nvibs,pvPtr);
		sumnormvec(paraPtr->newvibsize,pvPtr);
		clog << "relevent vib pops = ";
		for(int i=0;i<paraPtr->newvibsize;i++){
			clog << pvPtr[i] << " ";
		}
		clog << endl;
	}


	int getnewsize(const int oldsize,const double *vec){
		int newsize=0;
		double floor = 1e-3;
		for (int i=0;i<oldsize;i++){
			if(*vec>floor){
				newsize++;
				vec++;
			}
		}
		vec-=oldsize;
		return newsize;
	}


	void setpjvec(PARAMS *paraPtr){
		const int *jPtr = paraPtr->jPtr;
		double *pjPtr = paraPtr->pjPtr;
		const double *ejPtr = paraPtr->ejPtr;
		const int sizej = paraPtr->sizej;

		double frac, jj, mul, qj;
		for(int i=0;i<sizej;i++){

			((*jPtr%2)==0 ? mul=1 : mul=1); // I co2 think this is same as O2 but maybe opposite if electronic is symmetric for CO2 // even odd intensity ratio = mul 2:1 for n2, 7:5 for i2, 
			// 0:1 for o2 I think, for N2O this should be an even 1:1 ratio

			jj=static_cast<double>(*jPtr);
			frac = *ejPtr/(paraPtr->kTinau); // atomic energy units
			if(fabs(frac)<.1){
				*pjPtr=(2*jj+1)*(1+gsl_expm1(-frac));
			}else{
				*pjPtr=(2*jj+1)*exp(-frac);
			}
			*pjPtr*=mul;
			pjPtr++;
			ejPtr++;
			jPtr++;
		}
		pjPtr -=sizej;
		ejPtr -=sizej;
		jPtr -=sizej;
		qj=sumvec(sizej,pjPtr);
		scalevec(sizej,pjPtr,1/qj);
	} 

	void setelemsjm(const int m,const int sizej,double *aajmPtr,double *bbjmPtr,double *ccjmPtr){
		//  clog << "started setting 3jms ... " << endl;  
		static int jj;
		for(jj=0;jj<sizej;jj++){
			aajmPtr[jj]=0.0;
			bbjmPtr[jj]=0.0;
			ccjmPtr[jj]=0.0;
		}

		for(jj=m;jj<sizej;jj++) { // these coefficienst come from Arfken and Webber 4th ed. p753, eq.12.189 applied twice.
			aajmPtr[jj] = static_cast<double>((jj-m+1)*(jj+m+1)*(2*jj-1) + (2*jj+3)*(jj+m)*(jj-m))/static_cast<double>((2*jj-1)*(2*jj+1)*(2*jj+3));
			if (jj>1) { // this term is 0 when j=0 m=0, j=1 m=0, j=1 m=1, and j=1 m=-1
				bbjmPtr[jj] = sqrt(static_cast<double>((jj-m)*(jj+m)*(jj-m-1)*(jj+m-1)) / static_cast<double>((2*jj+1)*(2*jj-3)*gsl_pow_2(2*jj-1)) );  
			}
			ccjmPtr[jj] = sqrt(static_cast<double>((jj-m+1)*(jj+m+1)*(jj-m+2)*(jj+m+2)) / static_cast<double>((2*jj+1)*(2*jj+5)*gsl_pow_2(2*jj+3)));  
		}

	}

	// Calculate vibrational energies:

	void vibrationalenergies_ii(PARAMS *paraPtr){
		const int nterms=10;
		int v=0,n=0;
		double vv=0.;
		const double cc[]={214.5481, -0.616259, 7.507e-5, -1.263643e-4, 6.198129e-6, -2.0255975e-7, 3.9662824e-9, -4.6346554e-11, 2.9330755e-13, -7.61000e-16};
		const double *ccPtr=cc;
		for(v=0;v<paraPtr->nvibs;v++){
			vv=static_cast<double>(v);
			paraPtr->vibensPtr[v]=0;
			for(n=0;n<nterms;n++){
				paraPtr->vibensPtr[v]+=*ccPtr*gsl_pow_int(vv+0.5,n+1);
				ccPtr++;
			}
			ccPtr-=nterms;
			paraPtr->vibensPtr[v] /= icmPau;
		}
	}

	void vibrationalenergies_nn(PARAMS *paraPtr){
		const int n=min(11,paraPtr->nvibs);
		const double ens[]={1175.5, 3505.2, 5806.5, 8079.2, 10323.3, 12538.8, 14725.4, 16883.1, 19011.8, 21111.5, 23182.0};
		for(int v=0;v<n;v++){
			paraPtr->vibensPtr[v]=ens[v];
			paraPtr->vibensPtr[v] /=icmPau;
		}
	}
	void vibrationalenergies_oo(PARAMS *paraPtr){
		// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
		const int n=min(4,paraPtr->nvibs);
		const double env_0 = (1556.38991/2.0);
		const double dens[]={1556.38991, 1532.86724, 1509.5275};
		double ens[4];
		ens[0]=env_0;
		ens[1] = ens[0] + dens[0];
		ens[2] = ens[1] + dens[1];
		ens[3] = ens[2] + dens[2];
		for(int v=0;v<n;v++){
			paraPtr->vibensPtr[v]=ens[v];
			paraPtr->vibensPtr[v] /= icmPau;
		}
	}


	// Calculate rotational energies:
	void rotationalenergies_nno(PARAMS *paraPtr)
	{
		// From Bohlin et al., J. Raman Spect. v43 p604 (2012)
		// units are icm as usual
		std::vector<double> Bv {0.419011 , 0.419177 , 0.419969 , 0.419920 , 0.420125 , 0.420126 , 0.417255 , 0.419583 , 0.421079 , 0.420667 , 0.420671 , 0.417464 , 0.418372 , 0.415559 , 0.420618 , 0.420768 , 0.420772 , 0.421218 , 0.418147 , 0.418530 , 0.418531 , 0.415605 };
		std::vector<double> Dv {1.76e-7 , 1.78e-7 , 1.79e-7 , 2.49e-7 , 1.19e-7 , 1.18e-7 , 1.72e-7 , 2.11e-7 , 2.17e-7 , 1.61e-7 , 1.68e-7 , 1.74e-7 , 1.71e-7 , 1.75e-7 , 3.96e-7 , 0.17e-7 , 2.16e-7 , 2.72e-7 , 2.43e-7 , 1.20e-7 , 1.75e-7 , 1.63e-7 };
		std::vector<double> Hv {0.16e-13, -0.17e-13, -0.17e-13 , 29.55e-13, -29.50e-13, 0.95e-13, 1.46e-13, 12.22e-13, -3.59e-13, -9.91e-13, 30.15e-13, 1.07e-13, 2.17e-13, -0.13e-13, 140.28e-13, -142.92e-13, 4.83e-13, 1712.20e-13, 25.90e-13, -26.98e-13, 2.44e-13, 5.68e-13};

		assert(Bv.size() == Dv.size());
		assert(Bv.size() == Hv.size());
		const int vv = min((int)Bv.size(),paraPtr->currentv);
		for (unsigned i=0;i<paraPtr->sizej;i++){
			double jj= (double)(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1)-Dv[vv]*gsl_pow_2(jj)*gsl_pow_2(jj+1) + Hv[vv]*gsl_pow_3(jj)*gsl_pow_3(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}
	void vibrationalenergies_nno(PARAMS *paraPtr)
	{
		// From Bohlin et al., J. Raman Spect. v43 p604 (2012)
		// units are icm as usual
		std::vector<double> ens {0.0, 588.767870,588.767870 ,1168.13230 ,1177.74467 ,1177.74467 ,1284.90334 ,1749.06523 ,1749.06515 ,1766.91238 ,1766.91224 ,1880.26574 ,1880.26574 ,2223.75676 ,2322.57308 ,2331.12151 ,2331.12145 ,2356.25242 ,2461.99644 ,2474.79870 ,2474.79865 ,2563.33944};

		const int maxvibs = min(paraPtr->nvibs,(int)ens.size());
		for(unsigned v=0;v<maxvibs;v++){
			paraPtr->vibensPtr[v] = ens[v]/icmPau;
		}

	}
	void vibrationalenergies_oco(PARAMS *paraPtr)
	{
		// from Rothman and Young, J. Quonf. Spectrosc. Radia. Transfer Vol. 25, pp. 505-524. 1981
		std::vector<double> ens {0.0,667.379,1285.4087,1335.129,1388.1847,1932.472,2003.244,2076.855,2349.1433,	2548.373,2585.032,2671.146,2671.716,2760.735,2797.140,3004.012,3181.450,3240.564,3339.340,3340.501,3442.256,3500.590,3612.842,3659.271,3714.783	};
		const int maxvibs = min(paraPtr->nvibs,(int)ens.size());
		for(unsigned v=0;v<maxvibs;v++){
			paraPtr->vibensPtr[v] = ens[v]/icmPau;
		}
	}
	void rotationalenergies_oco(PARAMS *paraPtr)
	{
		// from Rothman and Young, J. Quonf. Spectrosc. Radia. Transfer Vol. 25, pp. 505-524. 1981
		std::vector<double> Bv {0.39021894,0.39064230,0.39048230,0.39167020,0.39018893,0.39073215,0.39238558,0.39041600,0.38714140,0.39110670,0.39193800,0.38954820,0.39308410,0.39153500,0.39059100,0.38759300,0.3910280,0.3926960,0.3900350,0.393908,0.3922100,0.3904610,0.38750493,0.38864000,0.38706227};
		std::vector<double> Dv {1.33373e-7,1.359e-7,1.57161e-7,1.389e-7,1.14952e-7,1.441e-7,1.403e-7,1.281e-7,1.33034e-7,1.7820e-7,1.390e-7,1.2630e-7,1.42e-7,1.44e-7,0.88e-7,1.349e-7,1.63e-7,1.51e-7,1.37e-7,1.44e-7,1.36e-7,1.14e-7,1.58150e-7,1.37445e-7,1.13570e-7};
		std::vector<double> Hv {0.16e-13,0.17e-13,2.33e-13,0.,1.91e-13,0.,0.,0.,0.17e-13,0.40e-13,0.,4.65e-13,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.74e-13,0.,1.10e-13};

		assert(Bv.size() == Dv.size());
		assert(Bv.size() == Hv.size());
		const int v = min(paraPtr->currentv,(int)Bv.size());
		for (unsigned i=0;i<paraPtr->sizej;i++){
			double j= (double)(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[v]*j*(j+1)-Dv[v]*gsl_pow_2(j)*gsl_pow_2(j+1) + Hv[v]*gsl_pow_3(j)*gsl_pow_3(j+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}

	}

	void rotationalenergies_oo(PARAMS *paraPtr){
		// from JOURNAL OF MOLECULAR SPECTROSCOPY 154,372-382 ( 1992) G. ROUILLE "High-Resolution Stimulated Raman Spectroscopy of O2"
		const int vv = min(3,paraPtr->currentv);
		const double Bv[] = {1.437676476, 1.42186454, 1.4061199, 1.39042};
		const double Dv[] = {4.84256e-6, 4.8418e-6, 4.8410e-6, 4.8402e-6};
		const double Hv = 2.8e-12;
		double jj;
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1)-Dv[vv]*gsl_pow_2(jj)*gsl_pow_2(jj+1) + Hv*gsl_pow_3(jj)*gsl_pow_3(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}
	void rotationalenergies_nn(PARAMS *paraPtr){
		// numbers come from Loftus and Kuprienie, J. Phys. Chem. ref. Data, Vol. 6, No. 1, 1977. p242
		const int vv = min(16,paraPtr->currentv);
		const double Bv[]={1.98957, 1.972, 1.9548, 1.9374, 1.9200, 1.9022, 1.8845, 1.8666, 1.8488, 1.8310, 1.8131, 1.7956, 1.7771, 1.7590, 1.7406, 1.7223};
		double Dv = 1e-6*5.75, jj;
		//  Dv*=10; // artificially increase the distortion rather than temperature.
		//  Dv *= 0.0; // artificially kill distortion
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1)-Dv*gsl_pow_2(jj)*gsl_pow_2(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}

	void rotationalenergies_nn_nodistortion(PARAMS *paraPtr){
		// numbers come from Loftus and Kuprienie, J. Phys. Chem. ref. Data, Vol. 6, No. 1, 1977. p242
		const int vv = min(16,paraPtr->currentv);
		const double Bv[]={1.98957, 1.972, 1.9548, 1.9374, 1.9200, 1.9022, 1.8845, 1.8666, 1.8488, 1.8310, 1.8131, 1.7956, 1.7771, 1.7590, 1.7406, 1.7223};
		double jj;
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			paraPtr->ejPtr[i] = Bv[vv]*jj*(jj+1);// in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}

	void rotationalenergies_ii(PARAMS *paraPtr){
		double Bv, Dv, vv, jj;
		vv = static_cast<double>(paraPtr->currentv);
		for (int i=0;i<paraPtr->sizej;i++){
			jj= static_cast<double>(paraPtr->jPtr[i]);
			Bv=3.7395e-2 - 1.2435e-4*(vv+0.5) + 4.498e-7*gsl_pow_2(vv+0.5) - 1.482e-8*gsl_pow_3(vv+0.5) - 3.64e-11*gsl_pow_4(vv+0.5);
			Dv=4.54e-9 + 1.7e-11*(vv+0.5) + 7e-12*gsl_pow_2(vv+0.5);
			paraPtr->ejPtr[i]=Bv*jj*(jj+1)-Dv*gsl_pow_2(jj)*gsl_pow_2(jj+1); // in cm^-1
			paraPtr->ejPtr[i] /= icmPau; // in atomic units
		}
	}


	int func (double t,const double y[],double f[],void *paraPtrvoid){
		PARAMS *paraPtr;
		paraPtr = (PARAMS *)paraPtrvoid;
		const int m=paraPtr->m;
		const int dim=paraPtr->dim;
		const int jstart = paraPtr->jstart;

		int realj;
		int *realjPtr=&realj;

		static double FF;
		FF = 0.0;

		static int i;

		for (i=0;i<dim;i++){
			f[i]=0.0;
		}

		if( inpulse(t,paraPtr,&FF) ){  // is coupling

			double aa,bb,cc;

			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				//      if (paraPtr->m==1 && paraPtr->jstart==0 && static_cast<int>(t) == 0){cout << *realjPtr << " ";}
				if (*realjPtr != -1){
					aa = (paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]) - (1+paraPtr->aajmPtr[*realjPtr])*FF;
					bb = -1.0 * (paraPtr->bbjmPtr[*realjPtr])*FF;
					cc = -1.0 * (paraPtr->ccjmPtr[*realjPtr])*FF;
					f[i] = aa*y[i+1];
					f[i+1] = -1.0*aa*y[i];
					if (i>0) {
						f[i] += bb*y[i-2+1];
						f[i+1] -= bb*y[i-2];
					}
					if (i<dim-2) {
						f[i] += cc*y[i+3];
						f[i+1] -= cc*y[i+2];
					}
				}
			}

		} else {  // no coupling

			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					f[i]= (paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])*y[i+1];
					f[i+1]= -(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])*y[i];
					// y' = E/i y; y = exp(-iEt)
				}
			}
		} 
		return GSL_SUCCESS;
	}

	int jac(double t,const double y[],double *dfdy,double dfdt[],void *paraPtrvoid){
		PARAMS *paraPtr;
		paraPtr = (PARAMS *)paraPtrvoid;
		const int m=paraPtr->m;
		const int dim = paraPtr->dim;
		const int jstart = paraPtr->jstart;

		static int realj;
		int *realjPtr=&realj;

		static int i;
		for (i=0;i<gsl_pow_2(dim);i++){
			dfdy[i]=0.0;
		}
		gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,dim,dim);
		gsl_matrix *mat = &dfdy_mat.matrix;


		for (int i=0;i<dim;i++){
			dfdt[i]=0.0;
		}

		double FF=0.0, dFFdt = 0.0;

		if( inpulse(t,paraPtr,&FF,&dFFdt) ){  // is coupling


			gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,dim,dim);
			gsl_matrix *mat = &dfdy_mat.matrix;


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					gsl_matrix_set(mat,i,i+1,(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])-(1+paraPtr->aajmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i,i-2+1,-(paraPtr->bbjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i,i+2+1,-(paraPtr->ccjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i+1,i,-((paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m])-(1+paraPtr->aajmPtr[*realjPtr])*FF));
					gsl_matrix_set(mat,i+1,i-2,(paraPtr->bbjmPtr[*realjPtr])*FF);
					gsl_matrix_set(mat,i+1,i+2,(paraPtr->ccjmPtr[*realjPtr])*FF);
				}
			}



			double dadt,dbdt,dcdt;


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					dadt=-(1+paraPtr->aajmPtr[*realjPtr])*dFFdt;
					dbdt=-(paraPtr->bbjmPtr[*realjPtr])*dFFdt;
					dcdt=-(paraPtr->ccjmPtr[*realjPtr])*dFFdt;

					dfdt[i] = dadt*y[i+1];
					dfdt[i+1] = -1.0*dadt*y[i];
					if (i>0) {
						dfdt[i] += dbdt*y[i-1];
						dfdt[i+1] -= dbdt*y[i-2];
					}
					if (i<dim-2) {
						dfdt[i] += dcdt*y[i+3];
						dfdt[i+1] -= dcdt*y[i+2];
					}

				}
			}
		} else { // no coupling


			for (i=0;i<dim;i+=2){
				setrealj(realjPtr,&i,paraPtr);
				if (*realjPtr != -1){
					gsl_matrix_set(mat,i,i+1,paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]);
					gsl_matrix_set(mat,i+1,i,-(paraPtr->ejPtr[*realjPtr]-paraPtr->ejPtr[jstart+m]));

				}
			}

		} 
		return GSL_SUCCESS;
	}


	void addtosignal(double *y,double *signal,double *imsignal,int tind,PARAMS *paraPtr){

		const int m=paraPtr->m;
		const int dim = paraPtr->dim;
		const int jstart = paraPtr->jstart;
		const int v = paraPtr->currentv;

		double realcossq=0.0, imagcossq=0.0;
		int realj;
		int *realjPtr=&realj;

		//  sqrnormalizey(y,paraPtr);

		for (int i=0;i<dim;i+=2){
			setrealj(realjPtr,&i,paraPtr);
			if (*realjPtr != -1){
				realcossq += ( gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]) ) * paraPtr->aajmPtr[*realjPtr];
				imagcossq += 0.0;
				if(i<(dim-2)){
					realcossq += (y[i] * y[i+2] + y[i+1] * y[i+3])  * paraPtr->ccjmPtr[*realjPtr];
					imagcossq += (-y[i+3] * y[i] + y[i+2] * y[i+1]) * paraPtr->ccjmPtr[*realjPtr];
				}
				if(i>2){
					realcossq += (y[i] * y[i-2] + y[i+1] * y[i-1]) * paraPtr->bbjmPtr[*realjPtr];
					imagcossq += (-y[i-1] * y[i] + y[i-2] * y[i+1]) * paraPtr->bbjmPtr[*realjPtr];
				}
			}
		}

		if(m==0){
			signal[tind] += realcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
			imsignal[tind] += imagcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
		} else {
			signal[tind] += 2*realcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
			imsignal[tind] += 2*imagcossq * paraPtr->pvPtr[v]  * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1);
		}
	}

// here we can add addtocos() rather than addtocossq().  That would include terms that look forward and backward only by one element //


void passtosignal(double *signal,double *imsignal,PARAMS *paraPtr){
	const int m = paraPtr->m;
	const int jstart = paraPtr->jstart;
	const int v = paraPtr->currentv;

	for (int tind=0;tind<paraPtr->ntsteps;tind++){
		if(m==0){
			signal[tind] += paraPtr->aajmPtr[jstart+m] * paraPtr->pvPtr[v] * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1); // same as adding 1/3 * pj
		} else {
			signal[tind] += 2 * paraPtr->aajmPtr[jstart+m] * paraPtr->pvPtr[v] * paraPtr->pjPtr[jstart+m]/(2*(jstart+m)+1); // same as adding 2/3 * pj
		}
	}
}

void addtopjnew(const double *y,PARAMS *paraPtr){
	const int m=paraPtr->m;
	const int v=paraPtr->currentv;
	const int dim = paraPtr->dim;
	const int jstart = paraPtr->jstart;
	double *pjnewPtr = paraPtr->pjnewPtr;
	const double *pjPtr = paraPtr->pjPtr;
	const double *pvPtr = paraPtr->pvPtr;
	double sqrmag;

	int realj;
	int *realjPtr=&realj;

	for (int i=0;i<dim;i+=2){
		setrealj(realjPtr,&i,paraPtr);
		sqrmag = gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
		if (m==0){
			pjnewPtr[*realjPtr] += sqrmag  * pvPtr[v] * pjPtr[jstart+m]/(2*(jstart+m)+1);
		} else {
			pjnewPtr[*realjPtr] += 2 * sqrmag * pvPtr[v] * pjPtr[jstart+m]/(2*(jstart+m)+1);
		}
	}
}

void passtopjnew(PARAMS *paraPtr){

	const int m=paraPtr->m;
	const int v=paraPtr->currentv;
	const int jstart = paraPtr->jstart;
	double *pjnewPtr = paraPtr->pjnewPtr;
	const double *pjPtr = paraPtr->pjPtr;
	const double *pvPtr = paraPtr->pvPtr;

	if (m==0){
		pjnewPtr[jstart+m] += pvPtr[v]* pjPtr[jstart+m]/(2*(jstart+m)+1);
	} else {
		pjnewPtr[jstart+m] += pvPtr[v]*2*pjPtr[jstart+m]/(2*(jstart+m)+1);
	}
}

void setrealj(int *realjPtr,const int *i,PARAMS *paraPtr){

	const int *m=&(paraPtr->m);
	const int *jstart = &(paraPtr->jstart);
	const int *maxj = &(paraPtr->maxj);

	*realjPtr = *i + *m; // + *jstart; is not right.
	*realjPtr += *jstart%2;

	if (*realjPtr<0 || *realjPtr>*maxj){
		*realjPtr=-1;
	}
}

void inity(double *y,PARAMS *paraPtr){
	const int *dim = &(paraPtr->dim);
	const int *jstart = &(paraPtr->jstart);
	//  const int *m = &(paraPtr->m);
	for(int i=0;i< *dim;i++){
		y[i]=0.0;                                // set ys to 0;
	}
	//  Now, we need to set the proper y[i] = 1.0;
	int indx;
	indx = 2 * ( *jstart / 2); // Grabs the integer part of j/2 so that jstart = 1 lands on index 0, and 3 on 2.
	y[indx] = 1.0;
}


void samplecoupling(double *y,PARAMS *paraPtr){
	int realj;
	clog << "\nlog10(population times 1e10):\n";
	for (int i=0;i<paraPtr->dim;i+=2){
		setrealj(&realj,&i,paraPtr); 
		double value = gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
		value *= 1e10;
		int dots = static_cast<int>(log10(value));
		clog << realj << ":\t";
		for (int i = 0; i<dots; i++)
			clog << ".";
		clog << "\n";
	}
	clog << "\n";
}

void samplefield(double t,double FF){
	cout << t*fsPau/1e3 << "\t" << FF << "\n";   
}

int writepops(PARAMS *paraPtr,string pjstartname,string pjnewname){
	//open output file
	ofstream out(pjstartname.c_str(),ios::out);
	if (!out){
		cerr << pjstartname << " could not be opened" << endl;
		return 1;
	}
	for (int i=0;i<paraPtr->sizej;i++){
		out << i << " " << paraPtr->pjPtr[i] << endl;
	}

	//open output file
	ofstream pjout(pjnewname.c_str(),ios::out);
	if (!pjout){
		cerr << pjnewname << " could not be opened" << endl;
		return 1;
	}
	double pjnewsum = 0;
	for (int i=0;i<paraPtr->sizej;i++){
		pjout << i << " " << paraPtr->pjnewPtr[i] << endl;
		pjnewsum +=  paraPtr->pjnewPtr[i];
	}

	//  clog << "total pop is " << pjnewsum << endl;

	return 0;
}

void sqrnormalizey(double *y,PARAMS *paraPtr){
	double pop=0.0;
	for (int i=0;i<paraPtr->dim;i+=2){
		pop +=  gsl_pow_2(y[i]) + gsl_pow_2(y[i+1]);
	}
	scalevec(paraPtr->dim,y,1.0/sqrt(pop));
}

int getmaxdim(PARAMS *paraPtr){
	const int sizej = paraPtr->sizej;
	const int m = paraPtr->m;
	int out = sizej - m;
	out += out%2;
	return out;
}


void freepropstep(double *t,const double * dt, double y[],void *paraPtrvoid){
	PARAMS * paraPtr = (PARAMS *)paraPtrvoid;
	int jstart = paraPtr->jstart;
	int m = paraPtr->m;
	gsl_vector_view yview = gsl_vector_view_array(y,paraPtr->dim);
	size_t cdim = (paraPtr->dim)/2;
	gsl_vector_view yrealview = gsl_vector_subvector_with_stride(&yview.vector,0,2,cdim);
	gsl_vector_view yimagview = gsl_vector_subvector_with_stride(&yview.vector,1,2,cdim);

	gsl_complex z;

	gsl_vector_view ejvec = gsl_vector_view_array(paraPtr->ejPtr,paraPtr->sizej);

	//    gsl_vector_view energiesPart = gsl_vector_subvector_with_stride( &ejvec.vector,
	//  								   m+(jstart%2),
	//  								   2,
	//  								   cdim );
	double phase;
	for (size_t i = 0; i< cdim-1; i++){
		z = gsl_complex_rect(gsl_vector_get(&yrealview.vector,i),
				gsl_vector_get(&yimagview.vector,i) );
		phase = gsl_vector_get( &ejvec.vector , i*2 +jstart%2+m );
		phase -= gsl_vector_get( &ejvec.vector, (jstart+m) );
		phase *= -*dt;
		z = gsl_complex_mul(z,gsl_complex_polar(1.0,phase));
		gsl_vector_set(&yrealview.vector,i,GSL_REAL(z));
		gsl_vector_set(&yimagview.vector,i,GSL_IMAG(z));
		// add phase = -(Ej-Ejstart)dt to coeffs
		// I got here, now I need to wrap the phase by the energy.
	}
	*t += *dt;
}

bool nearpulse(const double t,PARAMS * paraPtr){
	int p=0;
	for (p=0;p<paraPtr->npulses;p++){
		if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] - 3*paraPtr->tstepsize && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p] + 3*paraPtr->tstepsize){
			return true;
		}
	}
	return false;
}

bool inpulse(const double t,PARAMS * paraPtr){
	int p=0;
	for (p=0;p<paraPtr->npulses;p++){
		if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
			return true;
		}
	}
	return false;
}
bool inpulse(const double t,PARAMS * paraPtr,double *FF){
	bool isinpulse = false;
	int p=0;
	for (p=0;p<paraPtr->npulses;p++){
		if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
			isinpulse = true;
			*FF += paraPtr->strengths[p]* ( gsl_pow_2( cos(M_PI_2*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]) ) );
		}
	}
	return isinpulse;
}
bool inpulse(const double t,PARAMS * paraPtr,double *FF,double *dFFdt){
	bool isinpulse = false;
	int p=0;
	for (p=0;p<paraPtr->npulses;p++){
		if (t >= paraPtr->t0s[p]-paraPtr->Ctaus[p] && t <= paraPtr->t0s[p]+paraPtr->Ctaus[p]){
			isinpulse = true;
			*FF += paraPtr->strengths[p]* ( gsl_pow_2( cos(M_PI_2*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]) ) );
			*dFFdt += -paraPtr->strengths[p]/2 * ( M_PI/paraPtr->Ctaus[p] * sin(M_PI*(t-paraPtr->t0s[p])/paraPtr->Ctaus[p]));
		}
	}
	return isinpulse;
}



