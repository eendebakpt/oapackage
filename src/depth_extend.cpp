#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

#include "extend.h"
#include "lmc.h"

#include "mathtools.h"
#include "strength.h"
#include "tools.h"
#include "arraytools.h"
#include "arrayproperties.h"
#include "depth_extend.h"



#ifdef DOOPENMP
#include "omp.h"
#endif

#ifndef DOOPENMP
// disable omp pragma warnings
#ifdef WIN32
#pragma warning (disable : 4068 )
#else
// http://choorucode.com/2014/03/11/how-to-ignore-specific-warning-of-gcc/
//#pragma GCC diagnostic ignored "-Wno-unknown-pragmas"
//#pragma GCC diagnostic ignored "-Wno-pragmas"
#pragma GCC diagnostic ignored "-Wenum-compare"
#endif
#endif

using namespace std;

// extend arrays with dynamic selection
int extend_array_dynamic ( const arraylist_t &alist, arraydata_t &adx, OAextend &oaextend, arraylist_t &earrays, std::vector<std::vector<double> > &edata, int kfinal, double Dfinal, int directcheck, int verbose )
{
	const int rs=1;
	long ntotal=0;
	long n=0;
	long nlmc=0;


	double t0=get_time_ms();
	double DmaxDiscard = -1;
	long nmaxrnktotal = 0;	// total number of arrays with max rank generated

	// TODO: if we remove the static code in LMC, then we can make this loop parallel
	//printf("oaanalyse: openmp %d\n", omp_get_max_threads() );
	//#pragma omp parallel for
	for ( size_t i=0; i<alist.size();  i++ ) {
		dextend_t dextend;
		dextend.directcheck=directcheck;
		array_link al = alist[i];
		arraylist_t earrays0;
		std::vector<std::vector<double> > edata0;

		extend_array_dynamic ( al, adx, oaextend, earrays0, edata0, dextend, kfinal, Dfinal, rs, verbose );
		#pragma omp critical(extenddynamic)
		{
			appendArrays ( earrays0, earrays );
			edata.insert ( edata.end(), edata0.begin(), edata0.end() );


			ntotal += dextend.ntotal;
			nlmc += dextend.nlmc;
			n += dextend.n;
			DmaxDiscard =  max ( dextend.DmaxDiscard, DmaxDiscard );
			nmaxrnktotal += dextend.nmaxrnktotal;

			double gainf = double ( n ) / ( i+1 );
			double dt=get_time_ms()-t0;
			int nc=al.n_columns;
			if ( verbose>=2 )
				printf ( "## extend_array_dynamic array %d/%d (N %d, ncols %d->%d, dt %.1f [s]): found %ld/%ld/%ld extensions (total %ld/%ld/%ld, gain factor %.1f, saving factor %.4f)\n", ( int ) i, ( int ) alist.size(), adx.N, nc, nc+1, dt, dextend.n, dextend.nlmc, dextend.ntotal, ( long ) n, ( long ) nlmc, ( long ) ntotal, gainf, double ( n ) /nlmc );
		}

		// append
	}

	if ( verbose>=1 ) {
		printf ( "## extend_array_dynamic: total %ld/%ld/%ld\n", n, nlmc, ntotal );
		printf ( "## extend_array_dynamic: DmaxDiscard %f\n", DmaxDiscard );
		printf ( "## extend_array_dynamic: nmaxrnktotal %ld\n", nmaxrnktotal );
	}
	return n;
}

/// convert a list of D values to loss factors
std::vector<double> Dvalues2lossfactor ( const dextend_t &dextend, double Ccurrent, int kn )
{
	size_t nn = dextend.Deff.size();

	std::vector<double> L ( nn );
	for ( size_t ii=0; ii<nn; ii++ ) {
		//if (ii%nprint==0)
		//  printf("array %d/%d\n", ii, nn);
		double Ci = Dvalue2Cvalue ( dextend.Deff[ii], kn );
		if ( Ccurrent<1e-15 )
			L[ii]=1;
		else
			L[ii]=Ci/Ccurrent;
	}
	return L;

}



void testTranspose ( const array_link &al )
{
	array_link alt = al.transposed();
	std::vector<double> x = alt.GWLP();
	{
		printf ( "  transposed GWP " );
		printf_vector<double> ( x, "%.3f " );
		cout << endl;
	}
	std::vector<GWLPvalue> dopgmaT = projectionGWLPs ( alt );
	{
		printf ( "  transposed dopgwp " );
		display_vector<GWLPvalue> ( dopgmaT );
		cout << endl;
	}

	symmetry_group sgT ( dopgmaT, 0 );
	sgT.show();
}




template <class Type>
/// count number of occurences of an element in a std::vector
int countOccurences ( const std::vector<Type> &v, Type val )
{
	int n=0;
	for ( size_t i=0; i<v.size(); i++ ) {
		if ( v[i]==val ) {
			n++;
		}
	}
	return n;
}


int extend_array_dynamic ( const array_link &input_array, arraydata_t &adx, OAextend &oaextend, arraylist_t &earraysout, std::vector<std::vector<double> > &edata, dextend_t &dextend, int kfinal, double Dfinal, int rs, int verbose )
{
	OAextend oaextendx= oaextend;
	arraylist_t extensions;

	const int nprint = 4000;
	const int nprint2=20000;
	int extensioncol=input_array.n_columns;

	dextend.nmaxrnktotal = 0;

	int directcheck=dextend.directcheck;

	oaextendx.checkarrays=directcheck;
	oaextendx.use_row_symmetry=rs;
	oaextendx.nLMC=500000;
	oaextendx.setAlgorithmAuto ( &adx );
	oaextendx.extendarraymode=OAextend::extendarray_mode_t::APPENDEXTENSION;	// save memory by only appending the partial array

	if ( kfinal<input_array.n_columns+1 ) {
		printf ( "extend_array_dynamic: kfinal too small! kfinal %d, cols in extension %d\n", kfinal, input_array.n_columns+1 );
	}

	arraylist_t earrays;

	// make extensions
	if ( verbose>=3 )
		printf ( "extend_array_dynamic: verbose %d\n", verbose );

	const int lastcol=input_array.n_columns;
	int nex=extend_array ( input_array,  &adx, extensioncol, earrays, oaextendx );
	if ( verbose>=3 )
		printf ( "extend_array_dynamic: call to extend_array done\n" );

	int nn = earrays.size();
	dextend.resize ( nn );

	array_link tmparray ( input_array.n_rows, input_array.n_columns+1,-1 );
	std::copy ( input_array.array, input_array.array+input_array.n_columns*input_array.n_rows, tmparray.array );

	#pragma omp critical
	{
		if ( directcheck ) {
			cprintf ( verbose>=2, "extend_array_dynamic: directcheck done (nr. extensions %d, earrays %d)", nex, (int)earrays.size() );
		} else {
            int nlmc = 0;
			if ( verbose>=2 )
				printf ( "   LMC check %d values, extensioncol %d\n", nn, extensioncol );


			/* calculate extensions */
			LMCreduction_t reduction ( &adx );
			arraydata_t *ad = new arraydata_t ( adx.s, adx.N, adx.strength, extensioncol+1 );


			for ( size_t ii=0; ii< ( size_t ) nn; ii++ ) {
				if ( verbose>=2 ) {
					if ( ( ( ii%nprint==0 ) && verbose>=3 )  || ( ( ii%nprint2==0 ) && verbose==2 ) ) {
						printf ( "  LMC calculation for array %d/%d (current nlmc %d)\n", ( int ) ii, nn, nlmc );
						fflush ( stdout );
					}
				}


				tmparray.setcolumn ( lastcol, earrays[ii] );

				reduction.lastcol=-1;
				reduction.init_state=COPY; // INIT

				lmc_t v = LMCcheck ( tmparray, *ad, oaextendx, reduction );
				if ( verbose>=4 )
					printf ( "lmc_t %d, col %d\n", v, reduction.lastcol );
				dextend.lmctype[ii]=v;
				dextend.lastcol[ii]=reduction.lastcol;
				if ( v==LMC_MORE )
					nlmc++;
			}
			delete ad;

		}
	}

	if ( verbose>=2 )
		printf ( "   analysing %d values (D-efficiency value)\n", nn );
	flush_stdout();

	/* calculate D-efficiency values */
	for ( unsigned int ii=0; ii< ( size_t ) nn; ii++ ) {
		if ( ( ( ii%nprint==0 ) && verbose>=3 )  || ( ( ii%nprint2==0 ) && verbose==2 ) ) {
			printf ( "  D-efficiency calculation for array %d/%d \n", (int)ii, nn );
			flush_stdout();
		}
		tmparray.setcolumn ( lastcol, earrays[ii] );

		dextend.Deff[ii] = dextend_t::NO_VALUE;
		if ( dextend.Dcheck== dcalc_mode::DCALC_ALWAYS || dextend.lmctype[ii]==LMC_MORE ) {
			if (adx.is2level() )
				dextend.Deff[ii] = Defficiency ( tmparray, verbose>=2 );
			else {
				array_link altmp = hstack(tmparray, earrays[ii] );
				dextend.Deff[ii] =altmp.Defficiencies ()[0];
			}
		}

	}

	for ( unsigned int ii=0; ii< ( size_t ) nn; ii++ ) {
		if ( dextend.lmctype[ii]==LMC_MORE ) { // only use the LMC arrays
			if ( dextend.Deff[ii]>0 ) {
				dextend.nmaxrnktotal +=1;
			}
		}
	}
	/* filtering */

	int k = extensioncol;	// number of columns of the array being extended
	int kn=k+1;

	std::vector<double> rtmp;
	int r = array2rank_Deff_Beff ( input_array, &rtmp, 0 );
	double Dcurrent = rtmp[1];
	double Ccurrent = Dvalue2Cvalue ( Dcurrent, k );

	cprintf ( verbose>=3, "   current C value %f (%d columns), rank %d\n", Ccurrent, k, r );
	flush_stdout();
	std::vector<double> L = Dvalues2lossfactor ( dextend, Ccurrent, kn );

	double Lmax=vectormax<double> ( L, 1 );
	cprintf ( verbose>=3, "   loss factor: max %f at %ld\n", Lmax );


	/// calculate filter
	dextend.DefficiencyFilter ( Dfinal, k, kfinal, Lmax, verbose );


	//int ngood =std::accumulate ( dextend.filter.begin(),dextend.filter.end(),0 );
	int ngoodlmc = countOccurences ( dextend.lmctype, LMC_MORE );

	std::vector<int> ctype = dextend.filterArrays ( input_array, earrays, earraysout, edata,verbose );
	int ngoodcombined =std::accumulate ( ctype.begin(),ctype.end(),0 );

	double DmaxDiscard = -1;

	for ( size_t idx=0; idx<ctype.size(); idx++ ) {
		if ( ctype[idx] ) {
			std::vector<double> w;
			w.push_back ( dextend.Deff[idx] );
			w.push_back ( Lmax );
			edata.push_back ( w );
		} else {
			DmaxDiscard = max ( DmaxDiscard, dextend.Deff[idx] );
		}
	}

	dextend.ntotal=nn;
	dextend.nlmc = ngoodlmc;
	dextend.n = ngoodcombined;
	dextend.DmaxDiscard = DmaxDiscard;

	if ( verbose>=3 )
		printf ( "extend_array_dynamic: done\n" );
	flush_stdout();


	return 0;

}
