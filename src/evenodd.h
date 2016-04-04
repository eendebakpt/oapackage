/*! \file evenodd.h
 *  \brief Contains functions to generate even-odd designs
 *
 */

#pragma once

#ifdef OADEBUG
#else
#define DOOPENMP
#endif

#include "arraytools.h"
//#include "arrayproperties.h"
#include "lmc.h"
#include "extend.h"


/// structure to write arrays to disk, thread safe
struct arraywriter_t {
public:

	// since depth_extend is a depth first approach we need to store arrays with a different number of columns
	std::vector<arrayfile_t*> afiles;

	bool writearrays;

	/// number of arrays written to disk
	int nwritten;
	// long narraysmax;

	int verbose;

public:
	arraywriter_t() {
		writearrays=true;
		verbose=1;
	};

	~arraywriter_t() {
		flush();
		closeafiles();
	}

	void flush() {
		//printf("arraywriter_t: flush()\n");
		for ( size_t i=0; i<afiles.size(); i++ ) {
			arrayfile_t *af = afiles[i];
			if ( af!=0 ) {
				//printf("arraywriter_t: flush() %d\n", i);
				#pragma omp critical
				af->updatenumbers();
				af->flush();
			}
		}
	}
	// write array to disk
	void writeArray ( const array_link &A ) {

		// writing arrays with multiple threads at the same time is not supported
		#pragma omp critical
		{
			size_t i = A.n_columns;
			if ( writearrays ) {
				if ( i<afiles.size() ) {
					afiles[i]->append_array ( A );
				} else {
					fprintf ( stderr, "depth_extend_t: writeArray: problem: array file for %d columns was not opened\n", ( int ) i );
				}
				nwritten++;
			}
		}
	}

	// write a list of arrays to disk
	void writeArray ( const arraylist_t &lst ) {
		// NOTE: for slow filesystems we might want to cache these results

//#pragma omp critical
		{
			for ( size_t j=0; j<lst.size(); j++ ) {
				const array_link &A = lst[j];
				writeArray ( A );
			}
		}
	}

	// initialize the result files
	void initArrayFiles ( const arraydata_t &ad, int kstart, const std::string prefix, arrayfilemode_t mode = ABINARY_DIFF ) {
		afiles.clear();
		afiles.resize ( ad.ncols+1 );
		nwritten=0;

		for ( size_t i=kstart; i<= ( size_t ) ad.ncols; i++ ) {
			arraydata_t ad0 ( &ad, i );
			std::string afile = prefix + "-" + ad0.idstr() + ".oa";
			if ( verbose>=3 )
				printf ( "depth_extend_t: creating output file %s\n", afile.c_str() );

			int nb = arrayfile_t::arrayNbits ( ad );
			afiles[i] = new arrayfile_t ( afile, ad.N, i, -1, mode, nb );
		}
	}

	/// return the total number arrays
	int nArraysWritten() const {
		return nwritten;
	}

public:
	void closeafiles() {
		for ( size_t i=0; i< afiles.size(); i++ )
			delete afiles[i];
		afiles.clear();

	}

};


/// structure containing current position in search tree
struct depth_path_t {
	std::vector<int> ncurr;	// vector with current position
	std::vector<int> nmax;	// vector with target
	std::vector<int> necols;	// number of extension columns
	std::vector<int> ngecols;	// number of good extension columns
	int depthstart;

	depth_path_t() {
	}

	void updatePositionGEC ( int k,  int goodextensioncols ) {
		ngecols[k]=goodextensioncols;

	}
	void updatePosition ( int k, int c, int m, int extensioncols, int goodextensioncols ) {
		ncurr[k]=c;
		nmax[k]=m;
		necols[k]=extensioncols;
		ngecols[k]=goodextensioncols;

	}
	void show ( int depth, int maxentries=8 ) const {
		for ( int i=depthstart; i<=depth; i++ ) {

			if ( ( i-depthstart ) ==maxentries ) {
				printf ( "... " );
				break;
			}

			printf ( "%d: %d/%d (%d->%d) ", i, ncurr[i], nmax[i], necols[i], ngecols[i] );
		}
		printf ( "\n" );

	}
	void init ( int ncols, int _depthstart=9 ) {
		ncurr.resize ( ncols+1 );
		nmax.resize ( ncols+1 );
		necols.resize ( ncols+1 );
		ngecols.resize ( ncols+1 );
		depthstart=_depthstart;
	}
};

/// structure to count and show number of arrays generated, the structure is thread safe
struct counter_t {

	std::vector<int> nfound;	// vector with number of arrays found

	counter_t ( int n ) {
		nfound.resize ( n+1 );
	}

	void addNfound ( int col, int num ) {
		#pragma omp atomic
		this->nfound[col]+=num;
	}

	long nArrays() const {
		long na= std::accumulate ( this->nfound.begin(),this->nfound.end(),0 );;
		return na;
	}
	void addNumberFound ( int n, int k ) {
		#pragma omp critical (DEXTEND_NFOUND)
		{
			this->nfound[k]+=n;
		}
	}

	void clearNumberFound() {
		#pragma omp critical
		{
			for ( size_t k=0; k<this->nfound.size(); k++ ) {
				this->nfound[k]=0;
			}
		}
	}

	void addNumberFound ( const counter_t &de ) {
		#pragma omp critical
		{
			for ( size_t k=0; k<this->nfound.size(); k++ ) {
				this->nfound[k]+=de.nfound[k];
			}
		}
	}


	/// show information about the number of arrays found
	inline void showcountscompact() const {
		#pragma omp critical
		{
			//printf("depth_extend: counts ");
			//for ( size_t i=ad->strength; i<= ( size_t ) ad->ncols; i++ ) {
			//   printf ( " %d\n", this->nfound[i] );
			// }
			printf ( "depth_extend: counts " );
			display_vector ( this->nfound );
			printf ( "\n" );
		}
	}

	/// show information about the number of arrays found
	inline void showcounts ( const arraydata_t &ad ) const {
		printf ( "--results--\n" );
		for ( size_t i=ad.strength; i<= ( size_t ) ad.ncols; i++ ) {
			printf ( "depth_extend: column %ld: found %d\n", i, this->nfound[i] );
		}
	}

	/// show information about the number of arrays found
	inline void showcounts ( const char *str, int first, int last ) const {
		printf ( "--results--\n" );
		for ( size_t i=first; i<= ( size_t ) last; i++ ) {
			printf ( "%s: column %ld: found %d\n", str, i, this->nfound[i] );
		}
	}
};





/** Helper structure for dynamic extension
 *
 * In this structure we keep track of pointers to valid column extensions
 *
 */
struct depth_extend_sub_t {
public:
	std::vector<int> lmctype;
	/// last column changed in lmc check
	std::vector<int> lastcol;

	std::vector<double> strengthcheck;
	std::vector<int> valididx;
	// param
	int verbose;


	depth_extend_sub_t ( int nn=0 ) : verbose ( 0 )  {
		resize ( nn );
	};

	void resize ( int nn ) {
		this->lmctype.resize ( nn );
		this->lastcol.resize ( nn );
		this->strengthcheck.resize ( nn );
		// this->tmp.resize(nn);

		std::fill ( this->lmctype.begin(), this->lmctype.begin() + nn, LMC_MORE );
		std::fill ( this->lastcol.begin(), this->lastcol.begin() + nn, -1 );
	}


	inline size_t n() const {
		return lmctype.size();
	}

	std::vector<int> updateExtensionPointers ( int extcol ) {
		if ( verbose>=3 )
			printf ( "updateExtensionPointers: determine extensions that can be used at the next stage\n" );

		std::vector<int> pointers;
		//
		for ( size_t i=0; i<lmctype.size(); i++ ) {
			// we need proper strength
			if ( strengthcheck[i] ) {
				//printf("##  depth_extend_sub_t.updateExtensionPointers: i %zu, lastcol %d extcol %d \n", i, lastcol[i], extcol);
				if ( lastcol[i]>=extcol || lastcol[i]==-1 || extcol<5 ) {  // NOTE: extcol < 5 condition --> make generic
					// good candidate
					//      printf("  %d %d \n", lastcol[i], extcol);
					pointers.push_back ( valididx[i] );
				} else {
					// printf("##  %d %d \n", lastcol[i], extcol);

				}
			}
		}
		if ( verbose>=2 )
			printf ( "updateExtensionPointers: extcol %d, kept %ld/%ld pointers\n", extcol, pointers.size(), lmctype.size() );
		return pointers;
	}

/// initialize the new list of extension columns
	arraylist_t initialize ( const arraylist_t& alist, const arraydata_t &adf, const OAextend &oaextend );

/// select the arrays with are LMC and hence need to be written to disk
	inline  arraylist_t selectArraysZ ( const arraylist_t &alist ) const {
		if ( verbose>=2 )
			printf ( "depth_extend_sub_t: selectArrays: alist.size() %ld, lmctype %ld\n", alist.size(), lmctype.size() );
		arraylist_t ga;
		for ( size_t i=0; i<lmctype.size(); i++ ) {
			//size_t ii = valididx[i];
			if ( verbose>=3 )
				printf ( "  depth_extend_sub_t.selectArraysZ: array %ld: lmctype %d\n", i, lmctype[i] );
			if ( lmctype[i]>=LMC_EQUAL ) {
				//if (i>valididx.size())
				//  printf("  error: i %zu size %zu\n", i, valididx.size() );
				//if (valididx[i]>(int) alist.size())
				//  printf("  error: validx[i] %d alist size %zu\n",  valididx[i], alist.size());
				array_link ee = alist[i];

				ga.push_back ( ee );
				// printf ( "  selectArraysZ: selected array %zu/%zu\n", i, lmctype.size() );
			}
		}
		if ( verbose )
			printf ( "dextend_sub_t: selected %d/%d arrays\n", ( int ) ga.size(), ( int ) alist.size() );
		return ga;
	}

	inline  arraylist_t selectArraysXX ( const array_link &al, const arraylist_t &elist ) const {
		if ( verbose>=2 )
			printf ( "depth_extend_sub_t: selectArraysXX: alist.size() %ld, lmctype %ld\n", elist.size(), lmctype.size() );
		arraylist_t ga;
		for ( size_t i=0; i<n(); i++ ) {
			if ( verbose>=3 )
				printf ( "  selectArraysXX lmctype %d\n", lmctype[i] );
			if ( lmctype[i]>=LMC_EQUAL ) {
				array_link ee = hstacklastcol ( al, elist[valididx[i]] );

				ga.push_back ( ee );
			}
		}
		if ( verbose>=1 )
			printf ( "dextend_sub_t: selected %d/%d arrays\n", ( int ) ga.size(), ( int ) elist.size() );
		return ga;
	}

	inline  void info() const {
		size_t nl=0;
		for ( size_t t=0; t<lmctype.size(); t++ )
			nl += ( lmctype[t]>+LMC_EQUAL );
		//int ngood =std::accumulate(lmctype.begin(),lmctype.end(),0);//#include <numeric>
		if ( verbose ) {
			printf ( "lmc %ld/%d\n", nl, ( int ) lmctype.size() );
			printf ( "valididx size %ld\n", valididx.size() );
		}
	}
};

/** @brief Helper structure for dynamic extension
 *
 * This structure allows for writing the generated arrays to disk.
 * It also contains functions to print progress of the extension.
 * 
 * Multiple copies of this class are made, but they all share the same counter_t and arraywriter_t object.
 *
 */
struct depth_extend_t {
public:
	int verbose;
	OAextend oaextend;
	const arraydata_t *ad;

	int loglevelcol;
	/// print progress every x seconds
	double logtime;

	arraylist_t extension_column_list;	// list of possible extensions

	// shared by mutiple instances of dextend_t (could be made static)
	arraywriter_t *arraywriter;
	counter_t *counter;

	/// if set to true write arrays to disk
	int writearrays;

	int discardJ5;  	/// if true, then we discard the designs which have J5 maximal
	long discardJ5number;

	static double t0;	// time since start of calculation
	static double tp;	// time since last progress report

private:
	long narraysmax;
	depth_path_t searchpath;


public:

	// constructure function
	depth_extend_t ( const arraydata_t *ad_ , double _logtime=10000000, int _discardJ5 = -1 ) : verbose ( 1 ), ad ( ad_ ), discardJ5 ( _discardJ5 )  {
		loglevelcol=-1;
		t0=get_time_ms();
		tp=get_time_ms();

		logtime=_logtime;
		if ( ad==0 )
			printf ( "depth_extend_t: pointer to arraydata_t is zero!" );
		//nfound.resize ( ad->N );
		//writefiles=0;
		writearrays=1;
		narraysmax=LONG_MAX-1;

		arraywriter=0;
		counter=0;

		discardJ5number=0;

		searchpath.init ( ad->ncols );
	};

	depth_extend_t ( const depth_extend_t &de ) {
		//printf("depth_extend_t: copy constructor\n"); printf(" searchpath: "); de.searchpath.show(16);
		verbose=de.verbose;
		oaextend=de.oaextend;
		ad = de.ad;
		loglevelcol=de.loglevelcol;
		logtime=de.logtime;
		//cache_extensions=de.cache_extensions;
		extension_column_list=de.extension_column_list;
		arraywriter=de.arraywriter;
		writearrays=de.writearrays;
		narraysmax=de.narraysmax;
		searchpath=de.searchpath;
		discardJ5=de.discardJ5;
		discardJ5number=de.discardJ5number;

		counter=de.counter;
	}
	~depth_extend_t() {
		//closeafiles();
	}


public:

	void show() {
		printf ( "depth_extend_t: logtime %.1f [s]\n", logtime );
	}

	void setNarraysMax ( long n ) {
		this->narraysmax=n;
	}

	// helper function, thread safe
	void maxArrayCheck() {
		if ( arraywriter==0 )
			return;
		#pragma omp critical
		{
			if ( arraywriter->nwritten>this->narraysmax ) { /// HACK
				printf ( "dextend_t: number of arrays written: %d, quitting\n", arraywriter->nwritten );
				this->counter->showcounts ( *this->ad );
				this->arraywriter->closeafiles();
				//printfd ( "symmetry time: %.1f [s]\n", getdtsymm() );
				//printfd ( "lmc time: %.3f [s]\n", getdtlmc() );

				exit ( 0 );
			}
		}
	}


public:

	void showsearchpath ( int depth ) const {
		searchpath.show ( depth );
	}

	/// show information about the progress of the loop
	bool showprogress ( int showtime=1, int depth= 0, int forcelog=0 )  {
		{
			double dt0 = get_time_ms()-t0;
			double dt = get_time_ms()-tp;
			if ( ( dt>logtime ) || forcelog ) {
				#pragma omp critical
				tp=get_time_ms();
				this->arraywriter->flush();
				if ( showtime ) {
					int na= this->counter->nArrays();

#ifdef DOOPENMP
					printf ( "-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s), thread %d/%d\n", dt0, na, na/dt0, omp_get_thread_num(), omp_get_num_threads() );
#else
					printf ( "-- depth_extend: progress: %.1f [s], narrays %d (%.1f arrays/s)\n", dt0, na, na/dt0 );
#endif

					if ( depth>0 ) {
						printf ( "-- depth %2d: %.1f [s]: ", depth, dt0 );
						searchpath.show ( depth );
					}
				}
				//
				return true;
				//info();
			} else
				return false;
		}
	}
	inline  void info() const {
		printf ( "depth_extend: " );
		ad->show();
	}


	/// set the position in the dextend structure
	void setposition ( int k, int c, int m, int extensioncols=-1, int goodextensioncols=-1 ) {
		#pragma omp critical
		{
			searchpath.updatePosition ( k, c, m, extensioncols, goodextensioncols );
		}
	}
	/// set the position in the dextend structure
	void setpositionGEC ( int k, int goodextensioncols ) {
		#pragma omp critical
		{
			searchpath.updatePositionGEC ( k, goodextensioncols );
		}
	}
};

enum depth_alg_t {DEPTH_DIRECT, DEPTH_EXTENSIONS} ;


/// Helper structure for the even-odd depth extension
struct depth_extensions_storage_t {

	void resize ( size_t s ) {
		columnextensionsList.resize ( s );
		goodarrayslist.resize ( s );
		depthalglist.resize ( s );
		dextendsubList.resize ( s );
	}

	void set ( int ai, const arraylist_t &goodarrays, const arraylist_t &extension_column_list, depth_alg_t depthalg, const depth_extend_sub_t &dextendsub ) {
		this->goodarrayslist[ai]= ( goodarrays );
		this->columnextensionsList[ai]= ( extension_column_list );
		this->depthalglist[ai]= ( depthalg );
		this->dextendsubList[ai]= ( dextendsub );

	}
	//std::vector<arraylist_t>  extensions0list;
	std::vector<arraylist_t>  columnextensionsList;
	std::vector<arraylist_t>  goodarrayslist;
	std::vector<depth_alg_t> depthalglist;
	std::vector<depth_extend_sub_t> dextendsubList;
};



// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
