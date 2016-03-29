/*! \file evenodd.h
 *  \brief Contains functions to generate even-odd designs
 *
 */

#pragma once

#include "arraytools.h"
//#include "arrayproperties.h"

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


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
