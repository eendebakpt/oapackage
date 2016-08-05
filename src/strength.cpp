/*! \file strength.cpp
    \brief Contains code to perform a strength check on an OA

   Author: Pieter Eendebak <pieter.eendebak@gmail.com>
   Copyright: See LICENSE.txt file that comes with this distribution
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

#include "mathtools.h"
#include "strength.h"
#include "arraytools.h"


/**
 * @brief Constructor fuction
 */
extend_data_t::extend_data_t ( const arraydata_t *ad, colindex_t extcol ) : adata ( ad ), extcolumn ( extcol ), N ( ad->N ),  gidx ( 0 ), gstart ( 0 ), gsize ( 0 ), freqtable ( 0 ), freqtable_cache ( 0 ), range_low ( -1 ), range_high ( -1 )
{
#ifdef COUNTELEMENTCHECK
    this->elements=0;
#endif
    //OPTIMIZE: make some variables constant
    log_print ( DEBUG, "extend_data_t constructor: strength %d\n", adata->strength );

    extend_data_t *es = this;

    oaindextmin = get_oaindex ( es->adata->s, es->adata->strength-1, es->adata->N );

    /* set column combinations with extending column fixed */
    int ncolcombs;
    es->colcombs = set_colcombs_fixed ( es->lambda, es->nvalues, ncolcombs, ad->s, ad->strength, extcol, ad->N );
    es->lambda2lvl = es->adata->N/pow ( double(2), ad->strength ); // NOTE: use cross-platform interger power
    //printf("lambda2lvl: %d %d\n", es->lambda[0], es->lambda2lvl);
    es->indices = set_indices ( es->colcombs, ad->s, ad->strength, ncolcombs );	//sets indices for frequencies, does the malloc as well
    es->r_index = create_reverse_colcombs_fixed ( ncolcombs );
    es->ncolcombs = ncolcombs;

    // allocate table for frequency count
    this->freqtable = new_strength_freq_table ( es->ncolcombs, es->nvalues, this->freqtablesize );
    log_print ( DEBUG, "Allocated freq table of size %d\n", this->freqtablesize );

    // allocate table for frequency count cache system
    this->freqtable_cache = new strength_freq_table [this->N];
    for ( int zz=0; zz<this->N; zz++ ) {
        int dummy;
        this->freqtable_cache[zz] = new_strength_freq_table ( es->ncolcombs, es->nvalues, dummy );
    }

#ifdef FREQELEM
    // allocate table for frequency count cache system
    int x = this->N*ad->s[extcol];

    this->freqtable_elem = malloc2d<int> ( x, this->ncolcombs );
    this->element2freqtable = malloc2d<int> ( x, this->ncolcombs );
#endif

}


/**
 * @brief Destructor function
 */
extend_data_t::~extend_data_t()
{
    extend_data_t *es = this;
    free_colcombs_fixed ( es->colcombs, es->lambda, es->nvalues );
    free ( es->r_index->index );
    free ( es->r_index );

    free2d ( es->indices ); //, es->ncolcombs);

#ifdef COUNTELEMENTCHECK
    // OPTIMIZE: better to check these pointers are non-zero
    free ( this->elements );
#endif

    free2d_irr ( this->freqtable );

    // clear frequencies cache
    for ( rowindex_t zz=0; zz<this->N; zz++ )
        free2d_irr ( this->freqtable_cache[zz] );

    safedeletearray ( this->freqtable_cache );

    delete [] this->gidx;
    delete [] this->gstart;
    delete [] this->gsize;


#ifdef FREQELEM
    free2d ( this->freqtable_elem );
    free2d ( this->element2freqtable );
#endif

}


strength_freq_table new_strength_freq_table ( int ncolcombs, int *nvalues, int &nelements )
{
    strength_freq_table frequencies = malloc2d_irr<int> ( ncolcombs, nvalues, nelements );
    return frequencies;
}


/**
 * @brief Destructor
 * @param colcombs
 * @param lambda
 * @param nvalues
 */
void free_colcombs_fixed ( colindex_t **colcombs, int *lambda, int *nvalues )
{
    free2d<colindex_t> ( colcombs );
    free ( lambda );
    free ( nvalues );
}


/**
 * @brief Creat a reverse column combination structure to all columns
 * @param ncolcombs
 * @return
 */
rev_index *create_reverse_colcombs_fixed ( const int ncolcombs )
{
    register int	i, j;
    rev_index	*rev_colcombs;

    log_print ( DEBUG, "Creating reverse column combination index: ncolcombs %d\n", ncolcombs );

    rev_colcombs = ( rev_index* ) malloc ( sizeof ( rev_index ) );

    rev_colcombs[0].nr_elements = ncolcombs;
    rev_colcombs[0].index = ( int* ) malloc ( ncolcombs * sizeof ( int ) );
    for ( i=0; i<ncolcombs; i++ )
        rev_colcombs->index[i]=i;

    // print the resulting reverse index
    for ( j = 0; j < rev_colcombs->nr_elements; j++ )
        log_print ( DEBUG, "%i, ", rev_colcombs->index[j] );
    log_print ( DEBUG, "\n" );

    log_print ( DEBUG, "   create_reverse_colcombs: Done\n" );

    return rev_colcombs;
}


/**
 * Only include the column combinations with the final column involved
 *
 * @brief Create lookup table of value to column combination
 * @param colcombs
 * @param ncols
 * @param strength
 * @return
 */
rev_index *create_reverse_colcombs ( colindex_t **colcombs, const int ncols, const int strength )
{
    int		count, *tmp;
    rev_index	*rev_colcombs;
    const int	ncolcombs = ncombs ( ncols, strength );

    rev_colcombs = ( rev_index* ) malloc ( ncols * sizeof ( rev_index ) );
    tmp = new int [ncombs ( ncols, strength )];

    log_print ( DEBUG, "ncols = %i\n", ncols );

    for ( int k = 0; k < ncols; k++ ) {
        count = 0;
        for ( int i = 0; i < ncolcombs; i++ ) {
            for ( int j = 0; j < strength; j++ ) {
                if ( colcombs[i][j] == k ) {
                    tmp[count++] = i;
                    break;
                }
            }
        }
        rev_colcombs[k].nr_elements = count;
        log_print ( DEBUG, "For column %i found %i combis\n", k, count );
        rev_colcombs[k].index = ( int* ) malloc ( count * sizeof ( int ) );
        memcpy ( rev_colcombs[k].index, tmp, count * sizeof ( int ) );
    }
    delete [] tmp;

    // print the resulting reverse index
    if ( log_print ( DEBUG,"" ) ) {
        for ( int i = 0; i < ncols; i++ ) {
            log_print ( DEBUG, "%i:\t", i );
            for ( int j = 0; j < rev_colcombs[i].nr_elements; j++ )
                log_print ( DEBUG, "%i, ", rev_colcombs[i].index[j] );
            log_print ( DEBUG, "\n" );
        }
    }
    return rev_colcombs;
}


/**
 * @brief Calculate position of a row in the frequency table
 * @param N
 * @param activerow
 * @param ncols
 * @param colcomb
 * @param indices
 * @param array
 * @return
 */
inline int freq_position ( rowindex_t N, rowindex_t activerow, colindex_t ncols, colindex_t *colcomb, int *indices, carray_t *array )
{
    int freq_pos = 0;
    int oa_pos;
    for ( colindex_t j = 0; j < ncols; j++ ) { //use all columns involved in combination
        oa_pos = N * colcomb[j] + activerow;	//calculate position in OA for column involved
        freq_pos += array[oa_pos] * indices[j];	//calculate position in freq vector for combination involved
    }
    return freq_pos;
}

/* NOTE:
 *
 * The size of the freqtable is pretty larg compared to the number of nonzero elements
 * for a single row-column-value pair. Hence we use a sparse representation of this table where possible.
 *
 */
#ifdef FREQELEM


/** @brief Add row to frequency table using cache system
 *  Add row to frequency table using cache system. The cache systems keeps track for every row-value combination
 * for the active column what the corresponding t-tuples of values are for all relevant column combinations.
 *
 * @param es
 * @param activerow
 * @param freqtable
 */
void add_element_freqtable ( extend_data_t *es, rowindex_t activerow, carray_t *array, strength_freq_table freqtable )
{
    const array_t elem = array[es->extcolumn*es->N + activerow];
    int idx = es->adata->s[es->extcolumn]*activerow+elem;

    //int *postable = es->freqtable_elem[idx];
    int *postable3 = es->element2freqtable[idx];
    freq_t* freqtable0 = & ( freqtable[0][0] );
    for ( int z=0; z<es->ncolcombs; z++ ) {
        //int freqpos = postable[z]; freqtable[z][freqpos]++;

        //printf("%ld -> %ld\n", freqtable[z]+freqpos, postable2[z] );
        freqtable0[postable3[z]]++;
        //printf("row %d elem %d, z %d, freqpos %d\n", activerow, elem, z, freqpos);
    }
}

void add_element_freqtable_col ( extend_data_t *es, rowindex_t activerow, carray_t *arraycol, strength_freq_table freqtable )
{
    const array_t elem = arraycol[activerow];
    int idx = es->adata->s[es->extcolumn]*activerow+elem;

    //int *postable = es->freqtable_elem[idx];
    int *postable3 = es->element2freqtable[idx];
    freq_t* freqtable0 = & ( freqtable[0][0] );
    for ( int z=0; z<es->ncolcombs; ++z ) {
        //int freqpos = postable[z]; freqtable[z][freqpos]++;

        //printf("%ld -> %ld\n", freqtable[z]+freqpos, postable2[z] );
        freqtable0[postable3[z]]++;
        //printf("row %d elem %d, z %d, freqpos %d\n", activerow, elem, z, freqpos);
    }
}

/*
void add_element_freqtable_legacy(rowindex_t N, extend_data_t *es, rowindex_t activerow, carray_t *array, strength_freq_table freqtable)
{
	const array_t elem = array[es->extcolumn*es->N + activerow];
	int idx = es->adata->s[es->extcolumn]*activerow+elem;

	int *postable = es->freqtable_elem[idx];
	for(int z=0;z<es->ncolcombs;z++) {
		int freqpos = postable[z];
		freqtable[z][freqpos]++;
	}
}
*/
#else

/* *
 * Add a new row to the frequency count table.
 *
 * @param es
 * @param activerow
 * @param array
 * @param freqtable
 */
void add_element_freqtable ( extend_data_t *es, rowindex_t activerow, carray_t *array, strength_freq_table freqtable )
{
    log_print ( NORMAL, "add_element: es->r_index->nr_elements %d\n", es->r_index->nr_elements );

    int		freq_pos = 0;			//is position in frequency array
    int cur_combi = 0;		//is position in OA

    colindex_t strength = es->adata->strength;
    rowindex_t N = es->adata->N;


    for ( int i = 0; i < es->r_index->nr_elements; i++ ) { //find all t-tuples of columns with column p->col included
        cur_combi = es->r_index->index[i];	/* select a combination of t columns */
        freq_pos = freq_position ( N, activerow, strength, es->colcombs[cur_combi],  es->indices[cur_combi], array );

#ifdef OADEBUG
        if ( freq_pos < 0 ) {
            myprintf ( "add_element: this code should not be executed\n" );
        }
#endif
        ( freqtable[cur_combi][freq_pos] ) ++;	//add occurence of combination
    }
}

#endif

/**
 * @brief Print frequencies structure
 * @param frequencies
 * @param nelements
 * @param lambda
 * @param N
 */
void print_frequencies ( int **frequencies, const int nelements, const int *lambda, const int N )
{
    register int	i, j;

    for ( i = 0; i < nelements; i++ ) {
        myprintf ( "%i:\t", i );
        for ( j = 0; j < N/lambda[i]; j++ ) {
            myprintf ( "%2i ", frequencies[i][j] );
        }
        myprintf ( "\n" );
    }
    myprintf ( "\n" );
}

/**
 * @brief Initialize the table of t-tuple frequencies
 * @param es
 * @param array
 */
void init_frequencies ( extend_data_t *es, array_t *array )
{
    log_print ( DEBUG, "init_frequencies: extension column %d\n", es->extcolumn );
    array_t sc = es->adata->s[es->extcolumn];

    /* loop over rows */
    for ( rowindex_t row=0; row<es->adata->N; row++ ) {
        /* loop over all values of the rows */
        for ( array_t v=0; v<sc; v++ ) {
            int idx =row*sc+v;

            //log_print(NORMAL-1, "init_frequencies: setting idx %d\n", idx);
            /* clear data */
            memset ( es->freqtable_elem[idx], 0, es->ncolcombs * sizeof ( int ) );
            memset ( es->element2freqtable[idx], 0, es->ncolcombs * sizeof ( int ) );
            /* count */
            int cur_combi = 0;		//is position in OA

            colindex_t strength = es->adata->strength;
            rowindex_t N = es->adata->N;

            for ( int i = 0; i < es->ncolcombs; i++ ) { // find all columns column p->col is involved with
                cur_combi = es->r_index->index[i];	/* select a combination of t columns */
                array_t tmp = array[es->adata->N*es->extcolumn+row];
                array[es->adata->N*es->extcolumn+row]=v;
                int freq_pos = freq_position ( N, row, strength, es->colcombs[cur_combi],  es->indices[cur_combi], array );
                es->freqtable_elem[idx][i] = freq_pos;
                es->element2freqtable[idx][i] = & ( es->freqtable[i][freq_pos] )-& ( es->freqtable[0][0] );
                array[es->adata->N*es->extcolumn+row] = tmp;
            }

        }
    }
}

/** @brief Reset the table of frequencies of t-tuples and fill it for a selection of rows
 *
 * @param frequencies
 * @param es
 * @param p
 * @param currentcol
 * @param rowstart
 * @param rowlast
 * @param array
 */
void recount_frequencies ( int **frequencies, extend_data_t *es, colindex_t currentcol, rowindex_t rowstart, rowindex_t rowlast, carray_t *array )
{
    //myprintf("recount_frequencies: rowstart %d, rowlast %d, es->ncolcombs %d\n", rowstart, rowlast, es->ncolcombs);
    //reset old values
    for ( int i = 0; i < es->ncolcombs; i++ ) {
        memset ( frequencies[i], 0, es->nvalues[i] * sizeof ( int ) );
    }

    //recount new values
    for ( int i = rowstart; i <= rowlast; i++ ) {
        add_element_freqtable ( es, i, array, es->freqtable );
    }
}

#ifdef FULLPACKAGE

bool valid_element ( const extend_data_t *es, const extendpos *p, carray_t *array )
/*Put the possibliy correct value in p->value and it will be tested, given the frequency table, lambdas etc*/
{

#ifdef FREQELEM
    /* perform strength check using frequency element cache */
    const int idx =p->row*es->adata->s[es->extcolumn] + p->value;
    //int *freqpositions = es->freqtable_elem[idx];
    int *freqpositions2 = es->element2freqtable[idx];
    freq_t* freqtable0 = & ( es->freqtable[0][0] );

    // OPTIMIZE: we can improve the performance by changing the ordering in the frequency table. this makes a difference for cases with high strength and high number of columns
    for ( int z=0; z<es->ncolcombs; z++ ) {
        if ( 0 ) {
            //int freq_pos = freqpositions[z];
            //if(es->freqtable[z][freq_pos] + 1 > es->lambda[z])
            //	return false;
        } else {
            ;
            if ( freqtable0[freqpositions2[z]] + 1 > es->lambda[z] )
                return false;

        }
    }
    return true;
#else
    int		freq_pos;		//position in frequency array
    int	cur_combi;	//position in OA
    const array_t	tmp_oa = array[p->col * p->ad->N + p->row];
    array_t *acolp = &array[p->ad->N * p->col];

    acolp[p->row] = p->value;
    //for(int i = es->r_index->nr_elements-1; i >=0; i--)	//loop over all column combinations
    for ( int i = 0; i < es->r_index->nr_elements; i++ ) {	//loop over all column combinations
        cur_combi = es->r_index->index[i];
        // freq_pos is the value for this column combination
        freq_pos = freq_position ( p->ad->N, p->row, p->ad->strength, es->colcombs[cur_combi],  es->indices[cur_combi], array );

        if ( es->freqtable[cur_combi][freq_pos]  >= es->lambda[cur_combi] ) {
            acolp[p->row] = tmp_oa;
            return false;
        }
    }
    acolp[p->row] = tmp_oa;
    return true;
#endif
}

bool valid_element_2level ( const extend_data_t *es, const extendpos *p )
/*Put the possibliy correct value in p->value and it will be tested, given the frequency table, lambdas etc*/
{

#ifdef FREQELEM
    /* perform strength check using frequency element cache */
    const int idx =p->row*es->adata->s[es->extcolumn] + p->value;
    //int *freqpositions = es->freqtable_elem[idx];
    int *freqpositions2 = es->element2freqtable[idx];
    freq_t* freqtable0 = & ( es->freqtable[0][0] );

    // OPTIMIZE: we can improve the performance by changing the ordering in the frequency table. this makes a difference for cases with high strength and high number of columns
    for ( int z=0; z<es->ncolcombs; z++ ) {
        //int freq_pos = freqpositions[z];

        //OPTIMIZE: create a lambda minus one table
        //IDEA: use structure of freqtable to make better indexing

        if ( freqtable0[freqpositions2[z]] >= es->lambda2lvl ) { // (es->freqtable[z][freq_pos]  >= es->lambda2lvl)
            //cout << "colcomb:" << z << endl;
            return false;
        }
    }
    return true;
#else
    // not implemented
    printf ( "valid_element_2level: not implemented\n" );
#endif
}

#endif

bool strength_check ( const arraydata_t &ad, const array_link &al,  int verbose )
{
#ifdef OADEBUG
    myassert ( ad.ncols>=al.n_columns, "strength_check: array has too many columns" );
#endif

    if (ad.strength==0)
        return true;

    strength_check_t strengthcheck ( ad.strength );

    /* set column combinations with extending column fixed */
    int fixcol=al.n_columns-1;
    if ( verbose>=2 )
        myprintf ( "strength_check array: N %d, k %d, strength %d\n", ad.N, al.n_columns, ad.strength );
    strengthcheck.colcombs = set_colcombs_fixed ( strengthcheck.lambda, strengthcheck.nvalues, strengthcheck.ncolcombs, ad.s, ad.strength, fixcol, ad.N );

    //myprintf ( "nvalues: " ); print_perm ( nvalues, strengthcheck.ncolcombs );
    strengthcheck.indices = set_indices ( strengthcheck.colcombs, ad.s, ad.strength, strengthcheck.ncolcombs );	//sets indices for frequencies, does the malloc as well
    strengthcheck.create_reverse_colcombs_fixed();

//	strengthcheck.r_index = create_reverse_colcombs_fixed(ncolcombs);


    int val=true;
    strengthcheck.freqtable = new_strength_freq_table ( strengthcheck.ncolcombs, strengthcheck.nvalues, strengthcheck.freqtablesize );

    if ( verbose>=2 )
        strengthcheck.info();

    if ( verbose>=2 ) {
        myprintf ( "before:\n" );
        strengthcheck.print_frequencies ( );
    }
//   myprintf ( "  table of size %d\n", strengthcheck.freqtablesize );

    for ( int i=0; i<strengthcheck.ncolcombs; i++ ) {
        //myprintf ( "columns %d: ", i ); print_perm ( strengthcheck.colcombs[i], strength );

        assert ( ad.N<=MAXROWS );
        int valindex[MAXROWS];
        std::fill_n(valindex, ad.N, 0);
        for ( int t=0; t<ad.strength; t++ ) {
            colindex_t cc = strengthcheck.colcombs[i][t];
            int s = ad.s[cc];
            array_t *array_coloffset=al.array+cc*ad.N;
            for ( int r=0; r<ad.N; r++ ) {
                array_t val = array_coloffset[r];
                valindex[r] = valindex[r]*s+val;
            }
        }
        for ( int r=0; r<ad.N; r++ ) {
            int vi = valindex[r];
            strengthcheck.freqtable[i][vi]++;
        }

        for ( int j=0; j<strengthcheck.nvalues[i]; j++ ) {
            //    myprintf ( "strength: i %d, j %d: %d %d\n", i, j, strengthcheck.freqtable[i][j], nvalues[i] );
            if ( strengthcheck.freqtable[i][j]!=ad.N/strengthcheck.nvalues[i] ) {
                if ( verbose>=2 )
                    myprintf ( "no good strength: i %d, j %d: %d %d\n", i, j, strengthcheck.freqtable[i][j], strengthcheck.nvalues[i] );
                val=false;
                break;
            }
        }

        if ( val==false )
            break;
    }
    //myprintf ( "nvalues: " ); print_perm ( nvalues, strengthcheck.ncolcombs );
    if ( verbose>=2 ) {
        myprintf ( "table of counted value pairs\n" );
        strengthcheck.print_frequencies ( );
    }
//	print_frequencies(strengthcheck.freqtable, strengthcheck.ncolcombs, nvalues, ad.N);

    return val;
}



/*!
  check_divisibility checks if the number of runs is a multiple of any combination of the number of factors.
  The function only returns true or false, further alction should be done by the calling function.
  \brief Checks on the divisibility of the number of runs by the product of the levels in the factors for all t-tuples
  */
bool check_divisibility ( const arraydata_t *ad )
{
    const int	ncolcombs = ncombs ( ad->ncols, ad->strength );
    int		prod, *colcombs = 0;
    bool		ret = true;

    colcombs = ( int* ) malloc ( ad->strength * sizeof ( int ) );

    /* loop over all combinations of t-tuples of columns */
    for ( int i = 0; i < ncolcombs; i++ ) {
        if ( i == 0 )						//first combination needs to be set
            for ( int j = 0; j < ad->strength; j++ )
                colcombs[j] = j;
        else
            next_comb ( colcombs, ad->strength, ad->ncols );	//get next combination

        prod = 1;
        for ( int j = 0; j < ad->strength; j++ ) {
            prod *= ad->s[colcombs[j]];
            //cout << "prod: " << prod << endl;
        }
        if ( ad->N % prod != 0 ) {
            log_print ( SYSTEM, "Failed divisibility test!\n" );
            log_print ( SYSTEM, "Column combination: " );
            print_perm ( colcombs, ad->strength );
            log_print ( SYSTEM, "N %d, product of s[j] is %d and %i %% %i != %i\n", ad->N, prod, ad->N, prod, ad->N % prod );
            fflush ( NULL );
            ret = false;
            break;
        }
    }
    free ( colcombs );

    return ret;
}


colindex_t **set_colcombs_fixed ( int*& xlambda, int*& nvalues, int &ncolcombs, const array_t *s, const int strength, const int fixedcol, const int N )
{
    log_print ( DEBUG+1, "set_colcombs_fixed: strength %d, fixedcol %d, N%d\n", strength, fixedcol, N );
    register int	i,j;
    int		prod;
    colindex_t **colcombs = 0;

    int n = fixedcol;
    int k = strength - 1;	// we keep 1 column fixed, choose k columns
    ncolcombs = ncombs ( n, k );

    //printf("ncolcombs: %d, strength %d\n", ncolcombs, strength);
    colcombs = malloc2d<colindex_t> ( ncolcombs, strength );
    xlambda = ( int* ) malloc ( ncolcombs * sizeof ( int ) );
    nvalues = ( int* ) malloc ( ncolcombs * sizeof ( int ) );

    log_print ( DEBUG, "ncolcombs: %d\n", ncolcombs );

    /* set last entry to fixed column */
    for ( i=0; i< ncolcombs; i++ )
        colcombs[i][k] = fixedcol;

    for ( i = 0; i < k; i++ )	//set initial combination
        colcombs[0][i] = i;


    for ( i = 1; i < ncolcombs; i++ ) {
        memcpy ( colcombs[i], colcombs[i-1], k * sizeof ( colindex_t ) );
        next_combination<colindex_t> ( colcombs[i], k, n );

        prod = 1;
        for ( j = 0; j < strength; j++ ) {
            //myprintf("i %d j %d: %d\n", i, j, colcombs[i][j]);
            prod *= s[colcombs[i][j]];
        }
        nvalues[i] = prod;
        xlambda[i] = N/prod;
    }

    prod = 1;	//First lambda manually, because of copy-for-loop
    for ( i = 0; i < strength; i++ )
        prod *= s[colcombs[0][i]];
    nvalues[0] = prod;
    xlambda[0] = N/prod;

    return colcombs;
}


// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
