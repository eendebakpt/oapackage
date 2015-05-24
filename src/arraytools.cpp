#include <string>
#include <sstream>


#include "arraytools.h"
#include "bitarray/bit_array.h"

#include "tools.h"
#include "mathtools.h"
#include "arrayproperties.h"
#include "strength.h"

#include "lmc.h"

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;

#ifdef WIN32
#else
#include <fcntl.h>
int get_file_status ( FILE* f )
{
	int fd = fileno ( f );
	return fcntl ( fd, F_GETFL );

}

#endif

std::vector<int> array_transformation_t::rowperm() const
{
	std::vector<int> ww ( this->rperm, this->rperm+this->ad->N );
	return ww;
}

std::vector<int> array_transformation_t::colperm() const
{
	std::vector<int> ww ( this->cperm, this->cperm+this->ad->ncols);
	return ww;
}
std::vector<int> array_transformation_t::lvlperm ( int c ) const
{
	if ( c<0 || c>= this->ad->ncols ) {
		std::vector<int> ww;
		return ww;
	}
	std::vector<int> ww ( this->lperms[c], this->lperms[c]+this->ad->s[c] );
	return ww;
}


void array_transformation_t::setrowperm ( std::vector<int>rp )
{
	if ( ( int ) rp.size() !=this->ad->N ) {
		printf ( "array_transformation_t::setrowperm: argument has wrong dimensions\n" );
		return;
	}
	std::copy ( rp.begin(), rp.end(), this->rperm );
}
void array_transformation_t::setcolperm ( std::vector<int>colperm )
{
	if ( ( int ) colperm.size() !=this->ad->ncols ) {
		printf ( "array_transformation_t::setrowperm: argument has wrong dimensions\n" );
		return;
	}
	std::copy ( colperm.begin(), colperm.end(), this->cperm );
}
void array_transformation_t::setlevelperm ( int colindex, std::vector<int> lvlperm )
{
	if ( colindex<0 || colindex>=this->ad->ncols ) {
		printf ( "array_transformation_t::setrowperm: argument has wrong dimensions\n" );
		return;
	}
	if ( ( int ) lvlperm.size() !=this->ad->s[colindex] ) {
		printf ( "array_transformation_t::setrowperm: argument has wrong dimensions\n" );
		return;
	}
	std::copy ( lvlperm.begin(), lvlperm.end(), this->lperms[colindex] );
}




/**
 * @brief Print an array transformation to an output stream
 * @param out
 */
void array_transformation_t::show ( std::ostream & out ) const
{
	out << "array transformation: N " << ad->N;
	if ( this->isIdentity() ) {
		out << ": identity transformation" << std::endl;
	} else {
		out << std::endl;
		out << "column permutation: ";
		print_perm ( out, cperm, ad->ncols );

		out << "level perms:" << endl;
		for ( colindex_t c=0; c<ad->ncols; c++ ) {
			print_perm ( out, lperms[c], ad->s[c] );
		}

		out << "row permutation: ";
		print_perm ( out, rperm, ad->N );
	}
}

void array_transformation_t::show ( ) const
{
	this->show ( cout );
}

array_transformation_t::array_transformation_t ( )
{
	ad = 0;
	rperm = 0;
	lperms = 0;
	cperm= 0;
}

/**
 * @brief Create an OA transformation
 * @param adp
 */
array_transformation_t::array_transformation_t ( const arraydata_t *adp )
{
	//printf("array_transformation_t::array_transformation_t: constructor with arraydata_t pointer\n");
	ad = new arraydata_t ( *adp );
	init();
}

/// Copy construction
array_transformation_t::array_transformation_t ( const array_transformation_t &tt )
{

	//printf("array_transformation_t::array_transformation_t: constructor with array_transformation_t\n");
	ad = new arraydata_t ( * ( tt.ad ) );

	init();

	// copy data
	std::copy ( tt.rperm, tt.rperm+ad->N, rperm );
	std::copy ( tt.cperm, tt.cperm+ad->ncols, cperm );
	for ( colindex_t c=0; c<ad->ncols; c++ ) {
		std::copy ( tt.lperms[c], tt.lperms[c]+ ad->s[c], lperms[c] );
	}

}

void array_transformation_t::reset()
{
	init_perm ( this->rperm, this->ad->N );
	init_perm<colindex_t> ( this->cperm, ad->ncols );

	for ( colindex_t c=0; c<ad->ncols; c++ )
		init_perm ( lperms[c], ad->s[c] );

}

void array_transformation_t::init()
{
//printf("array_transformation_t::init\n");
	rperm = new_perm_init<rowindex_t> ( ad->N );
	cperm = new_perm_init<colindex_t> ( ad->ncols );

	lperms = new levelperm_t [ad->ncols];
	for ( colindex_t c=0; c<ad->ncols; c++ )
		lperms[c] = new_perm_init<array_t> ( ad->s[c] );
}

void array_transformation_t::free()
{
	if ( ad==0 )
		return;

//printf("array_transformation_t::free\n");
	delete_perm ( rperm );
//printf("array_transformation_t::free: delete colperm (ad %ld)\n", long(ad));
	delete_perm ( cperm );
//printf("array_transformation_t::free: delete lperms\n");
	for ( colindex_t c=0; c<ad->ncols; c++ ) {
		delete_perm ( lperms[c] );
	}
	delete [] lperms;

//printf("array_transformation_t::free: delete ad\n");
	if ( ad!=0 )
		delete ad;
}

/// Assignment operator
array_transformation_t& array_transformation_t::operator= ( const array_transformation_t &tt )
{
	// TODO: check we have transformations of similar type?

// printf("array_transformation_t::operator=: assignment operator\n");

	free();

// printf("array_transformation_t::operator=: init new data\n");

	ad = new arraydata_t ( * ( tt.ad ) );
	init();


	// printf("array_transformation_t::operator=: copy data\n");

	// copy data
	std::copy ( tt.rperm, tt.rperm+ad->N, rperm );
	std::copy ( tt.cperm, tt.cperm+ad->ncols, cperm );
	for ( colindex_t c=0; c<ad->ncols; c++ ) {
		std::copy ( tt.lperms[c], tt.lperms[c]+ ad->s[c], lperms[c] );
	}

	return *this;
}

/**
 * @brief array_transformation_t destructor
 */
array_transformation_t::~array_transformation_t()
{
	//this->print(cout);
	delete_perm ( rperm );
	delete_perm ( cperm );
	for ( colindex_t c=0; c<ad->ncols; c++ ) {
		delete_perm ( lperms[c] );
	}
	delete [] lperms;

	delete ad;
}

bool array_transformation_t::isIdentity() const
{
	//	printf("isIdentity:\n");
	for ( int i=0; i<ad->ncols; ++i ) {
		if ( cperm[i]!=i ) {
			return 0;
		}
	}
	//	printf("isIdentity: cols good\n");
	for ( int i=0; i<ad->N; ++i ) {
		if ( rperm[i]!=i ) {
			return 0;
		}
	}
	//	printf("isIdentity: rows good\n");
	for ( int c=0; c<ad->ncols; ++c ) {
		for ( int i=0; i<ad->s[c]; ++i ) {
			if ( lperms[c][i]!=i ) {
				return 0;
			}
		}
	}
	//	printf("isIdentity: yes\n");
	return 1;
}

/// apply transformation and show resulting array
void array_transformation_t::print_transformed ( carray_t *source ) const
{
	array_t *a = clone_array ( source, this->ad->N, this->ad->ncols );
	print_array ( a, this->ad->N, this->ad->ncols );
	destroy_array ( a );
}


/// apply transformation to an array
void array_transformation_t::apply ( array_t *sourcetarget )
{
	array_t *tmp = create_array ( ad->N, ad->ncols );
	copy_array ( sourcetarget, tmp, ad->N, ad->ncols );
	this->apply ( tmp, sourcetarget );
	destroy_array ( tmp );
}


/// initialize to a random transformation
void array_transformation_t::randomize()
{
	/* row permutation */
	random_perm ( rperm, ad->N );

	/* column permutation */
	for ( int x=0; x<ad->ncolgroups; x++ ) {
		random_perm ( cperm+ad->colgroupindex[x], +ad->colgroupsize[x] );
	}
	//cout << "random permutation: col perm "; print_perm(colperm, ad->ncols);

	/* level permutation */
	for ( colindex_t c=0; c<ad->ncols; c++ ) {
		random_perm ( lperms[c], ad->s[c] );

		//printf("randomize: col %d: \n", c); print_perm(lperms[c], ad->s[c] );
	}
}

/// initialize to a random transformation
void array_transformation_t::randomizecolperm()
{
	/* column permutation */
	for ( int x=0; x<ad->ncolgroups; x++ ) {
		random_perm ( cperm+ad->colgroupindex[x], +ad->colgroupsize[x] );
	}
}

/// initialize to a random transformation
void array_transformation_t::randomizerowperm()
{
	random_perm ( this->rperm, ad->N );
}
array_transformation_t array_transformation_t::inverse() const
{
	array_transformation_t A ( *this );


	invert_permutation ( this->rperm, this->ad->N, A.rperm );
	invert_permutation ( this->cperm, this->ad->ncols, A.cperm );

	/* level permutations */
	for ( colindex_t ci=0; ci<ad->ncols; ci++ ) {
//	  levelperm_t l1 = b.lperms[tmpcolpermai[ci]];
//	  levelperm_t l2 = a.lperms[ci];

		colindex_t cir = this->cperm[ci];
//	  colindex_t cir = A.colperm[ci];

		invert_permutation ( this->lperms[ci], this->ad->s[ci], A.lperms[cir] );

		if ( ci<0 ) {
			printf ( "ci %d: this->lperms[ci] ", ci );
			print_perm ( this->lperms[ci], this->ad->s[cir] );
			printf ( "cir %d: A->lperms[cir] ", cir );
			print_perm ( A.lperms[cir], this->ad->s[ci] );
		}
	}
	return A;
}

/**
 * Return factors of an array. These can be read from the array of the array is in LMC
 * form? Better is to use an oaconfig file...
 *
 */
void get_factors ( array_t *s, carray_t *array, const rowindex_t N, const colindex_t ncols )
{
	int	count = 0;

	array_t max;
	for ( colindex_t i = 0; i < ncols; i++ ) {
		max = 0;
		for ( rowindex_t j = 0; j < N; j++ )
			if ( array[count++] > max )
				max = array[count - 1];
		s[i] = max + 1;
	}
}

arraydata_t arraylink2arraydata ( const array_link &al, int extracols, int strength )
{
	int verbose=0;

	// create arraydatya
	int ncols0=al.n_columns;
	int N=al.n_rows;
	int ncols = ncols0+extracols;
	std::vector<int> s ( ncols );
	int ss=-1;
	for ( int ik=0; ik<ncols0; ik++ ) {
		array_t *xx=std::max_element ( al.array+N*ik, al.array+N* ( ik+1 ) );
		ss=*xx+1;
		s[ik]=ss;
	}
	for ( int ik=ncols0; ik<ncols; ik++ )
		s[ik]=ss; // repeat last factor value
	arraydata_t ad ( s, al.n_rows, strength, ncols );
	if ( verbose ) {
		printf ( "arraylink2arraydata: " );
		ad.show();
	}
	return ad;
}

/// helper function
void foldtest ( jstruct_t &js, const array_link &al, int jj, int verbose )
{
	const rowindex_t N = al.n_rows;
	const colindex_t k = al.n_columns;
	int *pp = new_comb_init<int> ( jj );
	//printf("array:\n"); print_array(al);

	array_t **tmpcol = malloc2d<array_t> ( jj+1, N );
	std::fill ( tmpcol[0], tmpcol[0]+N, 0 );

	int fold;
	int cb=0;
	int nc=ncombs ( al.n_columns, jj );
	for ( int x=0; x<nc; x++ ) {
		if ( verbose>=4 ) {
			printf ( "x %d: ", x );
			cout << printfstring ( "comb: " );
			print_perm ( pp, jj );
		}
		// update invalid columns
		for ( int ii=cb; ii<jj; ii++ ) {
			for ( int r=0; r<N; r++ ) {
				tmpcol[ii+1][r]=tmpcol[ii][r]+al.array[r+pp[ii]*N];
				//printf(    "         %d = %d + %d\n", tmpcol[ii+1][r], tmpcol[ii][r], al.array[r+pp[ii]*N]);
			}
			//printf("  tmp perm %d: ", ii+1); print_perm(tmpcol[ii+1], N);
		}

		// calculate j-value from last column
		int jval=0;
		int tmp=0;
		for ( int r=0; r<N; r++ ) {
			tmp=tmpcol[jj][r];
			tmp %= 2;
			// tmp *=2; tmp--;
			jval += tmp;
		}
		jval = - ( 2*jval - N );

		js.vals[x]=jval;
		fold = next_combination_fold ( pp, jj, k );
		cb=fold;

		if ( verbose>=2 ) {
			printf ( "comb %d: jval %d\n", x, jval );
		}
		//  cout << printfstring("   J value %d\n", jv);
	}

	free2d ( tmpcol, jj+1 );
	delete_comb ( pp );
}

// / write arrayfile to file in test format
//void write_arrayfile_text ( FILE *fid, carray_t *array, const int nrows, const int ncols, int index = 0 );

/** Return number of arrays with j_{2n+1}=0 for n<m */
std::vector<int> getJcounts ( arraylist_t *arraylist, int N, int k, int verbose )
{
	std::vector<int> aa ( k+1 );

	for ( int i=0; i< ( int ) arraylist->size(); i++ ) {
		if ( verbose ) {
			if ( i<100 || i%1000 == 0 )
				std::cout << "## analyzing array " << i << "/" << arraylist->size() << std::endl;
		}
		array_link al = arraylist->at ( i );

		int jl = 5;
		while ( jl<=al.n_columns ) {
			jstruct_t js ( al.n_rows, al.n_columns, jl );
			if ( 0 ) {
				arraylist_t *all = new arraylist_t();
				all->push_back ( al );
				std::vector<jstruct_t> xx = analyseArrays ( *all, verbose, jl );
				js=jstruct_t ( xx[0] );
				if ( verbose>=2 ) {
					printf ( " old: " );
					xx[0].show();
				}
				delete all;
			} else {
				foldtest ( js,al,  jl,verbose );
			}
			if ( !js.allzero() ) {
				if ( verbose>=3 ) {
					printf ( " new: " );
					js.show();
				}

				aa[jl]++;
				if ( verbose>=3 )
					js.show();
				break;
			}
			jl+=2;
			if ( verbose>=2 )
				printf ( "  jl %d\n", jl );
		}
		if ( jl>al.n_columns ) {
			aa[0]++;
		}


	}
	return aa;
}


/**
 * @brief Default copy constructor, no deep copy of the array for efficiency!!
 * @param A
 */
array_link::array_link ( const array_link &rhs )   //: n_rows(rhs.n_rows), n_columns(rhs.n_columns), index(rhs.index), array(rhs.array)
{
#ifdef CLEAN_ARRAY_LINK
	array=0;
	deepcopy ( rhs );
#else
	shallowcopy ( rhs );
#endif
}

array_link::array_link ( Eigen::MatrixXd &m )
{
	this->index=INDEX_DEFAULT;
	this->n_columns=m.cols();
	this->n_rows=m.rows();
	this->array = create_array ( this->n_rows, this->n_columns );
	for ( int i=0; i<this->n_columns*this->n_rows; i++ ) {
		this->array[i] = m ( i );
	}
}

/// create an array by permuting columns
array_link::array_link ( const array_link &rhs, const std::vector<int> &colperm ) : array ( 0 )
{
	init ( rhs.n_rows, rhs.n_columns );
	for ( int c=0; c<rhs.n_columns; c++ ) {
		int cp = colperm[c];
		std::copy ( rhs.array+n_rows*cp, rhs.array+n_rows* ( cp+1 ), this->array+c*n_rows );
	}
}

/**
 * @brief Element access
 */
array_t array_link::at ( const int index ) const
{
	if ( index<0 || index>this->n_rows*this->n_columns-1 )
		return -1;
	else
		return this->array[index];
}

array_t array_link::_at ( const int index ) const
{
	return this->array[index];
}
array_link array_link::clone() const
{

	array_link al ( *this );
	return al;
}


void array_link::setvalue ( int r, int c, int val )
{
	if ( ( r<0 ) || ( r >= this->n_rows ) || ( c<0 ) || ( c>=this->n_columns ) ) {
		fprintf ( stderr,  "array_link error: index out of bounds %d %d (%d %d)!!\n", r, c, this->n_rows, this->n_columns );
		return;
	}


	this->array[r+this->n_rows*c]= val;
}
void array_link::setvalue ( int r, int c, double val )
{
	if ( ( r<0 ) || ( r >= this->n_rows ) || ( c<0 ) || ( c>=this->n_columns ) ) {
		fprintf ( stderr,  "array_link error: index out of bounds %d %d (%d %d)!!\n", r, c, this->n_rows, this->n_columns );
		return;
	}


	this->array[r+this->n_rows*c]= val;
}

void array_link::_setvalue ( int r, int c, int val )
{
	this->array[r+this->n_rows*c]= val;
}


array_t array_link::_at ( const rowindex_t r, const colindex_t c ) const
{
	return this->array[r+this->n_rows*c];
}

/**
 * @brief Element access
 */
array_t array_link::at ( const rowindex_t r, const colindex_t c ) const
{
//#ifdef OADEBUG
	if ( ( r<0 ) || ( r >= this->n_rows ) || ( c<0 ) || ( c>=this->n_columns ) ) {
		fprintf ( stderr,  "array_link error: index out of bounds %d %d (%d %d)!!\n", r, c, this->n_rows, this->n_columns );
		return 0;
	}
//#endif

	return this->array[r+this->n_rows*c];
}

/**
 * @brief Element access
 */
array_t & array_link::at ( const rowindex_t r, const colindex_t c )
{
#ifdef OADEBUG
	if ( ( r<0 ) || ( r >= this->n_rows ) || ( c<0 ) || ( c>=this->n_columns ) ) {
		fprintf ( stderr,  "array_link error: index out of bounds %d %d (%d %d)!!\n", r, c, this->n_rows, this->n_columns );
		return this->array[0];
	}
#endif

	return this->array[r+this->n_rows*c];
}

/**
 * @brief Element access
 */
void array_link::setconstant ( array_t value )
{
	std::fill ( this->array, this->array+n_rows*n_columns, value );
}

//! Create array link with clone of an array
array_link::array_link ( const array_t *array, rowindex_t nrows, colindex_t ncols, int index_ ) : n_rows ( nrows ), n_columns ( ncols ), index ( index_ )
{
	//printf("array_link::constructor: from array (index %d)\n", index);
	this->array = clone_array ( array, nrows, ncols );

	//initswig();

}

//! Create array link from vector
array_link::array_link ( const std::vector<int> &v, rowindex_t nrows, colindex_t ncols, int index_ ) : n_rows ( nrows ), n_columns ( ncols ), index ( index_ )
{
	if ( v.size() != ( size_t ) nrows*ncols ) {
		fprintf ( stderr, "array_link: error size of vector does not match number of rows and columns\n" );
		return;
	}
	//printf("array_link::constructor: from array (index %d)\n", index);
	this->array = create_array ( nrows, ncols );
	std::copy ( v.begin(), v.begin() +nrows*ncols, this->array );

	//initswig();
}

//! Default constructor
array_link::array_link() : n_rows ( -1 ), n_columns ( 0 ), index ( INDEX_NONE ), array ( 0 )
{
	//printf("array_link::default constructor\n");
}

void array_link::init ( rowindex_t r, colindex_t c )
{
#ifdef CLEAN_ARRAY_LINK
	if ( array!=0 ) {
		destroy_array ( array );
		array=0;
	}
#endif
	n_rows=r;
	n_columns=c;

	this->array = create_array ( n_rows, n_columns );

	//initswig();

}


//! Default destructor
array_link::~array_link()
{
	//	printf("~array_link: index %d\n", this->index);
	/* we do not destroy the array for efficiency in sorting, this should be done manually! */
#ifdef CLEAN_ARRAY_LINK
	if ( array!=0 ) {
		destroy_array ( array );
		array=0;
	}
#endif
}

//! Create array link with array
array_link::array_link ( rowindex_t nrows, colindex_t ncols, int index_ ) : n_rows ( nrows ), n_columns ( ncols ), index ( index_ )
{
	this->array = create_array ( nrows, ncols );

//	initswig();

}

//! Create array link with array
array_link::array_link ( rowindex_t nrows, colindex_t ncols, int index_, carray_t *data ) : n_rows ( nrows ), n_columns ( ncols ), index ( index_ )
{
// printf("array_link::array_link: nrows %d, ncols %d\n", nrows, ncols);
	this->array = create_array ( nrows, ncols );
	this->setarraydata ( data, nrows*ncols );

//	initswig();

}


/// return md5 sum of array representation (as represented with 32bit int datatype in memory)
std::string array_link::md5() const
{
	if ( sizeof ( array_t ) ==4 ) {
		// we have int as the representation type
		return ::md5 ( this->array, this->n_columns*this->n_rows*sizeof ( array_t ) );
	} else {
		short int nn = this->n_columns*this->n_rows;
		short int *x = new short int[nn];
		std::copy ( array, array+nn, x );
		std::string m = ::md5 ( x, nn*sizeof ( short int ) );
		delete x;
		return m;
	}
}

/**
 * @brief Print an array to a stream (formatted)
 * @param stream
 * @param A
 * @return
 */
std::ostream &operator<< ( std::ostream & stream, const array_link &A )
{
	write_array_format ( stream, A.array, A.n_rows, A.n_columns );
	return stream;
}

array_link array_link::selectFirstColumns ( int ii ) const
{
	// printf("  %d\n", this->n_columns);
	mycheck ( ii>=0, "array_link::selectFirstColumn: ii<0\n" );
	array_link d ( this->n_rows, ii, -1 );
	for ( int i=0; i<ii; i++ )
		std::copy ( this->array+i*this->n_rows, this->array+ ( i+1 ) *this->n_rows, d.array+i*this->n_rows );
	return d;
}

array_link array_link::selectLastColumns ( int ii ) const
{
	//printf("array_link::selectLastColumns  %d\n", this->n_columns);
	mycheck ( ii>=0, "array_link::selectFirstColumn: ii<0\n" );
	array_link d ( this->n_rows, ii, INDEX_DEFAULT );
	for ( int j=0; j<ii; j++ ) {
		int i=this->n_columns-ii+j;
		std::copy ( this->array+i*this->n_rows, this->array+ ( i+1 ) *this->n_rows, d.array+j*this->n_rows );
	}
	return d;
}

array_link array_link::selectColumns ( std::vector<int> c ) const
{
	array_link d ( this->n_rows, c.size(), INDEX_DEFAULT );
	for ( int j=0; j< ( int ) c.size(); j++ ) {
		int i=c[j];
		mycheck ( i>=0, "array_link::selectColumns: i<0\n" );
		mycheck ( i<this->n_columns, "array_link::selectColumns: i>ncols\n" );
		std::copy ( this->array+i*this->n_rows, this->array+ ( i+1 ) *this->n_rows, d.array+j*this->n_rows );
	}
	return d;
}

array_link array_link::deleteColumn ( int index ) const
{
	// printf("  %d\n", this->n_columns);
	mycheck ( this->n_columns>=2, "array_link::deleteColumn: array has <2 columns\n" );
	array_link d ( this->n_rows, this->n_columns-1, -1 );
	if ( index<0 )
		index=this->n_columns-1;
	for ( int i=0; i<index; i++ )
		std::copy ( this->array+i*this->n_rows, this->array+ ( i+1 ) *this->n_rows, d.array+i*this->n_rows );
	for ( int i=index+1; i<this->n_columns; i++ )
		std::copy ( this->array+i*this->n_rows, this->array+ ( i+1 ) *this->n_rows, d.array+ ( i-1 ) *this->n_rows );
	return d;
}


array_link array_link::transposed() const
{
	// printf("  %d\n", this->n_columns);
	array_link d ( this->n_columns, this->n_rows, -1 );
	for ( int i=0; i<this->n_columns; i++ )
		for ( int j=0; j<this->n_rows; j++ )
			d.array[i+j*d.n_rows] = 	this->array[j+i*this->n_rows];

	return d;
}

array_link array_link::randomperm() const
{
	arraydata_t arrayclass = arraylink2arraydata ( *this );
	array_transformation_t trans ( &arrayclass );
	trans.randomize();
	return trans.apply ( *this );
}


/** Return example array */
array_link exampleArray ( int idx, int verbose )
{

	std::string dstr = "";
	
	switch ( idx ) {
	default
			:
		printf ( "exampleArray: no such index %d", idx );
	case 0: {
		dstr ="array in OA(8,2, 2^2)";
		array_link al ( 8,2, 0 );
		std::vector<int> s;
		s.push_back ( 2 );
		s.push_back ( 2 );
		arraydata_t ad ( s, 8, 2, 2 );
		al.create_root ( ad );
		return al;
		break;
	}
	case 1: {
		// array 3 in OA(16, 2, 2^5)
		dstr = "array 3 in OA(16, 2, 2^5)" ;
		if ( verbose )
			printf ( "exampleArray: %s\n", dstr.c_str() );
		array_link al ( 16,5, 0 );
		int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		             1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0,
		             0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0,
		             0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0
		            };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 2: {
			dstr="array 6 in OA(16, 2, 2^6)";
		if ( verbose )
			printf ( "exampleArray: %s\n", dstr.c_str() );

		array_link al ( 16,6, 0 );
		int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		             1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
		             0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
		             1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0,
		             0, 1, 0, 1
		            };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 3: {
		// array 7 in OA(32, 3, 2^9)
		dstr="array 7 in OA(32, 3, 2^9)";
if ( verbose )
			printf ( "exampleArray: %s\n", dstr.c_str() );

		array_link al ( 32,7, 0 );
		int tmp[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
		             1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
		             1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
		             1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
		             1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		             1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
		             1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
		             0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
		             0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0,
		             1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1,
		             0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
		             1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1,
		             1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1
		            };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 4: {
		if ( verbose )
			printf ( "exampleArray: array 4 in OA(16, 2, 2^7)\n" );

		// array 4 in OA(16, 2, 2^7)
		array_link al ( 16,7, 0 );
		int tmp[] = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		              1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
		              0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
		              0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1,
		              0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0
		            };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 5: {
		if ( verbose )
			printf ( "exampleArray: array 0 in OA(24, 2, 4 3 2^a)\n" );

		// array 0 in OA(24, 2, 4 3 2^a)
		array_link al ( 24,5, 0 );
		int tmp[] = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3,
		              3, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1,
		              2, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0,
		              1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1,
		              1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0,
		              0, 1, 0, 1, 1
		            };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}

	case 6: {
		if ( verbose )
			printf ( "exampleArray: array in OA(4, 2, 2^a)\n" );

		// array in OA(4, 2, 2^a)
		array_link al ( 4,3, 0 );
		int tmp[] = { 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1 };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 7: {
		if ( verbose )
			printf ( "exampleArray: array 0 in OA(4, 2, 2^a)?\n" );

		//
		array_link al ( 4,3, 0 );
		int tmp[] = { 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0 };
		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 8: {
		if ( verbose )
			printf ( "exampleArray: array in OA(40, 3, 2^7)?\n" );

		//
		array_link al ( 40,7, 0 );
		int tmp[] = 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
		                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
		                 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		                 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
		                 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
		                 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
		                 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0,
		                 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0,
		                 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1,
		                 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0,
		                 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1,
		                 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1,
		                 0, 0, 1, 1
		             };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 10: {
		if ( verbose )
			printf ( "exampleArray: array in OA(9, 3^2)\n" );

		//
		array_link al ( 9,3, 0 );
		int tmp[] = 	{0,0,0,1,1,1,2,2,2,0,1,2,1,1,2,0,0,2,2,0,2,0,2,1,0,1,1     };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	case 11: {
		if ( verbose )
			printf ( "exampleArray: D-optimal array in OA(44, 2^8)\n" );

		//
		array_link al ( 44,8, 0 );
		int tmp[] = 	{1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
		                 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
		                 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0,
		                 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
		                 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1,
		                 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1,
		                 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0,
		                 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1,
		                 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1,
		                 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1,
		                 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0,
		                 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1,
		                 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0,
		                 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1,
		                 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1,
		                 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0
		             };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}

	case 12: {
		if ( verbose )
			printf ( "exampleArray: even-odd array OA(64, 2^13)\n" );

		//
		array_link al ( 64, 13, 0 );
		int tmp[] = 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		                 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0,
		                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		                 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
		                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
		                 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
		                 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
		                 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
		                 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1,
		                 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1,
		                 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
		                 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
		                 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
		                 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
		                 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
		                 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1,
		                 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1,
		                 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0,
		                 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
		                 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0,
		                 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1,
		                 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1,
		                 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0,
		                 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1,
		                 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0,
		                 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0,
		                 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1,
		                 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1,
		                 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1,
		                 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1,
		                 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
		                 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1,
		                 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1,
		                 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0,
		                 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1,
		                 1, 0, 0, 1
		             };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}

	case 13: {
		if ( verbose )
			printf ( "exampleArray: array in OA(25, 2^5)\n" );

		//
		array_link al ( 24,5, 0 );
		int tmp[] = 	{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
       1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0,
       0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0,
       1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0,
       1, 0, 1, 1, 0 };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}
	


	case 14: {
			dstr= "design in D(28, 2^5), D-efficiency is low" ;
		if ( verbose )
			printf ( "exampleArray: %s\n", dstr.c_str() );

		array_link al ( 28,5, 0 );
		int tmp[] = 	{1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1,
       0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1,
       1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1,
       1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0,
       0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1,
       1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0,
       0, 1 };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}

		case 15: {
			dstr= "design in D(56, 2^10), D-efficiency is low" ;
		if ( verbose )
			printf ( "exampleArray: %s\n", dstr.c_str() );

		array_link al ( 56,10, 0 );
		int tmp[] = 	{0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0,
       0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1,
       0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,
       0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0,
       0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1,
       0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1,
       1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0,
       0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0,
       0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1,
       0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0,
       0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0,
       0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0,
       0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1,
       1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1,
       0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0,
       1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0,
       0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0,
       1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
       0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
       0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1,
       0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1,
       0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
       1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0,
       0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1,
       0, 0, 0, 0, 1, 1, 0, 1
 };

		al.setarraydata ( tmp, al.n_rows*al.n_columns );
		return al;
		break;
	}

	
	
	}

	return array_link ( 1,1,-1 );
}

array_link array_link::reduceDOP() const
{
	array_link d = reduceDOPform ( *this, 0 );
	return d;
}

array_link array_link::reduceLMC() const
{
	int strength=this->strength();
	arraydata_t ad = arraylink2arraydata ( *this, 0, strength );
	LMCreduction_t reduction ( &ad );
	reduction.mode=OA_REDUCE;

	OAextend oaextend;
	//oaextend.setAlgorithmAuto(&ad);
	oaextend.setAlgorithm ( MODE_ORIGINAL, &ad );

	reduction.init_state=SETROOT;

	if ( 0 ) {
		printf ( "\n## array_link::reduceLMC()\n" );
		oaextend.info();
		printf ( "array_link::reduceLMC: strength %d, reduction.init_state %d \n", strength, reduction.init_state ); //al2.showarray();

		reduction.show();
	}

	LMCcheck ( *this, ad, oaextend, reduction );

	return reduction.getArray();

}

symmetry_group array_link::row_symmetry_group() const
{

	std::vector<mvalue_t<int> > rr;
	for ( int i=0; i<this->n_rows; i++ ) {
		mvalue_t<int> m;
		for ( int k=0; k<this->n_columns; k++ )
			m.v.push_back ( this->at ( i, k ) );
		rr.push_back ( m );
	}
	symmetry_group sg ( rr, true, 0 );
	return sg;
}


double array_link::nonzero_fraction() const
{
	double nz=0;
	long nn=this->n_columns*this->n_rows;
	for ( int i=0; i<nn; i++ )
		if ( this->array[i]!=0 )
			nz++;
	return nz/nn;
}

bool array_link::is2level() const
{
	int N = this->n_rows;
	for ( int r=0; r<this->n_rows; r++ ) {
		for ( int c=0; c<this->n_columns; c++ ) {
			if ( this->array[r+c*N]<0 )
				return false;
			if ( this->array[r+c*N]>1 )
				return false;
		}
	}
	return true;
}

void array_link::showproperties() const
{
	printf ( "array: %d rows, %d cols\n", this->n_rows, this->n_columns );
	printf ( "  strength %d, rank %d\n", this->strength(), this->rank() );
	printf ( "  D-efficiency %.3f\n", this->Defficiency() );
	std::vector<double> gwlp = this->GWLP();
	printf ( "  GWLP " );
	display_vector ( gwlp );
	printf ( "\n" );

//  int t = this->strength();
	return;
}

long array_link::data()
{
	//return static_cast<long>(array);
	return ( ( long ) array );
}
/*
void array_link::initswig()
{
	//printf("initswig! C side\n");
}
*/
std::string array_link::showarrayS() const
{
	std::stringstream ss;
	ss << "array: \n";
	write_array_format ( ss, array, this->n_rows, this->n_columns );
	return ss.str();
}

void array_link::showarraycompact() const
{
	for ( int r=0; r<this->n_rows; r++ ) {
		for ( int c=0; c<this->n_columns; c++ ) {
			printf ( "%d", this->_at ( r,c ) );
		}
		printf ( "\n" );
	}
}

void array_link::showarray() const
{
	printf ( "array: \n" );
	//write_array_format ( std::cout, array, this->n_rows, this->n_columns );	// does not work with ipython...
	write_array_format ( array, this->n_rows, this->n_columns );
//write_array(stdout, this->array, this->n_rows, this->n_columns);
}

void array_link::create_root ( const arraydata_t &ad )
{
	if ( ! ( ad.N<=this->n_rows ) ) {
		printf ( "array_link::create_root: number of columns too small for root of size %d\n", ad.N );
		return;
	}
	// printf("ad.strength %d, this->n_columns %d\n", ad.strength, this->n_columns );

	if ( ! ( ad.strength<=this->n_columns ) ) {
		printf ( "array_link::create_root: number of columns (%d) too small for specified strength %d\n", this->n_columns, ad.strength );
		return;
	}
	assert ( ad.strength<=this->n_columns );
	std::fill ( this->array, this->array+this->n_rows*this->n_columns, 0 );
	::create_root ( ( this->array ), &ad );
}

/*!
  Creates the root of an OA according to the strength, the number of rows and the levels per column.
  \brief Create Root
  \param array Pointer to OA, where the root is placed in
  \param ad Parameter with data of the array
  */
void create_root ( array_t *array, const arraydata_t *ad )
{
	int steps = ad->N;
	for ( colindex_t i = 0; i < ad->strength; i++ ) { /* loop over all root columns */
		// printf("create_root: i %d\n", i);
		steps /= ad->s[i];			//adjust nr of steps per collumn
		if ( steps==0 ) {
			//printf("create_root: steps=0, updating to 1\n");
			steps=1;
		}
		array_t l = 0;
		int coloffset = i * ad->N;
		for ( int j = 0; j < ad->N; j+= steps ) {	//big steps
			for ( int k = 0; k < steps; k++ ) { //small steps
				//printf("create_root: i j k %d %d %d\n", i, j,k);
				array[coloffset + j + k] = l;
			}
			l++;
			if ( l == ad->s[i] )
				l = 0;	//reset value, it has looped
		}
	}
	//    printf("create_root: done\n");
}

/*!
  Creates the root of an OA. The root is appended to the current list of arrays
  \brief Create Root
  \param ad Pointer to arraydata_t structure
  \param solutions List of current solutions
  \sa create_root
  */
void create_root ( const arraydata_t *ad, arraylist_t &solutions )
{
	array_link	cur_solution = array_link ( ad->N, ad->strength, 1 );
	create_root ( cur_solution.array, ad );
	solutions.push_back ( cur_solution );
}



std::vector<int> numberModelParams ( const array_link &al, int order=2 )

{
	int k = al.n_columns;
	std::vector<int> n ( order+1 );
	n[0]=1;
	n[1]=al.n_columns;

	if ( order>2 ) {
		printf ( "numberModelParams: not implemented for order > 2\n" );
	}
	arraydata_t arrayclass = arraylink2arraydata ( al, 0, 2 );
	std::vector<int> s = arrayclass.getS();
	std::vector<int> df = s;
	std::transform ( df.begin(), df.end(),  df.begin(), std::bind2nd ( std::minus<int>(), 1.0 ) );

	/* main effect contrasts */
	int mesize = std::accumulate ( df.begin(),df.end(),0 );

	/* 2fi*/
	int n2fi=0;
	for ( int ii=0; ii<k-1; ii++ ) {
		for ( int jj=ii+1; jj<k; jj++ ) {
			n2fi +=	df[ii]*df[jj];
		}
	}

	n[1]=mesize;
	n[2]=n2fi;
	return n;
}

MatrixFloat array_link::getModelMatrix ( int order, int intercept ) const
{
	int verbose=0;
	int N = this->n_rows;
	std::pair<MatrixFloat,MatrixFloat> mmx = array2eigenModelMatrixMixed ( *this, 0 );

	//std::cout << mmx.first;


	std::vector<int> np= numberModelParams ( *this, order );

	MatrixFloat intcpt =MatrixFloat::Zero ( N, 1 );
	intcpt.setConstant ( 1 );

	if ( order==2 ) {
		if ( verbose>=2 )
			printf ( "array_link::getModelMatrix %d+%d+%d\n", np[0], np[1], np[2] );
	}
	if ( order==0 ) {
		return intcpt;
	}
	if ( order==1 ) {
		MatrixFloat mm ( N, np[0]+np[1] );
		mm << intcpt, mmx.first;
		return mm;
	}
	if ( order==2 ) {
		MatrixFloat mm ( N, np[0]+np[1]+np[2] );
		//eigenInfo(mm, "mm");
		//std::cout << "gr\n" << mmx.first << std::endl << mmx.second << std::endl;
		mm << intcpt, mmx.first, mmx.second;
		//std::cout << "#### gr\n" << mm << std::endl;

		return mm;
	}

	printf ( "array_link::getModelMatrix: order > 2 not supported!\n" );
	return intcpt;
}


double array_link::CL2discrepancy() const
{
	return ::CL2discrepancy ( *this );

}

bool  array_link::foldover() const
{
	std::vector< double > g = this->GWLP();
	for ( size_t i=1; i<g.size(); i+=2 ) {
		if ( g[i]!=0 )
			return false;
	}
	return true;
}

std::vector< double > array_link::GWLP ( int truncate, int verbose ) const
{
	return ::GWLP ( *this, verbose, truncate );
}

/// convert array to second order interaction matrix in Eigen format
MatrixFloat array2eigenX2 ( const array_link &al )
{
	int k = al.n_columns;
	int n = al.n_rows;
	int m = 1 + k + k* ( k-1 ) /2;

	MatrixFloat mymatrix = MatrixFloat::Zero ( n,k* ( k-1 ) /2 );

	// init interactions
	int ww=0;
	for ( int c=0; c<k; ++c ) {
		for ( int c2=0; c2<c; ++c2 ) {
			int ci = c*n;
			int ci2 = c2*n;

			for ( int r=0; r<n; ++r ) {
				mymatrix ( r, ww ) = ( al.array[r+ci]+al.array[r+ci2] ) %2;
			}
			ww++;
		}
	}

	mymatrix.array() -= .5;
	mymatrix.array() *= 2;

	return mymatrix;
}


Eigen::VectorXd dummy()
{
	printf ( "dummy: create VectorXd\n" );
	Eigen::VectorXd v ( 3 );
	v[0]=1;
	v[1]=2;
	v[2]=4;
	return v;
}
Eigen::MatrixXd dummy2()
{
	printf ( "dummy2: create MatrixXd\n" );
	fflush ( stdout );
	Eigen::MatrixXd m ( 2,2 );
	m ( 0,0 ) = 3;
	m ( 1,0 ) = 2.5;
	m ( 0,1 ) = -1;
//	printf("dummy2: return MatrixXd\n"); fflush(stdout);
	return m;
}

// helper function for Python interface
void eigen2numpyHelper ( double* pymat1, int n, const Eigen::MatrixXd &m )
{
	//for(size_t i=0; i<n; i++) {
	//    pymat[i]=this->array[i];
	//}
	//eigenInfo ( m );
	printf ( "pymat1 %ld\n", ( long ) pymat1 );
	std::copy ( m.data(),m.data() +m.size(), pymat1 );
}

void eigenInfo ( const MatrixFloat m, const char *str, int verbose )
{
	if ( verbose==1 )
		printf ( "%s: %dx%d\n", str, ( int ) m.rows(), ( int ) m.cols() );
	if ( verbose==2 )
		printf ( "%s: %dx%d\n", str, ( int ) m.rows(), ( int ) m.cols() );

}
void eigenInfoF ( const Eigen::MatrixXf m, const char *str, int verbose )
{
	if ( verbose==1 )
		printf ( "%s: %dx%d\n", str, ( int ) m.rows(), ( int ) m.cols() );
	if ( verbose==2 )
		printf ( "%s: %dx%d\n", str, ( int ) m.rows(), ( int ) m.cols() );

}

MatrixFloat array2eigenME ( const array_link &al, int verbose )
{
	std::pair<MatrixFloat,MatrixFloat> mm = array2eigenModelMatrixMixed ( al, verbose );
	return mm.first;

}

// code from Eric Schoen, adapted to work for arrays of strength < 1
std::pair<MatrixFloat, MatrixFloat> array2eigenModelMatrixMixed ( const array_link &al, int verbose )
{
	//verbose=2;
	if ( verbose>=2 )
		printf ( "array2eigenModelMatrixMixed: start" );

	int N = al.n_rows;
	int k =  al.n_columns;
	arraydata_t arrayclass = arraylink2arraydata ( al, 0, 0 );
	std::vector<int> s = arrayclass.getS();

	std::vector<int> df = s;
	//printf("df: %d\n", df.size() );
	std::transform ( df.begin(), df.end(),  df.begin(), std::bind2nd ( std::minus<int>(), 1.0 ) );
	//printf("df: %d\n", df.size() );

	if ( verbose>=2 ) {
		arrayclass.show();
		printf ( "array2eigenME: N %d, k %d\n", N, k );
		printf ( "df " );
		printf_vector ( df, "%d " );
		printf ( "\n" );

	}
	MatrixFloat  AA = al.getEigenMatrix();

	/* main effect contrasts */
	if ( verbose>=2 )
		printfd ( "main effects\n" );

	int mesize = std::accumulate ( df.begin(),df.end(),0 );

	std::vector<int> np = numberModelParams ( al, 2 );
	MatrixFloat ME = MatrixFloat::Zero ( N, mesize );

	int meoffset=0;
	for ( int c=0; c<k; c++ ) {
		int md = df[c];
		MatrixFloat Z = MatrixFloat::Zero ( N, md+1 );	// large tmp buffer

		for ( int ii=0; ii<md+1; ii++ ) {
			for ( int r=0; r<N; r++ ) {
				Z ( r, ii ) =AA ( r, c ) > ( ii-1 );
			}
		}

		// make Helmert contrasts (these are automatically orthogonal)
		for ( int r=0; r<N; r++ ) {
			//printf("r: %d\n", r);
			int v = AA ( r,c );
			Z ( r, 0 ) = 1;
			if ( v>0 ) {
				Z ( r, v ) =v;
			}
			for ( int q=1; q<v; q++ )
				Z ( r,q ) =0;
			for ( int q=v+1; q<md+1; q++ )
				Z ( r,q ) =-1;

		}

		if ( verbose >=2 ) {
			eigenInfo ( Z , "Z first stage " );
			std::cout << Z << std::endl;

		}

		for ( int ii=0; ii<md; ii++ ) {
			if ( verbose>=3 ) {
				printf ( "array2eigenME: calculate ii %d\n", ii );

				eigenInfo ( Z.block ( 0,0,N,ii+1 ), "Zx" );
			}

			MatrixFloat tmp = Z.block ( 0,0,N,ii+1 ).transpose() * Z.block ( 0,0,N,ii+1 );
			MatrixFloat tmp2 = Z.block ( 0,0,N,ii+1 ).transpose() * Z.block ( 0,ii+1,N,1 ); // right part
			if ( verbose>=3 ) {
				eigenInfo ( tmp, "tmp" );
				std::cout << tmp << std::endl;
				eigenInfo ( tmp2, "tmp2" );
				std::cout << tmp2 << std::endl;
			}
			MatrixFloat b = tmp.colPivHouseholderQr().solve ( tmp2 );

			b *= 0;

			//Eigen::MatrixXd b =  tmp2; // should be ! tmp\tmp2;
			if ( verbose>=3 ) {
				eigenInfo ( Z.block ( 0,0,N,ii+1 ) , "Z.block(0,0,N,ii+1) " );
				eigenInfo ( b, "b" );
				std::cout << b << std::endl;
			}
			Z.col ( ii+1 ) -= Z.block ( 0,0,N,ii+1 ) * b;

			tmp=Z.col ( ii+1 ).transpose() * Z.col ( ii+1 );

			//eigenInfo ( tmp, "denom" );
			ME.col ( meoffset+ii ) =sqrt ( double ( N ) ) *Z.col ( ii+1 ) / sqrt ( double ( tmp ( 0,0 ) ) );
//			ME{jj}(:,ii)=sqrt(N)*Z{jj}(:,ii+1)/sqrt(Z{jj}(:,ii+1)'*Z{jj}(:,ii+1));

		}

		if ( verbose>=2 ) {
			eigenInfo ( Z, "Z" );
			std::cout << Z << std::endl;
		}
		meoffset+=md;



	}

	/* 2fi */
	if ( verbose>=2 )
		printf ( "2fi\n" );

	MatrixFloat tfi = MatrixFloat::Zero ( N, np[2] );

	if ( verbose>=2 )
		printfd ( "create 2fi\n" );
	int tel=0;
	int n = al.n_columns;
	int po=0, qo=0; // offsets
	for ( int ii=0; ii<n-1; ii++ ) {
		int n1 = df[ii];
		po = std::accumulate ( df.begin(), df.begin() +ii, 0 );

		for ( int jj=ii+1; jj<n; jj++ ) {
			int n2 = df[jj];

			qo = std::accumulate ( df.begin(), df.begin() +jj,0 );

			if ( verbose>=3 )
				printfd ( "create 2fi: p0=o %d, qo %d\n", po, qo );

			for ( int pp=0; pp<n1; pp++ ) {
				for ( int qq=0; qq<n2; qq++ ) {
					//if (verbose) printfd("create 2fi: \n");	eigenInfo(ME.col(pp), "ME.col(pp)");

					tfi.col ( tel ) =ME.col ( pp+po ).cwiseProduct ( ME.col ( qq+qo ) );
					tel++;
				}
			}
			qo +=n2;
		}
		po+=n1;
	}


	if ( verbose>=2 )
		printf ( "done\n" );


	return std::pair<MatrixFloat,MatrixFloat> ( ME, tfi );
}

MatrixFloat array2eigenX1 ( const array_link &al, int intercept )
{

	int k = al.n_columns;
	int n = al.n_rows;

	intercept = intercept>0;
	MatrixFloat mymatrix = MatrixFloat::Zero ( n, k+intercept );

	int ww=0;
	if ( intercept ) {

		// init first column
		for ( int r=0; r<n; ++r ) {
			mymatrix ( r, ww ) = 1;
		}
		ww+=1;
	}
	// init array
	for ( int c=0; c<k; ++c ) {
		int ci = c*n;
		for ( int r=0; r<n; ++r ) {
			mymatrix ( r, ww+c ) = al.array[r+ci];
		}
	}

	mymatrix.array() -= .5;
	mymatrix.array() *= 2;

	return mymatrix;
}

MatrixFloat array2eigenX1 ( const array_link &al )
{
	int k = al.n_columns;
	int n = al.n_rows;

	MatrixFloat mymatrix = MatrixFloat::Zero ( n, 1+k );

	// init first column
	int ww=0;
	for ( int r=0; r<n; ++r ) {
		mymatrix ( r, ww ) = 1;
	}

	// init array
	ww=1;
	for ( int c=0; c<k; ++c ) {
		int ci = c*n;
		for ( int r=0; r<n; ++r ) {
			mymatrix ( r, ww+c ) = al.array[r+ci];
		}
	}

	mymatrix.array() -= .5;
	mymatrix.array() *= 2;

	return mymatrix;
}

//Eigen::MatrixXd array2eigenME ( const array_link &al, int verbose );

Eigen::MatrixXi array2eigenModelMatrixInt ( const array_link &al )
{
	int k = al.n_columns;
	int n = al.n_rows;
	int m = 1 + k + k* ( k-1 ) /2;

	if ( n*k>0 ) {
		assert ( *std::max_element ( al.array, al.array+al.n_columns*al.n_rows ) <2 );
	}

	//MatrixFloat mymatrix = MatrixFloat::Zero ( n,m );
	//eigenFloat *data = mymatrix.data();
	
	// create data in integer type (we are working with 2-level arrays, convert them later */
	Eigen::MatrixXi mymatrix = Eigen::MatrixXi::Zero ( n,m );
	int *data = mymatrix.data();

	
	// init first column
	int ww=0;
	for ( int r=0; r<n; ++r ) {
		mymatrix ( r, ww ) = 1;
	}

	// init array
	ww=1;
	for ( int c=0; c<k; ++c ) {
		int ci = c*n;

		std::copy(al.array+ci, al.array+ci+n, data+(ww+c)*n);
		for ( int r=0; r<n; ++r ) {
			//mymatrix ( r, ww+c ) = al.array[r+ci];		
		}
	}

	// init interactions
	ww=k+1;

	for ( int c=0; c<k; ++c ) {
		for ( int c2=0; c2<c; ++c2 ) {
			int ci = c*n;
			int ci2 = c2*n;

			for ( int r=0; r<n; ++r ) {
				//mymatrix ( r, ww ) = ( al.array[r+ci]+al.array[r+ci2] ) %2;
				data[r+ww*n] = ( al.array[r+ci]+al.array[r+ci2] ) %2;
				//if (mymatrix(r,ww) != ( al.array[r+ci]+al.array[r+ci2] ) %2 ) exit(1);
			}
			ww++;
		}
	}

	mymatrix.array() *=2; mymatrix.array() -= 1;
	//mymatrix.array() -= .5; mymatrix.array() *= 2;

	return mymatrix;
	//return mymatrix;
}

MatrixFloat array2eigenModelMatrix ( const array_link &al )
{
	return array2eigenModelMatrixInt(al).cast<eigenFloat>();
}

double array_link::DsEfficiency ( int verbose ) const
{
	if ( ! this->is2level() ) {
		return this->Defficiencies() [1];
	}

	const array_link &al = *this;
	int k = al.n_columns;
	int k1 = al.n_columns+1;
	int n = al.n_rows;
	int m = 1 + k + k* ( k-1 ) /2;
	int N = n;

	//Eigen::MatrixXd X1 = array2eigenX1(al);
	MatrixFloat X2 = array2eigenX2 ( al );
	MatrixFloat X = array2eigenModelMatrix ( al );

#define SCALEN
#ifdef SCALEN
	MatrixFloat tmp = ( X.transpose() *X/n );
	double f1 = tmp.determinant();
	double f2 = ( X2.transpose() *X2/n ).determinant();
	double Ds = 0;
	if ( fabs ( f1 ) <1e-15 ) {
		if ( verbose )
			printf ( "DsEfficiency: f1 < 1e-15, setting Ds to zero\n" );
	} else {
		Ds=pow ( ( f1/f2 ), 1./k1 );
	}
	if ( verbose )
		printf ( "f1 %e, f2 %e, Ds %f\n", f1, f2, Ds );
#else
	MatrixFloat tmp = ( X.transpose() *X );
	double f1 = tmp.determinant();
	double f2 = ( X2.transpose() *X2 ).determinant();
	double Ds = pow ( ( f1/f2 ), 1./k1 ) /n;
	if ( verbose )
		printf ( "scale after: f1 %f, f2 %f, Ds %f\n", f1, f2, Ds );
#endif

	return Ds;
}


std::vector<double> array_link::Defficiencies ( int verbose, int addDs0 ) const
{
	const array_link &al = *this;

	arraydata_t arrayclass = arraylink2arraydata ( al );

	return ::Defficiencies ( al, arrayclass, verbose,  addDs0 );

}


double array_link::Defficiency() const
{
	if ( ! this->is2level() ) {
		return this->Defficiencies() [0];
	}

	return ::Defficiency ( *this );
}
double array_link::VIFefficiency() const
{
	return ::VIFefficiency ( *this );
}
double array_link::Aefficiency() const
{
	printf ( "warning: definition of Aefficiency has changed!\n" );
	return ::Aefficiency ( *this );
}

double array_link::Eefficiency() const
{
	return ::Eefficiency ( *this );
}

std::vector<int> array_link::Fvalues ( int jj ) const
{
	jstruct_t js ( *this, jj );
	std::vector<int> FF=js.calculateF();
	return FF;
}

std::vector<double> array_link::PECsequence() const
{
	return ::PECsequence ( *this );
}

std::vector<int> array_link::Jcharacteristics ( int jj ) const
{
	return ::Jcharacteristics ( *this, jj, 0 );
}

int array_link::rank() const
{
	return arrayrank ( *this );
}

int array_link::strength() const
{
	int t=-1;
	for ( int i=0; i<=this->n_columns; i++ ) {
		int b = strength_check ( *this, i );
		// printf("strength_check %d->%d\n", i, b);
		if ( b )
			t=i;
		else
			break;
	}
	return t;
}

int array_cmp ( carray_p A, carray_p B, const rowindex_t r, const colindex_t c )
{
	int count = 0;
	for ( int cpos=0; cpos<c; cpos++ ) {
		for ( int rpos=0; rpos<r; rpos++ ) {
			if ( A[count] < B[count] ) {
				return -1;
			}
			if ( A[count] > B[count] ) {
				return 1;
			}
			count++;
		}
	}
	return 0;
}

int array_diff ( carray_p A, carray_p B, const rowindex_t r, const colindex_t c, rowindex_t &rpos, colindex_t &cpos )
{
	int count = 0;
	for ( cpos=0; cpos<c; cpos++ ) {
		for ( rpos=0; rpos<r; rpos++ ) {
			if ( A[count] < B[count] ) {
				return -1;
			}
			if ( A[count] > B[count] ) {
				return 1;
			}
			count++;
		}
	}
	return 0;
}

/// create new arraydata_t object
arraydata_t::arraydata_t ( const array_t *s_, rowindex_t N_, colindex_t t, colindex_t nc ) : N ( N_ ), ncols ( nc ), strength ( t ), order ( ORDER_LEX ), colgroupindex ( 0 ), colgroupsize ( 0 )
{
	s = new array_t[nc];
	memcpy ( ( void * ) s, ( const void * ) s_, sizeof ( array_t ) *nc );
	complete_arraydata();
} ;
arraydata_t::arraydata_t ( const std::vector<int> s_, rowindex_t N_, colindex_t t, colindex_t nc ) : N ( N_ ), ncols ( nc ), strength ( t ), order ( ORDER_LEX ), colgroupindex ( 0 ), colgroupsize ( 0 )
{
	if ( ( int ) s_.size() <nc ) {
		fprintf ( stderr, "arraydata_t: warning: in constructor size s < number of columns nc)\n" );
		nc=s_.size();
	}

	//printf("arraydata_t: construct from std::vector\n");
	s = new array_t[nc];
	// for(int i=0; i<nc; i++) s[i]=s_[i];
	std::copy ( s_.begin(), s_.begin() +nc, s );
	complete_arraydata();
} ;

/// instantiate function
template void array_link::setarraydata ( const short int* tmp, int n );
template void array_link::setarraydata ( const int* tmp, int n );
template void array_link::setarraydata ( const long* tmp, int n );
template void array_link::setarraydata ( const std::vector<short int> tmp, int n );
//template void array_link::setarraydata ( const std::vector<int> tmp, int n );
template void array_link::setarraydata ( const std::vector<long> tmp, int n );

//   template void func<int>(int param);

arraydata_t::arraydata_t ( array_t s_, rowindex_t N_, colindex_t t, colindex_t nc ) : N ( N_ ), ncols ( nc ), strength ( t ), order ( ORDER_LEX ), colgroupindex ( 0 ), colgroupsize ( 0 )
{
	if ( s_<1 || s_>100 ) {
		fprintf ( stderr, "arraydata_t: level factors should be > 0 and < 100\n" );
	}
	s = new array_t[nc];
	for ( int i=0; i<nc; i++ )
		s[i]=s_;
	complete_arraydata();
} ;


arraydata_t::arraydata_t ( const arraydata_t *adp, colindex_t newncols )
{
	if ( adp==0 ) {
		printf ( "arraydata_t: error: pointer to arraydata_t is invalid\n" );
		return;
	}

	N=adp->N;
	strength=adp->strength;
	ncols=newncols;
	order=adp->order;
	colgroupindex=0;
	colgroupsize=0;

	if ( ncols>adp->ncols ) {
		printf ( "arraydata_t: error: number of columns %d is too large (object %d cols)\n", ncols, adp->ncols );
		s = new array_t[ncols];
		for ( size_t i=0; i< ( size_t ) ncols; i++ )
			s[i]=2;
		memcpy ( s, adp->s, sizeof ( array_t ) *adp->ncols );
	} else {
		s = new array_t[ncols];
		memcpy ( s, adp->s, sizeof ( array_t ) *ncols );
	}
	complete_arraydata();

	//printf("hack: created arraydata_t:\n"); this->show(3);
}

/// @brief copy constructor
arraydata_t::arraydata_t ( const arraydata_t &adp ) : N ( adp.N ), ncols ( adp.ncols ), strength ( adp.strength ), order ( adp.order ), colgroupindex ( 0 ), colgroupsize ( 0 )
{
	s = new array_t[ncols];
	memcpy ( s, adp.s, sizeof ( array_t ) *ncols );
	complete_arraydata();
}
arraydata_t::~arraydata_t()
{
	//printf("~arraydata_t\n");
	delete [] s;
	delete [] colgroupindex;
	delete [] colgroupsize;
}

void arraydata_t::writeConfigFile ( const char *file ) const
{
	arraydata_t ad = *this;
	//***open config file***
	colindex_t ncols = ad.ncols;
	std::ofstream outFile;

	outFile.open ( file );
	if ( !outFile ) {
		std::cerr << "Unable to open file " << file << std::endl;
		//throw -1;
		//exit(1); // terminate with error
	}

	/* read design specifications: runs, strength, number of factors */
	std::string str;
	outFile << "runs " << ad.N << std::endl;
	outFile << "strength " << ad.strength << std::endl;
	outFile << "nfactors " << ad.ncols << std::endl;

	for ( int j = 0; j < ncols; j++ ) {
		outFile << ad.s[j];
		if ( j<ncols-1 )
			outFile << " ";
		else
			outFile << std::endl;
	}
	outFile.close();
}

/**
 * @brief Show array data
 */
void arraydata_t::show ( int verbose ) const
{
	std::cout << showstr() << std::endl;
	if ( verbose>=2 ) {
		for ( int i=0; i<this->ncolgroups; i++ ) {
			printf ( " colgroupindex[%d] %d\n", i, this->colgroupindex[i] );
			printf ( " colgroup size[%d] %d\n", i, this->colgroupsize[i] );
		}
	}
}

/**
 * @brief Show array data
 */
std::string arraydata_t::showstr() const
{
	std::stringstream ss;
	ss << printfstring ( "arrayclass: N %d, k %d, strength %d, s ", this->N, this->ncols, this->strength );
	print_perm ( ss, this->s, this->ncols );
	std::string s = ss.str();
	s=s.substr ( 0, s.size()-1 );
	s += printfstring ( ", order %d", this->order );
	return s;
}

/**
 * @brief Show array data
 */
std::string arraydata_t::latexstr ( int cmd, int series ) const
{
	std::stringstream ss;
	if ( cmd ) {
		//ss << printfstring ( "\\OA(%d, %d, ", this->N, this->strength );
		ss << printfstring ( "\\oadesign{%d}{%d}{", this->N, this->strength );
	} else
		ss << printfstring ( "\\mathrm{OA}(%d, %d, ", this->N, this->strength );

	for ( int i=0; i<this->ncolgroups; ++i ) {
		//int s = this->s[this->colgroupindex[i]];
		//printf("this->colgroupindex[i] %d \n", this->colgroupindex[i]);
		int cgi=this->colgroupindex[i];
		int s=this->s[cgi];
		if ( series>0 && i==this->ncolgroups-1 &&  this->colgroupsize[i]>1 ) {
			ss << printfstring ( "%d^{a}", s );
		} else {
			ss << printfstring ( "%d^{%d}", s, this->colgroupsize[i] );
		}
	}
	if ( cmd )
		ss << printfstring ( "}" );
	else
		ss << printfstring ( ")" );
	return ss.str();
}

/**
 * @brief Return identifier string
 */
std::string arraydata_t::idstrseriesfull() const
{
	std::string fname = "";
	fname += itos ( N );

//    char xstr="abcdefg";

	for ( int i=0; i<this->ncolgroups; ++i ) {
		int cgi=this->colgroupindex[i];
		int s=this->s[cgi];
		if ( i==0 )
			fname += '.';
		else
			fname += '-';
		if ( i==this->ncolgroups-1 &&  this->colgroupsize[i]>1 ) {
			fname += printfstring ( "%da", s );
		} else {
			fname += printfstring ( "%d^%d", s, this->colgroupsize[i] );
		}
	}
	fname += "-t" + printfstring ( "%d", this->strength );
	return fname;
};

/// return full identifier string
std::string arraydata_t::fullidstr ( int series ) const
{
	if ( series )
		return this->idstrseriesfull(); // + "-t" + printfstring ( "%d", this->strength );
	else
		return this->idstr() + "-t" + printfstring ( "%d", this->strength );
}

/**
 * @brief Return identifier string
 */
std::string arraydata_t::idstr() const
{
	std::string fname = "";
	fname += itos ( N );

	for ( int i = 0; i < ncols; i++ ) {
		if ( i==0 )
			fname += '.';
		else
			fname += '-';
		fname += itos ( s[i] );
	}
	return fname;
};


/**
 * @brief Calculate derived data such as the index and column groups from a design
 */
void arraydata_t::complete_arraydata()
{
	const int verbose=0;

	if ( verbose ) {
		for ( int i=0; i<this->ncols; i++ )
			printf ( "complete_arraydata: k %d, s %d\n", i, s[i] );
	}
	if ( this->strength>this->ncols ) {
		printf ( "arraydata_t: warning strength %d > ncols %d, reducing strength\n", this->strength, this->ncols );
		this->strength=this->ncols;
	}
	if ( this->strength<1 ) {
		if ( verbose>=2 ) {
			printf ( "arraydata_t: warning strength < 1\n" );
		}
		//this->strength=1;
	}
	arraydata_t *ad = this;
	this->calcoaindex ( ad->strength );

	/* calculate column group structure */
//    colindex_t *dummy;
	std::vector<int> xx ( ad->s, ad->s+ad->ncols );

	symmetry_group sg ( xx, 0 );
	//printf("-- arraydata_t::complete_arraydata\n"); sg.show(2);

	ad->ncolgroups = sg.ngroups;
	//printf("ncolgroups %d\n", ad->ncolgroups);
	ad->colgroupindex = new colindex_t[ad->ncolgroups+1];
	std::copy ( sg.gstart.begin(), sg.gstart.end(), ad->colgroupindex );
	ad->colgroupsize = new colindex_t[ad->ncolgroups+1];
	std::copy ( sg.gsize.begin(), sg.gsize.end(), ad->colgroupsize );

// check
	int nbits = 8*sizeof ( rowsort_value_t );
	rowsort_value_t val=1;
	for ( int i=0; i<ad->ncols; i++ ) {
		if ( val != 0 && ( std::numeric_limits<rowsort_value_t>::max() / ( rowsort_value_t ) ad->s[i] ) < val ) {
			// multiplication would exceed range of unsigned
			printf ( "error: LMC checks for %d columns would lead to integer overflow\n", i );
			std::cout << "  column: "<< i << ": max rowsort value " << val << printfstring ( " (ncols %d, nbits %d)", ad->ncols, nbits ) << std::endl;
			std::cout << printfstring ( "      (ncols %d, nbits %d, ad->s[i] %d)", ad->ncols, nbits, ad->s[i] ) <<  std::endl;

		}
		val *= ad->s[i];

	}
	//std::cout << "max rowsort value: "<< val << printfstring(" (ncols %d, nbits %d)", ad->ncols, nbits) <<  std::endl;

//ad->ncolgroups = symm_group_index<array_t, colindex_t>(ad->s, ad->ncols, dummy, ad->colgroupindex, ad->colgroupsize);
	//delete [] dummy;
}

array_link arraydata_t::randomarray ( int strength, int ncols ) const
{

	if ( ncols==-1 )
		ncols=this->ncols;
	array_link al ( this->N, this->ncols, -1 );
	al.setconstant ( 0 );

	//al.show(); printf("----\n"); al.showarray();
	for ( int i=0; i<this->ncols; i++ ) {
		int coloffset = this->N*i;
		array_t s = this->getfactorlevel ( i );

		int step = floor ( double ( N ) /s );
		//printf("randomarray: col %d: s %d, step %d\n", i, s, step );
		if ( strength==1 ) {
			for ( int j=0; j<s; j++ ) {
				std::fill ( al.array+coloffset+step*j, al.array+coloffset+step* ( j+1 ), j );
			}
			random_perm ( al.array+coloffset, this->N );

		} else {
			for ( int r=0; r<N; r++ ) {
				al.array[r+coloffset]=fastrand() % s;
			}
		}
	}
	//printf("----\n"); al.showarray();
	return al;
}

array_link arraydata_t::create_root() const
{
	array_link al ( this->N, this->strength, -1 );
	al.create_root ( *this );
	return al;
}

bool arraydata_t::is2level() const
{
	for ( int i=0; i<this->ncols; i++ ) {
		if ( this->s[i] != 2 )
			return false;
	}
	return true;
}

bool arraydata_t::ismixed() const
{
	colindex_t s=this->s[0];
	for ( int i=0; i<this->ncols; i++ ) {
		if ( s!= this->s[i] )
			return true;
	}
	return false;
}

void arraydata_t::set_colgroups_jj ( const symmetry_group &sg, int jj )
{
	delete [] colgroupsize;
	delete [] colgroupindex;

//   mycheck ( sg.n == this->ncols,  "arraydata_t::set_colgroups: invalid size" );
	this->ncolgroups = sg.ngroups+1;
	this->colgroupindex = new colindex_t[this->ncolgroups+1];
	std::copy ( sg.gstart.begin(), sg.gstart.end(), this->colgroupindex );
	this->colgroupindex[sg.ngroups]=jj;
	this->colgroupsize = new colindex_t[this->ncolgroups+1];
	std::copy ( sg.gsize.begin(), sg.gsize.end(), this->colgroupsize );
	this->colgroupsize[sg.ngroups]=this->ncols-jj;
}


void arraydata_t::set_colgroups ( const std::vector<int> splits )
{
	delete [] colgroupsize;
	delete [] colgroupindex;

//   mycheck ( sg.n == this->ncols,  "arraydata_t::set_colgroups: invalid size" );
	this->ncolgroups = splits.size();
	this->colgroupindex = new colindex_t[this->ncolgroups+1];
	std::copy ( splits.begin(), splits.end(), this->colgroupindex );
	this->colgroupsize = new colindex_t[this->ncolgroups+1];
	for ( int j=0; j< ( ncolgroups-1 ); j++ )
		this->colgroupsize[j]=this->colgroupindex[j+1]-this->colgroupindex[j];
	this->colgroupsize[ncolgroups-1]=this->ncols-this->colgroupindex[ncolgroups-1];
}

void arraydata_t::set_colgroups ( const symmetry_group &sg )
{
	delete [] colgroupsize;
	delete [] colgroupindex;

	mycheck ( sg.n == this->ncols,  "arraydata_t::set_colgroups: invalid size sg.n %d, ncols %d", sg.n, this->ncols );
	this->ncolgroups = sg.ngroups;
	this->colgroupindex = new colindex_t[this->ncolgroups+1];
	std::copy ( sg.gstart.begin(), sg.gstart.end(), this->colgroupindex );
	this->colgroupsize = new colindex_t[this->ncolgroups+1];
	std::copy ( sg.gsize.begin(), sg.gsize.end(), this->colgroupsize );
}

/**
 * @brief Complete arraydata but treat last column as single column group
 */
void arraydata_t::complete_arraydata_splitn ( int ns )
{
	delete [] colgroupsize;
	delete [] colgroupindex;

	for ( int i=0; i< ( this->ncols-1 ); i++ ) {
		if ( this->s[i]!=this->s[i+1] ) {
			printf ( "complete_arraydata_splitn: s[%d]=%d, s[%d]=%d, NOT IMPLEMENTED\n", i, this->s[i], i+1, this->s[i+1] );
			assert ( 0 );
		}
	}


	arraydata_t *ad = this;
	int combs = 1;
	for ( int i=0; i<ad->strength; i++ )
		combs *= ad->s[i];
	ad->oaindex = ad->N/combs;
	ad->ncolgroups=2;

	//ad->colgroupindex = new colindex_t[ad->N];
	colgroupindex = new colindex_t[ad->ncolgroups+1];
	colgroupsize = new colindex_t[ad->ncolgroups+1];
	colgroupindex[0]=0;
	colgroupindex[1]=ns;
	colgroupsize[0]=ns;
	colgroupsize[1]=ad->ncols-ns;


}

/**
 * @brief Complete arraydata but treat last column as single column group
 */
void arraydata_t::complete_arraydata_fixlast()
{
	delete [] colgroupsize;
	delete [] colgroupindex;

	complete_arraydata();

	if ( colgroupsize[ncolgroups-1]==1 ) {
		/* last column already a single column group */
		return;
	} else {
		colgroupsize[ncolgroups-1]--;
		colgroupindex[ncolgroups]=ncols-1;
		colgroupsize[ncolgroups]=1;
		ncolgroups = ncolgroups+1;
	}
}


/* analyse arrays */

jstruct_t::jstruct_t()
{
// printf("jstruct_t()\n");
	this->nc=0;
	this->A=-1;
}

std::vector<int> jstruct_t::Fval ( int strength ) const
{
	int x=pow ( ( double ) 2, strength+1 );	// TODO: replace by integer power
	int nn = floor ( ( double ) N/x ) +1;
	std::vector<int> Fv ( nn );
	for ( int i=0; i<nn; i++ ) {
		Fv[i]=N-x*i;
	}
	return Fv;
}

std::vector<int> jstruct_t::calculateF ( int strength ) const
{
	int Nmax=N;

	int x=pow ( double ( 2 ), strength+1 ); // TODO: replace by integer power
	int nn = floor ( ( double ) N/x ) +1;
	std::vector<int> F ( nn );

	for ( int i=0; i<nc; i++ ) {
		int fi = ( N-abs ( vals[i] ) ) /x;
		//if (fi>=nn)
		//  printf("  k %d, val %d, x %d, fi %d\n",  k, vals[i], x, fi);
		F[fi]++;
	}
//		printf("  x\n");
	return F;
}

void jstruct_t::calc ( const array_link &al )
{
	int *pp = new_perm_init<int> ( jj );
	//int ncolcombs = ncombs ( k, jj );

	for ( int x=0; x<this->nc; x++ ) {
		int jv = jvalue ( al, jj, pp ); // TODO: this is slow, we already have better code for this
		this->vals[x]=jv;
		next_comb_s ( pp, jj, k );
		//  cout << printfstring("   J value %d\n", jv);
	}

	delete_perm ( pp );
}

void jstruct_t::calcj4 ( const array_link &al )
{
	// printf ( "jstruct_t::calcj4\n" );
	assert ( jj==4 );
	int *pp = new_perm_init<int> ( jj );
	int ncolcombs = ncombs ( k, jj );
	const int nr = al.n_rows;
	const int N = nr;

	// fill double column table
	array_link dtable ( nr, al.n_columns*al.n_columns, -1 );

	// loop over all column pairs
	for ( int i=0; i<al.n_columns; i++ ) {
		for ( int j=0; j<al.n_columns; j++ ) {
			int idx=i+j*al.n_columns;
			// loop over all rows of original array
			for ( int x=0; x<nr; x++ ) {
				//dtable.array[x+idx*dtable.n_rows] = al.atfast ( x, i ) + al.atfast ( x, j );
				dtable.array[x+idx*dtable.n_rows] = al.array[x+al.n_rows*i] + al.array[x+al.n_rows*j];
			}

		}
	}

	//printf("dtable\n"); dtable.showarray();

	for ( int x=0; x<this->nc; x++ ) {
		//int jv = jvalue ( al, jj, pp );  this->vals[x]=jv;
		const int idx1= pp[0]+pp[1]*al.n_columns;
		const int idx2= pp[2]+pp[3]*al.n_columns;
		int jv=0;
		if ( 0 ) {
			for ( int xr=0; xr<nr; xr++ ) {
				int tmp = dtable.atfast ( xr, idx1 ) +dtable.atfast ( xr, idx2 );
				tmp %= 2;
				jv += tmp;
				// tmp *=2; tmp--;jv -= tmp;
			}
		} {
			const array_t *o1 = dtable.array+dtable.n_rows*idx1;
			const array_t *o2 = dtable.array+dtable.n_rows*idx2;
			for ( int xr=0; xr<nr; xr++ ) {
				int tmp = o1[xr]  + o2[xr]; // dtable.atfast ( xr, idx2 );
				tmp %= 2;
				jv += tmp;
			}
//			jv %= 2;
		}
		jv = 2*jv-N;
		this->vals[x]=jv;
		//printf(" val %d -> %d\n", jvalue ( al, jj, pp ), jv); printf("  perm "); print_perm(pp, jj);

		next_comb_s ( pp, jj, k );
		//xi++;
		//  cout << printfstring("   J value %d\n", jv);
	}

	delete_perm ( pp );
}


jstruct_t::jstruct_t ( const array_link &al, int jj )
{
	const int k=al.n_columns;
	const int N = al.n_rows;


	this->init ( N, k, jj );
	if ( jj==4 && 1 )
		this->calcj4 ( al );
	else
		this->calc ( al );

	// calculate A value
	this->calculateAberration();

}

jstruct_t::jstruct_t ( const int N_, const int k_, const int jj_ )
{
	init ( N_, k_, jj_ );
}

void jstruct_t::init ( int N_, int k_, int jj_ )
{
	this->N = N_;
	this->k = k_;
	this->jj = jj_;

	if ( this->jj<0 ) {
		printf ( "jstruct_j: J-characteristics for J<0 are undefined, setting J to 0\n" );
		jj=0;
	}
	if ( this->jj>20 ) {
		printf ( "jstruct_j: J-characteristics for J>20 are not supported, setting J to 0\n" );
		jj=0;
	}

	this->nc = ncombs ( k_, jj_ );
	vals = std::vector<int> ( nc );
//  printf("jstruct_t(N,k,jj): vals %d\n", this->vals);
	this->A=-1;
}

jstruct_t::jstruct_t ( const jstruct_t &js )
{
	N = js.N;
	k = js.k;
	jj = js.jj;
	nc = js.nc;
	A=js.A;
	vals = std::vector<int> ( nc );
	std::copy ( js.vals.begin(), js.vals.begin() +nc, vals.begin() );
	// printf("jstruct_t(): copy constructor: js.vals %d, new %d\n", js.vals, vals);
}

jstruct_t &jstruct_t::operator= ( const jstruct_t &rhs )
{
	//printf("jstruct_t::operator=\n");
	this->N = rhs.N;
	this->k = rhs.k;
	this->jj = rhs.jj;
	this->nc = rhs.nc;

	this->A = rhs.A;
	vals = std::vector<int> ( nc );
	std::copy ( rhs.vals.begin(), rhs.vals.begin() +nc, vals.begin() );
	// printf("jstruct_t(): copy constructor: js.vals %d, new %d\n", js.vals, vals);

	return *this;
}


jstruct_t::~jstruct_t()
{
	//printf("~jstruct_t(): vals %d\n", this->vals);
}

std::string jstruct_t::showstr()
{
	std::stringstream ss;
	ss << "jstruct_t: " << printfstring ( "N %d, jj %d", N, jj );
	return ss.str();
}

void jstruct_t::show()
{
	cout << "jstruct_t: " << printfstring ( "N %d, jj %d, values ", N, jj );
	for ( int x=0; x<this->nc; x++ ) {
		cout << printfstring ( " %d", vals[x] );
	}
	cout << std::endl;
}

void jstruct_t::showdata()
{
	for ( int x=0; x<this->nc; x++ ) {
		cout << printfstring ( " %d", vals[x] );
	}
	cout << std::endl;
}


/** Calculate J-value for an array
 *
 * The array should be a 2-factor array
 *
 * @param array_link Array
 * @param J Number of columns
 * @param pp Indices of columns to use
 * @return j-value
 */
int jvalue ( const array_link &ar,const int J, const int *pp )
{
	int jval=0;

	for ( rowindex_t r=0; r<ar.n_rows; r++ ) {
		int tmp=0;
		for ( int i=0; i<J; i++ ) {
			tmp+=ar.at ( r,pp[i] );
		}
		tmp %= 2;
		tmp *=2;
		tmp--;
		jval -= tmp;
	}
	return ( jval );
}

/** Analyse a list of arrays
 *
 * Currently only j-values are calculated
 *
 */
vector<jstruct_t> analyseArrays ( const arraylist_t &arraylist,  const int verbose, const int jj )
{
	//vector<jstruct_t> results(arraylist.size());

	if ( verbose ) {
		printf ( "analyseArrays (j-values): %ld arrays, jj %d\n", ( long ) arraylist.size(), jj );
	}

	vector<jstruct_t> results;
	results.reserve ( arraylist.size() );
	jstruct_t *js ;

	//printf("analyseArrays: start\n");
	for ( unsigned int ii=0; ii<arraylist.size(); ii++ ) {

		const array_link ll = arraylist.at ( ii );
		const int k=ll.n_columns;

		//const int N = ll.n_rows;
		int *pp = new_perm_init<int> ( jj );
		int ncolcombs = ncombs ( k, jj );

		js = new jstruct_t ( ll, jj );

		if ( verbose>=3 ) {
			cout << printfstring ( "array %d: abberation %.3f j-values ", ii, js->A );
			print_perm ( cout, js->vals, js->nc );
		}
		if ( verbose>=2 ) {
			std::vector<int> FF=js->calculateF();
			printf ( "F%d (high to low): ", jj );
			display_vector ( FF );
			std::cout << std::endl;
		}
		delete_perm ( pp );
		//delete [] js;

		results.push_back ( *js );
		///results[ii]=*js; printf("assign\n");
		delete js;
	}

	// printf("analyseArrays: return value\n");
	return results;
}




/* reading and writing of arrays */

void write_array_latex ( FILE *fid, carray_t *array, const int nrows, const int ncols );

/**
 * @brief Write an array to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void write_array ( FILE *fid, carray_t *array, const int nrows, const int ncols )
{
	int count;
	stringstream ss ;

	for ( int j = 0; j < nrows; j++ ) {
		count = j;
		for ( int k = 0; k < ncols; k++ ) {
			const char *s = ( k<ncols-1 ) ? " ":"\n";
			ss << array[count] << s;
			count += nrows;
		}
	}
	fprintf ( fid, "%s", ss.str().c_str() );
}

/**
 * @brief Write array file in LaTeX format
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void write_array_latex ( FILE *fid, carray_t *array, const int nrows, const int ncols )
{
	int count;
	stringstream ss ;


	ss << "\\begin{tabular}{";
	for ( int x=0; x<ncols; x++ )
		ss << 'c';
	ss << "}" << endl;

	for ( int j = 0; j < nrows; j++ ) {
		count = j;
		for ( int k = 0; k < ncols; k++ ) {
			const char *s = ( k<ncols-1 ) ? " & ":" \\\\ \n";
			ss << array[count] << s;
			count += nrows;
		}
	}
	ss << "\\end{tabular}" << endl;

	fprintf ( fid, "%s", ss.str().c_str() );
}


/// dummy constructor
/*   arrayfile_t::arrayfile_t()
    {
        filename="";
        iscompressed=0;
        nfid=0;
        gzfid=0;
        nrows=-1;
        ncols=-1;
        narrays=-1;
        narraycounter=0;
        nbits=-1;
        mode=AERROR;

        verbose=1;
    }
    */
bool arrayfile_t::isbinary() const
{
	return ( this->mode==ABINARY || this->mode==ABINARY_DIFF || this->mode==ABINARY_DIFFZERO );
}


void arrayfile_t::append_array ( const array_link &a, int specialindex )
{
	int r;
	if ( ! this->isopen() ) {
		return;
	}

	arrayfile_t *afile = this;
	int index;
	if ( specialindex==-1 )
		index=a.index;
	else
		index=specialindex;

	switch ( afile->mode ) {
	case arrayfile::ATEXT:
		fprintf ( afile->nfid, "%i\n", index );
		write_array ( afile->nfid, a.array, a.n_rows, a.n_columns );
		break;
	case arrayfile::ALATEX:
		fprintf ( afile->nfid, "%i\n", index );
		write_array_latex ( afile->nfid, a.array, a.n_rows, a.n_columns );
		break;
	case arrayfile::ABINARY: {
		r = fwrite ( &index, 1, sizeof ( int32_t ), afile->nfid );
		if ( r!=sizeof ( int32_t ) *1 ) {
			printfd ( "error during write to file\n" );
		}
		//printf("   %d %d\n", a.n_rows, a.n_columns);
		afile->write_array_binary ( a.array, a.n_rows, a.n_columns );
		break;
	}
	case arrayfile::ABINARY_DIFF: {
		fwrite ( &index, 1, sizeof ( int32_t ), afile->nfid );
		afile->write_array_binary_diff ( a );
		break;
	}
	case arrayfile::ABINARY_DIFFZERO: {
		//fwrite ( &index, 1, sizeof ( int32_t ), afile->nfid );
		afile->write_array_binary_diffzero ( a );
		break;
	}
	default
			:
		std::cout << "warning: arrayfile_t::append_array: no such mode " << afile->mode << std::endl;
		break;
	}
	narraycounter++;
}

/// return true if file is open
int arrayfile_t::isopen() const
{
#ifdef USEZLIB
	return ( this->nfid != 0 || this->gzfid != 0 );
#else
	return ( this->nfid != 0 );
#endif
}

int arrayfile_t::append_arrays ( const arraylist_t& arrays, int startidx )
{
	arraylist_t::const_iterator	it;

	assert ( this->nfid );

	for ( it = arrays.begin(); it != arrays.end(); it++ ) {

		//cout << "writing arrays in mode " << afile->mode << endl;

		this->append_array ( *it, startidx );
		startidx++;
	}

	return 0;

//   return ::append_arrays ( this, arrays, startidx );
}


/**
 * @brief Free all the arrays in a list of arrays
 *
 * The memory allocation by the arrays is not freed in the default destructor for efficiency reasons.
 *
 * @param solutions
 * @return
 */
int free_sols ( arraylist_t &solutions )
{
	arraylist_t::iterator	it;

	log_print ( DEBUG, "Freeing %d solutions\n", solutions.size() );
	for ( it = solutions.begin(); it != solutions.end(); it++ ) {
		destroy_array ( it->array );
		it->array = NULL;
	}
	solutions.clear();
	return 0;
}


int append_arrays ( FILE *fid, arraylist_t &arrays, int startidx = 0 )
{
	arraylist_t::iterator	it;

	for ( it = arrays.begin(); it != arrays.end(); it++ ) {
		fprintf ( fid, "%i\n", startidx++ );
		write_array ( fid, it->array, it->n_rows, it->n_columns );
	}

	return 0;
}

template <class TypeIn, class TypeOut>
/// Write array to binary blob of selected datatype
void writeblob ( const TypeIn *src, int n, FILE *fid )
{
	TypeOut *dst = new TypeOut [n];
//		printf(" src %ld, dst %ld, n %d\n", long(src), long(dst), n);
	//printf("dst[0] %d, sizeof %d\n", dst[0], (int)sizeof(dst[0]));


	for ( int i=0; i<n; i++ ) {
		//printf("src[%d] %d, \n", i, src[i]);

		//printf("dst[%d] %d, \n", i, dst[i]);
		dst[i] = src[i];
	}
	fwrite ( ( const void * ) dst, sizeof ( TypeOut ), n, fid );

	delete [] dst;
}

template <class TypeIn, class TypeOut>
/// Convert binary blob to datatype
void readblob ( TypeOut *dst, int n, FILE *fid )
{
	TypeIn *src = new TypeIn [n];
	int r = fread ( ( void * ) src, sizeof ( TypeIn ), n, fid );

	for ( int i=0; i<n; i++ )
		dst[i] = src[i];

	delete [] src;
}
#ifdef USEZLIB
template <class TypeIn, class TypeOut>
/// Convert binary blob to datatype
void readblob ( TypeOut *dst, int n, gzFile gzfid )
{
	TypeIn *src = new TypeIn [n];
	gzread ( gzfid, src, sizeof ( TypeIn ) *n );

	for ( int i=0; i<n; i++ )
		dst[i] = src[i];

	delete [] src;
}
#endif

/**
 * @brief Read a binary array from a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::read_array_binary ( array_t *array, const int nrows, const int ncols )
{
//	  printf("  arrayfile_t::read_array_binary: gr %d\n", this->nbits);
	switch ( this->nbits ) {

	case 1: {

		// construct bit array
		BIT_ARRAY* bitarr = bit_array_create ( nrows*ncols );
		word_addr_t num_of_words = nwords ( nrows*ncols );
		int result = afread ( bitarr->words,num_of_words,sizeof ( word_t ) );
		assert ( result );

		// fill bit array
		for ( int i=0; i<nrows*ncols; i++ ) {
			// TODO: this can be made much faster by parsing multiple bits
			array[i] = bit_array_get_bit_nocheck ( bitarr, i );
		}
		bit_array_free ( bitarr );

		//printf("1-bit read not implemented yet\n");

		break;
	}
	case 8:
#ifdef USEZLIB
		if ( this->iscompressed )
			readblob<char, array_t> ( array, nrows*ncols, this->gzfid );
		else
#endif
			readblob<char, array_t> ( array, nrows*ncols, this->nfid );
		break;
	case 32:
#ifdef USEZLIB
		if ( this->iscompressed )
			readblob<int32_t, array_t> ( array, nrows*ncols, this->gzfid );
		else
#endif
			readblob<int32_t, array_t> ( array, nrows*ncols, this->nfid );
		break;
	default
			:
		printf ( " no such number of bits supported\n" );
		break;
	}
}

/**
 * @brief Read an array from a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void read_array ( FILE *fid, array_t *array, const int nrows, const int ncols )
{
	const int maxbuf = 1024;
	char buf[maxbuf];
	int count;
	for ( int j = 0; j < nrows; j++ ) {
		count = j;

		char *r = fgets ( buf, maxbuf, fid );
		stringstream ss ( stringstream::in | stringstream::out );
		ss << buf;

		for ( int k = 0; k < ncols; k++ ) {
			ss >> array[count];
			count += nrows;
		}

	}
}


/**
 * @brief Create array file
 * @param fname
 * @param rows
 * @param cols
 * @param narrays
 * @return
 */
arrayfile_t* create_arrayfile ( const char *fname, int rows, int cols, int narrays, arrayfile::arrayfilemode_t mode, int nbits )
{
	//printf("create array file: mode %d\n", mode);
	std::string s = fname;
	arrayfile_t *afile = new arrayfile_t ( s, rows, cols, narrays, mode, nbits );

	return afile;
}


void arrayfile_t::closefile()
{
	if ( verbose>=2 )
		printf ( "arrayfile_t::closefile(): nfid %ld\n", ( long ) nfid ) ;

	if ( ! this->isopen() ) {
		if ( verbose>=2 )
			printf ( "arrayfile_t::closefile(): file already closed\n" ) ;
		return;
	}
	//   printf ( "arrayfile_t: closefile(): filename %s, nfid %ld, narrays %d, narraycounter %d, this->rwmode %d\n", filename.c_str(), ( long ) nfid, narrays, narraycounter, this->rwmode );

	if ( verbose>=3 ) {
		printf ( "arrayfile_t::closefile: narrays: %d\n", narrays );
		printf ( "arrayfile_t::closefile: narraycounter: %d\n", narraycounter );
		printf ( "arrayfile_t::closefile: rwmode: %d\n", rwmode );
	}

	if ( narraycounter>=0 && narrays==-1 && ( this->rwmode==WRITE || this->rwmode==READWRITE ) && this->isbinary() ) {
		if ( verbose>=2 )
			printf ( "arrayfile_t: closing binary file, updating numbers %d->%d\n", narrays, narraycounter );
		if ( verbose>=3 )
			printf ( "arrayfile_t: nfid %ld\n", long ( nfid ) );
		if ( nfid != 0 ) {
			long pos = ftell ( nfid );
			//printfd("seek to from %ld to 4\n", pos);
			int r = fseek ( nfid, 4*sizeof ( int32_t ), SEEK_SET );
			//printfd("seek result %d\n", r);
			r = this->afwrite ( &narraycounter, sizeof ( int32_t ), 1 );
			if ( verbose>=2 )
				printf ( "   arrayfile_t::closefile: result of write: %d\n", ( int ) r );
			fseek ( nfid, pos, SEEK_SET ); // place back pointer
		}

	}


	// close file handles
	if ( this->nfid!=0 ) {
		fclose ( this->nfid );
		this->nfid=0;
	}
#ifdef USEZLIB
	if ( this->gzfid!=0 ) {
		gzclose ( this->gzfid );
		this->gzfid=0;
	}
#endif
//delete afile;
}

arrayfile_t::~arrayfile_t()
{
#ifdef SWIG
	swigcheck();
#endif

	if ( verbose>=2 )
		printf ( "arrayfile_t: destructor: filename %s, nfid %ld, narraycounter %d, this->rwmode %d\n", filename.c_str(), ( long ) nfid, narraycounter, this->rwmode );


	closefile();
	/*
	// close file handles
	if ( nfid!=0 ) {
	    fclose ( nfid );
	    nfid=0;
	}
	#ifdef USEZLIB
	if ( gzfid!=0 ) {
	    gzclose ( gzfid );
	    gzfid=0;
	}
	#endif
	*/
}

#ifdef HAVE_BOOST
#include "boost/filesystem.hpp"
#endif

/// return true if the specified file exists
bool file_exists ( const char *filename )
{
#ifdef HAVE_BOOST
	return boost::filesystem::exists ( filename );
#else
	FILE *fid = fopen ( filename, "r" );
	if ( fid==0 )
		return false;
	else {
		fclose ( fid );
		return true;
	}
#endif
}
// arrayfile_t::arrayfile_t(const char *fname) {
// 	const std::string s = fname;
// 	this->arrayfile_t(s);
// }

arrayfile_t::arrayfile_t ()
{
	this->verbose=0;

	//this->verbose=2;
	//printf("arrayfile_t::arrayfile_t ()\n");
	this->nfid=0;
	this->filename="";
#ifdef USEZLIB
	this->gzfid=0;
#endif

	this->mode=ATEXT;

	this->rwmode=READ;

	this->nrows=0;
	this->ncols=0;
	this->narrays=0;
	this->narraycounter=0;

}

void arrayfile_t::createfile ( const std::string fname, int nrows, int ncols, int narrays, arrayfilemode_t m , int nb )
{
	this->closefile();

	this->verbose=0;
	this->filename=fname;

	this->iscompressed=0;
#ifdef USEZLIB
	this->gzfid=0;
#endif
	this->nrows=nrows;
	this->ncols=ncols;
	this->narrays=narrays;
	this->narraycounter=0;
	this->mode=m;
	this->rwmode=WRITE;
	this->nbits = nb;	// only used in binary mode
#ifdef OADEBUG2
	printf ( "arrayfile_t: opening %s in mode %d, nb %d\n", fname.c_str(), this->mode, nb );
#endif

	if ( strcmp ( fname.c_str(), "<standard output>" ) ==0 ) {
		this->nfid = stdout;
	} else {
		this->nfid = fopen ( fname.c_str(), "w+b" );
		if ( this->nfid==0 ) {
			printf ( "arrayfile_t: ERROR: problem opening %s\n", fname.c_str() );
			return;
		}
	}
	writeheader();

}

arrayfile_t::arrayfile_t ( const std::string fnamein, int verbose )
{

	int warngz = verbose>=1;
#ifdef SWIG
	swigcheck();
#endif

	narraycounter=-1;	// make invalid

	this->verbose=verbose;

	std::string fname = fnamein;
	this->rwmode=arrayfile::READ;
	this->filename = fname;
	this->narrays=-1;	// initialize to -1

	if ( verbose>=3 )
		printfd ( "start\n" );

	std::string gzname = fname+".gz";
	//printf("arrayfile_t::arrayfile_t: %s -> %s\n", fname.c_str(), gzname.c_str() );
	if ( ! file_exists ( fname.c_str() ) && file_exists ( gzname.c_str() ) ) {
		if ( verbose && warngz ) {
			printf ( "  file %s does not exist, but gzipped file does\n", fname.c_str() );
		}
		this->filename = gzname;
		fname=gzname;
	}

#ifdef USEZLIB
	// simple check for compressed files
	if ( verbose>=2 )
		printfd ( "zlib file %s\n", fname.c_str() );
	if ( fname.substr ( fname.find_last_of ( "." ) + 1 ) == "gz" ) {
		this->iscompressed=1;
		this->gzfid = gzopen ( fname.c_str(), "rb" );
		this->nfid=0;
		if ( verbose>=2 )
			printf ( "   opened file |%s| in compressed mode: nfid %ld, gzfid %ld\n", fname.c_str(), ( long ) this->nfid, ( long ) this->gzfid );
	} else {
		this->iscompressed=0;
		this->nfid = fopen ( fname.c_str(), "r+b" );
		//printfd(" opened file: nfid %d (mode r+b)\n", this->nfid);
		if ( 0 ) {
			fseek ( nfid, 0, SEEK_SET );
			char buf[4]= {0,1,3,4};
			int r = fwrite ( buf, 1, 4, nfid );
			printfd ( "written %d bytes...\n", r );
		}
		this->gzfid=0;
		if ( verbose>=2 )
			printf ( "   opened file |%s|: nfid %ld, gzfid %ld\n", fname.c_str(), ( long ) this->nfid, ( long ) this->gzfid );
	}
#else
	this->iscompressed=0;
	this->nfid = fopen ( fname.c_str(), "r+b" );
	// printfd(" opened file: nfid %d (mode r+b)\n", this->nfid);
	this->gzfid=0;
#endif

	if ( verbose>=2 )
		printf ( " nfid %ld, gzfid %ld, isopen %d \n", ( long ) nfid, ( long ) gzfid, this->isopen() );

	if ( ! this->isopen() ) {
		if ( verbose ) {
			printf ( "problem opening file %s (iscompressed %d)\n", fname.c_str(), this->iscompressed );
		}
		this->closefile();
		return;
	}
	int magic;
	int result = afread ( &magic,1,sizeof ( int32_t ) );
	if ( result==0 ) {
		this->closefile();
		return;
	}
	if ( this->nfid ) {
		fclose ( this->nfid );
		this->nfid=0;
	}
#ifdef USEZLIB
	if ( this->gzfid ) {
		gzclose ( this->gzfid );
		this->gzfid=0;
	}
#endif

	if ( verbose>=2 )
		printf ( "arrayfile_t: reading array file %s, magic %d\n", fname.c_str(), magic );


// read the header
	if ( magic==65 ) {
		// we have a file in binary format
#ifdef USEZLIB
		if ( iscompressed ) {
			this->gzfid = gzopen ( fname.c_str(), "rb" );
		} else {
			this->nfid = fopen ( fname.c_str(), "r+b" );
		}
#else
		this->nfid = fopen ( fname.c_str(), "r+b" );
#endif

		int result = afread ( &magic,1,sizeof ( int32_t ) );
		result = afread ( &this->nbits,sizeof ( int32_t ),1 );
		assert ( result==1 );
		result = afread ( &this->nrows,sizeof ( int32_t ),1 );
		assert ( result==1 );
		result = afread ( &this->ncols,sizeof ( int32_t ),1 );
		assert ( result==1 );
		result = afread ( &this->narrays,sizeof ( int32_t ),1 );

		int binary_mode;

		int reserved;
		result = afread ( &binary_mode,sizeof ( int32_t ),1 );
		assert ( result==1 );
		result = afread ( &reserved,sizeof ( int32_t ),1 );
		result = afread ( &reserved,sizeof ( int32_t ),1 );

		//printf("arrayfile_t: constructor: binary_mode %d\n", binary_mode);

		switch ( binary_mode ) {
		case 1001:
			this->mode=arrayfile::ABINARY;
			break;
		case 1002:
			this->mode=arrayfile::ABINARY_DIFF;
			break;
		case 1003:
			this->mode=arrayfile::ABINARY_DIFFZERO;
			break;
		case 0:
			if ( verbose ) {
				printf ( "  arrayfile_t::arrayfile_t: legacy file format, file %s!\n", this->filename.c_str() );
			}
			this->mode=arrayfile::ABINARY;
			break;
		default
				:
			this->mode=arrayfile::AERROR;
			fprintf ( stderr,  "arrayfile_t::arrayfile_t : error opening binary file (binary_mode %d)\n", binary_mode );
			break;
		}

		if ( this->mode == arrayfile::ABINARY_DIFFZERO && this->nbits!=1 ) {
			printf ( "arrayfile_t: error for mode ABINARY_DIFFZERO we need 1 bit: bits %d\n", this->nbits );
			this->mode=AERROR;
			this->closefile();
			return;
		}

		if ( result!=1 )
			fprintf ( stderr,  "open binary file: wrong count in afread! %d\n", result );
	} else {
		if ( iscompressed ) {
			if ( verbose ) {
				printf ( "   compressed file: file or corrupt or gzipped text file, cannot read compressed files in TEXT mode..\n" );
			}
			this->mode=AERROR;
			this->closefile();
			return;
		} else {
			this->nfid = fopen ( fname.c_str(), "rb" );
			this->gzfid=0;
		}

		this->mode=arrayfile::ATEXT;

		char buf[1];
		buf[0]=-1;
		int r = fread ( buf, sizeof ( char ), 1, this->nfid );
		if ( buf[0] < 48 || buf[0] > 57 || r<0 ) {
			// printf("   read char %d\n", int(buf[0]));
			if ( verbose>=1 )
				fprintf ( stderr,  "   problem opening file %s (iscompressed %d)\n", fname.c_str(), this->iscompressed );
			this->closefile();
			return;

		}
		r = fseek ( this->nfid,0, SEEK_SET );

		r = fscanf ( this->nfid, "%i %i %i\n", &this->ncols, &this->nrows, &this->narrays );
		this->nbits = 0;
		//printf("arrayfile_t: text mode\n");
		if ( verbose>=2 ) {
			printf ( "arrayfile_t: text mode: header %d %d %d\n",  this->ncols, this->nrows, this->narrays );
		}
		if ( this->nrows<0 || this->nrows>20000 || this->ncols<0 || this->ncols>10000 || this->narrays<0 ) {
			if ( verbose>=1 )
				fprintf ( stderr, "   problem opening file %s (iscompressed %d)\n", fname.c_str(), this->iscompressed );
			this->nfid=0;
			this->gzfid=0;
		}


	}

	narraycounter=0;

	if ( verbose>=2 )
		printf ( "   arrayfile_t::arrayfile_t: mode %d, nrows %d, ncols %d, narrays %d\n", mode, this->nrows, this->ncols, this->narrays );

}


/**
 * @brief Write an array in binary mode to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::write_array_binary_diff ( const array_link &A )
{
	myassertdebug ( this->rwmode == WRITE, "error: arrayfile_t not in write mode" );

	myassertdebug ( A.maxelement() <= 1, "arrayfile_t::write_array_binary_diff: array is not binary" );

	int ngood=0;

	rowindex_t N = this->nrows;
	const int num=N*sizeof ( array_t );

	if ( 0 ) {
		diffarray.show();
		diffarray.showarray();
		printf ( "-------------\n" );
	}
	if ( 0 ) {
		A.show();
		A.showarray();
		printf ( "-------------\n" );
	}

//    printf("write_array_binary_diff: i %d\n", i);

	for ( int i=0; i<diffarray.n_columns; i++ ) {
// printf("write_array_binary_diff: i %d: %d\n", i, memcmp ( this->diffarray.array+N*i, A.array+N*i, num));
		if ( ! memcmp ( this->diffarray.array+N*i, A.array+N*i, num ) ) {
			ngood++;

		} else
			break;
	}

	int32_t nwrite=A.n_columns-ngood;
	//printf("write_array_binary_diff: %d good columns, writing %d to disk\n", ngood, nwrite);
	array_link rest = A.selectLastColumns ( nwrite );

	int n = fwrite ( ( const void* ) &nwrite,sizeof ( int32_t ),1, nfid );

	this->write_array_binary ( rest );

// update with previous array
	this->diffarray = A;
}

/**
 * @brief Write an array in binary mode to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::write_array_binary_diffzero ( const array_link &A )
{
	myassertdebug ( this->rwmode == WRITE, "error: arrayfile_t not in write mode" );
	int ngood=0;

	rowindex_t N = this->nrows;
	const int num=N*sizeof ( array_t );

	for ( int i=0; i<diffarray.n_columns; i++ ) {
//printf("write_array_binary_diffzero: i %d: %d\n", i, memcmp ( this->diffarray.array+N*i, A.array+N*i, num));
		if ( ! memcmp ( this->diffarray.array+N*i, A.array+N*i, num ) ) {
			ngood++;

		} else
			break;
	}

	int16_t nwrite=A.n_columns-ngood;
	//printf("write_array_binary_diffzero: nwrite %d, ngood %d\n", nwrite, ngood);
	//A.show();
	array_link rest = A.selectLastColumns ( nwrite );

	int n = fwrite ( ( const void* ) &nwrite,sizeof ( int16_t ),1, nfid );

	//printf("diffarray:\n"); this->diffarray.show();

	if ( diffarray.n_columns>0 ) {

		array_link diffrest = this->diffarray.selectLastColumns ( nwrite );

		array_link z = rest-diffrest;
		for ( int i=0; i<z.n_columns*z.n_rows; i++ )
			z.array[i] = ( 2+z.array[i] ) % 2;

//   printf("write_array_binary_diff: %d good columns, writing %d to disk\n", ngood, nwrite);
		//  printf ( "diffzero: writing rest: fraction non-zero %.2f\n", z.nonzero_fraction() );
		//z.transposed().showarray();


		//printf("write_array_binary_diffzero: writing array: "); z.show();
		this->write_array_binary ( z );
	} else {
		array_link z = rest;
		//printf("write_array_binary_diffzero: writing array (i): "); z.show();
		this->write_array_binary ( z );
	}

// update with previous array
	this->diffarray = A;
}

/// Write an array in binary mode to a file
void arrayfile_t::write_array_binary ( const array_link &A )
{
	write_array_binary ( A.array, A.n_rows, A.n_columns );
}


/**
 * @brief Write an array in binary mode to a file
 * @param fid
 * @param array
 * @param nrows
 * @param ncols
 */
void arrayfile_t::write_array_binary ( carray_t *array, const int nrows, const int ncols )
{
#ifdef OADEBUG
	int m=0;
	for ( int i=0; i<nrows*ncols; ++i )
		if ( array[i]>m )
			m=array[i];
	if ( this->nbits==1 ) {
		if ( m>1 )
			printf ( "ERRROR!\n" );
	}
#endif

//printf("arrayfile_t::write_array_binary: nbits %d\n", this->nbits);

	// TODO: the type of array should be taken into account?
	if ( sizeof ( array_t ) ==sizeof ( int32_t ) && this->nbits==32 )
		int n = fwrite ( ( const void* ) array,sizeof ( int32_t ),nrows*ncols,nfid );
	else {
		switch ( this->nbits ) {
		case 8:
			writeblob<array_t, char> ( array, nrows*ncols,nfid );
			break;
		case 32:
			writeblob<array_t, int32_t> ( array, nrows*ncols, nfid );
			break;
		case 1: {
			//printf("writing array in binary mode\n");

			// construct bit array
			BIT_ARRAY* bitarr = bit_array_create ( nrows*ncols );
			// fill bit array
			for ( int i=0; i<nrows*ncols; i++ ) {
				if ( array[i] )
					bit_array_set_bit ( bitarr, i );
				else
					bit_array_clear_bit ( bitarr, i );
			}
			word_addr_t num_of_words = nwords ( nrows*ncols ); //printf("num_of_words: %d\n", (int)num_of_words);

			fwrite ( bitarr->words,num_of_words,sizeof ( word_t ), this->nfid );

#ifdef OADEBUG2
			printf ( "1-bit write: %ld bytes for array with %d elements\n", num_of_words*sizeof ( word_t ), nrows*ncols );
#endif
			bit_array_free ( bitarr );
		}

		break;
		default
				:
			printf ( "error: number of bits not supported\n" );
			break;
		}
	}
}




arrayfile_t::arrayfile_t ( const std::string fname, int nrows, int ncols, int narrays, arrayfilemode_t m, int nb )
{
	//closefile();
	this->verbose=0;
	this->nfid=0;
	this->narrays=-1;
	this->narraycounter=-1;
	this->rwmode=READ;
	this->mode=ATEXT;

#ifdef USEZLIB
	this->gzfid=0;
#endif
	createfile ( fname, nrows, ncols, narrays, m, nb );

}

size_t arrayfile_t::afread ( void * ptr, size_t sz, size_t cnt )
{
#ifdef USEZLIB
	size_t r;
	if ( iscompressed ) {
		assert ( this->gzfid );
		r=gzread ( this->gzfid, ptr, sz*cnt );
		r/=sz;
	} else {
		r=fread ( ptr, sz, cnt, this->nfid );
	}
	//printf( "r: %d,, n %d\n", r, n);
	if ( r==0 ) {
		//printf("could not read from file");
	}
#else
	size_t r;
	r=fread ( ptr, sz, cnt, this->nfid );
#endif
	//assert(r==n);
	return r;
}


int arrayfile_t::seek ( int pos )
{
	if ( mode!= ABINARY ) {
		printf ( "error: seek only possible for binary files\n" );
		return -1;
	}

	if ( pos<0 || pos>1000000000 ) {
		printf ( "error: position specified is invalid\n" );
		return -1;
	}
	if ( nfid )
		fseek ( nfid, this->headersize() + this->barraysize() *pos, SEEK_SET );
#ifdef USEZLIB
	if ( gzfid )
		gzseek ( gzfid, this->headersize() + this->barraysize() *pos, SEEK_SET );
#endif

	narraycounter = pos;
	return pos;

}


int arrayfile_t::read_array_binary_zero ( array_link &a )
{
	const int dbg=0;

	// no index is written or read (to save disk space
	a.index=array_link::INDEX_NONE;
	int index=array_link::INDEX_NONE;

	if ( dbg )
		printf ( " arrayfile_t::read_array_binary_zero: gp\n" );

	if ( a.n_columns!=diffarray.n_columns && diffarray.n_columns!=-1 && diffarray.n_columns!=0 ) {
		fprintf ( stderr, "arrayfile_t::read_array_binary_zero: error different number of columns %d %d\n", diffarray.n_columns, a.n_columns );
		return array_link::INDEX_ERROR;
	}

	if ( dbg )
		printf ( " arrayfile_t::read_array_binary_zero: reading nrest\n" );

	int16_t nrest;
	int result = afread ( &nrest, sizeof ( int16_t ), 1 );
	if ( result!=1 ) {
		// error reading array, we could have reached the end of the file
		if ( this->narrays==-1 ) {
		} else {
			fprintf ( stderr, "arrayfile_t::read_array_binary_zero: error reading array: index %d, result %d\n", index, result );
		}
		index=array_link::INDEX_ERROR;
		return index;
	}

	int nrows = a.n_rows;
	int ncols=a.n_columns;
	int ngood=ncols-nrest;
	if ( dbg ) {
		diffarray.show();
		a.show();
		printf ( " arrayfile_t::read_array_binary_zero: xxx nrows %d, ngood %d\n", nrows, ngood );
	}
	if ( diffarray.n_columns>0 ) {
		copy_array ( diffarray.array, a.array, nrows, ngood );
	} else {
		// diffarray not initialized yet...
	}
	if ( dbg )
		printf ( "arrayfile_t::read_array: nrows %d, ncols %d,  nrest %d, ngood %d\n", nrows, ncols, nrest, ngood );
	this->read_array_binary ( a.array+nrows*ngood, nrows, nrest );

	if ( dbg ) {
		printf ( "arrayfile_t::read_array:read_array_binary done: %d %d\n", a.n_columns, diffarray.n_columns );
		a.showarray();
	}
	if ( a.n_columns==diffarray.n_columns ) {
		// update array
		if ( dbg )
			printf ( "  here: %d\n", a.n_columns );

		int N = a.n_rows;
		for ( int i=0; i<nrest; i++ ) {
			for ( int j=0; j<N; j++ ) {
				int idx=j+ ( i+ngood ) *N;
				a.array[idx] += diffarray.array[idx];
				a.array[idx] = a.array[idx] % 2;
			}
		}
		if ( dbg )
			printf ( "  here: %ld %ld\n", ( long ) diffarray.array, ( long ) a.array );

	}
	diffarray=a;
	if ( dbg )
		printf ( "  arrayfile_t::read_array:read_array_binary_zero: xxx\n" );
	return index;
}

int arrayfile_t::read_array ( array_link &a )
{
	int32_t index;
	switch ( this->mode ) {
	case arrayfile::ATEXT:
		index = this->read_array ( a.array, a.n_rows, a.n_columns );
		a.index = index;
		break;
	case arrayfile::ABINARY:
		index = this->read_array ( a.array, a.n_rows, a.n_columns );
		a.index = index;
		break;
	case arrayfile::ABINARY_DIFFZERO:
		index = this->read_array_binary_zero ( a );
		a.index = index;
		break;
	case arrayfile::ABINARY_DIFF: {
		int result = afread ( &index, sizeof ( int32_t ), 1 );
		//printf("   arrayfile::ABINARY_DIFF %d %d\n", result, index);
		if ( result!=1 ) {
			// error reading array, we could have reached the end of the file
			if ( this->narrays==-1 ) {
			} else {
				fprintf ( stderr,  "error reading array: index %d, result %d\n", index, result );
			}
			index=-1;
			a.index=index;
			break;

		}
		if ( index==-1 ) {
			printf ( "arrayfile_t::read_array: updating index (-1 is invalid)\n" );
			index=0;
		}
		int32_t nrest;
		result = afread ( &nrest, 1,sizeof ( int32_t ) );

		int nrows = a.n_rows;
		int ncols=a.n_columns;
		int ngood=ncols-nrest;
		copy_array ( diffarray.array, a.array, nrows, ngood );
		//a=diffarray;

		if ( 0 )
			printf ( "arrayfile_t::read_array: nrows %d, ncols %d,  nrest %d, ngood %d\n", nrows, ncols, nrest, ngood );
		this->read_array_binary ( a.array+nrows*ngood, nrows, nrest );
		a.index=index;
		diffarray=a;
		//printf("hack index: %d %d\n", index, diffarray.index);
	}
	break;
	default
			:
		index=-1;
		break;
	}

	return index;
}

int arrayfile_t::read_array ( array_t* array, const int nrows, const int ncols )
{
	int index=-10;
	if ( nrows!=this->nrows )
		printf ( "arrayfile_t::read_array: nrows unequal %d %d\n", nrows, this->nrows );
	if ( ncols!=this->ncols )
		printf ( "arrayfile_t::read_array: ncols unequal\n" );

	//printf("arrayfile_t::read_array: mode %d\n", this->mode);

	switch ( this->mode ) {
	case arrayfile::ATEXT: {
		int r = fscanf ( nfid, "%d\n", &index );
		//printf("index %d\n", index);
		::read_array ( nfid, array, nrows, ncols );
		break;
	}
	case arrayfile::ABINARY: {
		int result = afread ( &index, sizeof ( int32_t ), 1 );
		if ( result!=1 ) {
			index=-1;
			printf ( "arrayfile_t::read_array: error: could not read array index: result %d, index %d\n", result, index );
			break;
		}
		if ( index==-1 ) {
			printf ( "updating index (-1 is invalid)\n" );
			index=0;
		}

		//printf("arrayfile_t::read_array: index 2 %d\n", index);
		this->read_array_binary ( array, nrows, ncols );
	}
	break;
	default
			:
		fprintf ( stderr,  "arrayfile_t::read_array: error: no such mode %d\n", this->mode );
		break;
	}

	return index;
}


void arrayfile_t::writeheader()
{
	assert ( this->nfid );

	if ( this->mode == arrayfile::ABINARY || this->mode == arrayfile::ABINARY_DIFF || this->mode == arrayfile::ABINARY_DIFFZERO ) {
		//printf("  arrayfile_t::writeheader(): binary mode\n");
		int magic=65;
		int32_t reserved=0;

		if ( this->mode == arrayfile::ABINARY_DIFFZERO ) {
			// TODO: why this check?
			assert ( this->nbits==1 );
		}
		fwrite ( ( const void* ) & magic,sizeof ( int ),1,this->nfid );
		fwrite ( ( const void* ) & this->nbits,sizeof ( int ),1,this->nfid );
		fwrite ( ( const void* ) & this->nrows,sizeof ( int ),1,this->nfid );
		fwrite ( ( const void* ) & this->ncols,sizeof ( int ),1,this->nfid );
		fwrite ( ( const void* ) & this->narrays,sizeof ( int ),1,this->nfid );

		if ( this->mode == arrayfile::ABINARY )


			reserved=1001;
		else {
			if ( this->mode == arrayfile::ABINARY_DIFF )
				reserved=1002;
			else
				reserved=1003;
		}
		fwrite ( ( const void* ) & reserved,sizeof ( int ),1,this->nfid );
		reserved=9999;
		fwrite ( ( const void* ) & reserved,sizeof ( int ),1,this->nfid );
		fwrite ( ( const void* ) & reserved,sizeof ( int ),1,this->nfid );
	} else {
		fprintf ( this->nfid, "%i %i %i\n", this->ncols, this->nrows, this->narrays );
	}
}

int nArrays ( const char *fname )
{
	arrayfile_t af ( fname, 0 );
	int n = af.narrays;
	af.closefile();
	return n;
}

arraylist_t  readarrayfile ( const char *fname, int verbose, int *setcols )
{
	arraylist_t v;
	readarrayfile ( fname, &v, verbose, setcols );
	return v;
}

/**
 * @brief Read all arrays in a file and store then in an array list
 * @param fname
 * @param arraylist
 * @return
 */
int readarrayfile ( const char *fname, arraylist_t * arraylist, int verbose, colindex_t *setcols, rowindex_t *setrows, int *setbits )
{
	//verbose=3;
	if ( arraylist==0 )
		arraylist = new arraylist_t;

	arrayfile_t *afile = new arrayfile_t ( fname, verbose );
	if ( setcols!=0 )
		( *setcols ) = afile->ncols;
	if ( setrows!=0 )
		( *setrows ) = afile->nrows;
	if ( setbits!=0 )
		( *setbits ) = afile->nbits;

	//if ( verbose )
	//   log_print ( NORMAL, "init_restart: number of columns: %i, number of rows: %i, number of arrays: %i\n", afile->ncols, afile->nrows, afile->narrays );

	if ( ! afile->isopen() ) {
		if ( verbose ) {
			printf ( "readarrayfile: problem with file %s\n", fname );
			printf ( " readarrayfile %d\n", ( int ) arraylist->size() );
		}

		return 0;
	}

	int narrays = afile->narrays;
	if ( afile->narrays<0 )
		narrays = arrayfile_t::NARRAYS_MAX;

	array_link *alink;
	long i;
	int index;
	for ( i = 0; i < narrays; i++ ) {

		if ( i%10000==0 && verbose ) {
			log_print ( QUIET, "readarrayfile: loading arrays: %d/%d\n", i, afile->narrays );
		}

		if ( verbose>=3 )
			printf ( "readarrayfile: i %ld\n", i );
		alink = new array_link ( afile->nrows, afile->ncols, i+1 );
		//index=afile->read_array ( alink->array, afile->nrows, afile->ncols );
		index=afile->read_array ( *alink );
		// printf("hack %d\n", index);

		// NOTE: for array files we have used -1 to read arrays to the end of file,
		if ( index<0 ) {
			printf ( "readarrayfile: index %d, problem?\n", index );
			delete alink;
			break;

		}
		if ( verbose>=4 ) {
			alink->showarray();
		}
		arraylist->push_back ( *alink );
		delete alink;
	}

	delete afile;
	return i;
}


/**
 * @brief Write all arrays in a list to file
 * @param fname
 * @param arraylist
 * @param mode
 * @return
 */
int writearrayfile ( const char *fname, const arraylist_t *arraylist, arrayfile::arrayfilemode_t mode, int nrows, int ncols )
{
	int nb=8;	// default: char

	if ( arraylist->size() ==0 ) {
		if ( mode==arrayfile::ABINARY_DIFFZERO ) {
			// special case: diffzero should always have 1 bit files
			nb=1;
		}
		if ( nrows<=0 )
			printf ( "writearrayfile: warning: empty list, using nrows %d, ncols %d\n", nrows, ncols );
	} else {
		nrows=arraylist->at ( 0 ).n_rows;
		ncols = arraylist->at ( 0 ).n_columns;

		nb = arrayfile_t::arrayNbits ( arraylist->at ( 0 ) );
	}

	arrayfile_t *afile = new arrayfile_t ( fname, nrows, ncols, arraylist->size(), mode, nb );

	if ( ! afile->isopen() ) {
		printf ( "writearrayfile: problem with file %s\n", fname );
		return 0;
	}

	int i = afile->append_arrays ( *arraylist,1 ); // append_arrays ( afile, *arraylist, 1 );
	afile->finisharrayfile();
	delete afile;

	return i;
}

#include <errno.h>

int writearrayfile ( const char *fname, const array_link &al, arrayfile::arrayfilemode_t mode )
{
	arraylist_t s;
	s.push_back ( al );
	return writearrayfile ( fname, &s, mode );
}

/// append a single array to an array file. creates a new file if no file exists
int appendarrayfile ( const char *fname, const array_link al )
{
	int dverbose=0;
	int nb=8;	// default: char
	int nrows=-1, ncols=-1;
	arrayfilemode_t mode  = arrayfile::ABINARY;

	arraylist_t arraylist;
	arraylist.push_back ( al );

	arrayfile_t *afile = new arrayfile_t ( fname );

	if ( dverbose )
		printf ( "\n### appendarrayfile: opened array file: %s\n", afile->showstr().c_str() );




	if ( ! afile->isopen() ) {
		printf ( "appendarrayfile: creating new array file %s\n", fname );


		{
			nrows=al.n_rows; //
			ncols=al.n_columns;

			nb = arrayfile_t::arrayNbits ( al );
		}

		afile = new arrayfile_t ( fname, nrows, ncols, -1, mode, nb );
	} else {
		if ( dverbose ) {
			printfd ( "file is at position %d\n", ftell ( afile->nfid ) );
		}
		if ( 0 ) {
			//afile->nfid = fopen(fname, "r+b");

			fseek ( afile->nfid, 0, SEEK_SET );
			printf ( "  " );
			printfd ( "file is at position %d\n", ftell ( afile->nfid ) );

			char buf[4]= {1,2,3,4};
			int p = fseek ( afile->nfid, 0, SEEK_END );
			printf ( "  p %d, afile->nfid %ld\n", p, ( long ) afile->nfid );
			printf ( " fseek errno: %s\n", strerror ( errno ) );
			int rr=fread ( buf, 1, 4, afile->nfid );
			printf ( "rr %d, something went wrong with fread()? %s\n", rr, strerror ( errno ) );
			int pp = ftell ( afile->nfid );
			fseek ( afile->nfid, 0, SEEK_CUR );
			rr=fwrite ( buf, 1, 4, afile->nfid );
			printf ( "rr %d, something went wrong with fwrite()! %s\n", rr, strerror ( errno ) );
			printf ( "  written rr %d\n", rr );
			afile->seek ( afile->narrays );
			exit ( 0 );
		}
	}

	if ( 0 ) {
		printfd ( "hack...\n" );
		fseek ( afile->nfid, 0, SEEK_SET );
		char buf[4]= {0,1,3,4};
		int r = fwrite ( buf, 1, 4, afile->nfid );
		printfd ( "written %d bytes...\n", r );
		exit ( 0 );
	}

	//afile->setVerbose(3);
	//printf("afile state: %s\n", afile->showstr().c_str() );

	//printf("before append: afile->narrays %d, narraycounter %d\n", afile->narrays, afile->narraycounter);
	afile->seek ( afile->narrays );
	//printf("after seek: afile->narrays %d, narraycounter %d, ftell %ld\n", afile->narrays, afile->narraycounter, ftell(afile->nfid) );

	//int IsReadOnly = afile->nfid->_flag;

#ifdef WIN32
#else
	// debugging code
	if ( 0 ) {
		int fs = get_file_status ( afile->nfid );
		printf ( "# file status %d (O_RDONLY %d, O_RDWR %d, O_APPEND %d)\n", fs, O_RDONLY , O_RDWR, O_APPEND );
		printf ( "# file status mask: (O_RDONLY %d, O_WRONLY %d, O_RDWR %d, O_APPEND %d)\n", O_RDONLY & fs , O_WRONLY & fs, O_RDWR & fs, O_APPEND & fs );
		printf ( "# file status mod 4: %d\n", fs % 4 );
	}
#endif


	if ( 0 ) {
		char buf[4]= {1,2,3,4};
		int p = fseek ( afile->nfid, 0, SEEK_SET );
		printf ( "  p %d\n", p );
		int rr=fwrite ( buf, 1, 4, afile->nfid );
		printf ( "  written rr %d\n", rr );
		afile->seek ( afile->narrays );
	}

	if ( ! afile->isopen() ) {
		printf ( "writearrayfile: problem with file %s\n", fname );
		return 0;
	}

	int i = afile->append_arrays ( arraylist,1 ); // append_arrays ( afile, *arraylist, 1 );
	//printf("after append: afile->narrays %d, narraycounter %d\n", afile->narrays, afile->narraycounter);
	//printf("afile state: %s\n", afile->showstr().c_str() );
	afile->narrays=-1;
	afile->rwmode=READWRITE;

	afile->finisharrayfile();
	delete afile;

	return i;
}

void  selectArrays ( const std::string filename,   std::vector<int> &idx, arraylist_t &rl, int verbose )
{
	arrayfile_t af ( filename, 0 );
	array_link al ( af.nrows, af.ncols, -1 );
	if ( af.mode==ABINARY ) {
		for ( std::vector<int>::iterator it = idx.begin(); it<idx.end(); ++it ) {
			if ( verbose )
				printf ( "selectArrays: idx %d\n", *it );
			af.seek ( *it );
			af.read_array ( al );
			rl.push_back ( al );
		}
	} else {
		// check whether list is sorted
		indexsort vv ( idx );
		//vv.sortdescending(idx);

		if ( vv.issorted() ) {
			int cpos=0;
			for ( std::vector<int>::iterator it = idx.begin(); it<idx.end(); ++it ) {
				int pos =  *it;
				if ( verbose )
					printf ( "selectArrays: idx %d, cpos %d\n", pos, cpos );
				if ( pos<cpos ) {
					printf ( "selectArrays: arrayfile in text mode and negative seek, aborting!!! %d %d\n", pos, cpos );
					// TODO: implement this
					return;
				}
				int nsk=pos-cpos;
				if ( verbose )
					printf ( "selectArrays: skipping %d arrays\n", nsk );
				for ( int j=0; j< ( nsk ); j++ )  {
					af.read_array ( al );
					if ( verbose>=3 ) {
						std::vector<double> tmp = al.GWLP();
						printf ( "  gwlp: " );
						display_vector ( tmp );
						printf ( "\n" );
					}
					cpos++;
				}
				af.read_array ( al );
				cpos++;
				rl.push_back ( al );
			}
		} else {
			if ( verbose>=2 )
				printf ( "selectArrays: text file with unsorted indices, UNTESTED CODE!\n" );
			if ( verbose )
				printf ( "selectArrays: no sorted indices!\n" );
			if ( verbose>=2 ) {
				cout << "idx: ";
				display_vector<int> ( idx );
				std::cout << std::endl;
			}
			// not sorted!
			int nn=vv.indices.size();
			std::vector<int> sidx = vv.sorted ( idx );
			if ( verbose>=2 ) {
				std::cout << "sidx: ";
				display_vector<int> ( sidx );
				std::cout << std::endl;
			}
			//vv.sortdescending(idx);
			// std::vector<int> sidx2 = vv.sorted(idx);
			//cout << "sidx2: "; display_vector<int>(sidx2); cout << endl;
			arraylist_t tmp;
			if ( verbose>=2 ) {
				std::cout << "vv.indices: ";
				display_vector<int> ( vv.indices );
				std::cout << std::endl;
			}
			selectArrays ( filename,  sidx, tmp, 0 );
			if ( ( int ) tmp.size() !=nn ) {
				printf ( "ERROR!\n" );
				return;
			}
			for ( int j=0; j<nn; j++ ) {
				int x=vv.indices[j];
				// TODO: or reverse entry?
				rl.push_back ( tmp[x] );
			}
		}
	}
}

array_link selectArrays ( const std::string filename, int ii )
{
	arrayfile_t af ( filename );
	array_link al ( af.nrows, af.ncols, -1 );
	if ( ii<0 ) {
		fprintf ( stderr,  "selectArrays: error: negative index\n" );
		return al;
	}
	if ( af.mode==ABINARY ) {
		af.seek ( ii );
		af.read_array ( al );
	} else {
		for ( int i=0; i<ii; i++ )
			af.read_array ( al );
		af.read_array ( al );
	}
	return al;
}

void  selectArrays ( const arraylist_t &al,   std::vector<int> &idx, arraylist_t &rl )
{
	for ( std::vector<int>::iterator it = idx.begin(); it<idx.end(); it++ ) {
		rl.push_back ( al.at ( *it ) );
	}
}
void  selectArrays ( const arraylist_t &al,   std::vector<long> &idx, arraylist_t &rl )
{
	for ( std::vector<long>::iterator it = idx.begin(); it<idx.end(); it++ ) {
		if ( ( *it ) >=0 && ( ( *it ) < ( int ) al.size() ) )
			rl.push_back ( al.at ( *it ) );
		else {
			printf ( "selectArrays: index %ld out of bounds!\n", *it );
		}
	}
}

arraylist_t  selectArrays ( const arraylist_t &al,   std::vector<int> &idx )
{
	arraylist_t rl;
	for ( std::vector<int>::iterator it = idx.begin(); it<idx.end(); ++it ) {
		int val = *it;
		if ( val>=0 && val<= ( int ) al.size() )
			rl.push_back ( al.at ( val ) );
		else
			fprintf ( stderr,  "selectArrays: error: index out of bounds: index %d, size %zu\n", val, al.size() );
	}
	return rl;
}

arraylist_t  selectArrays ( const arraylist_t &al,   std::vector<long> &idx )
{
	arraylist_t rl;
	for ( std::vector<long>::iterator it = idx.begin(); it<idx.end(); ++it ) {
		int val = *it;
		if ( val>=0 && val<= ( int ) al.size() )
			rl.push_back ( al.at ( val ) );
		else
			fprintf ( stderr,  "selectArrays: error: index out of bounds: index %d, size %zu\n", val, al.size() );
	}
	return rl;
}

array_link array_link::operator+ ( const array_link &b ) const
{
	assert ( this->equalsize ( b ) );
	array_link tmp = ( *this );
	for ( int i=0; i<tmp.n_columns*tmp.n_rows; i++ )
		tmp.array[i] += b.array[i];
	return tmp;
}

array_link array_link::operator- ( const array_link &b ) const
{
	assert ( this->equalsize ( b ) );
	array_link tmp = ( *this );
	for ( int i=0; i<tmp.n_columns*tmp.n_rows; i++ )
		tmp.array[i] -= b.array[i];
	return tmp;
}


/// stack to arrays together
array_link hstack ( const array_link &al, const array_link &b )
{
	assert ( al.n_rows==b.n_rows );
	array_link v ( al.n_rows, al.n_columns+b.n_columns, array_link::INDEX_NONE );
	std::copy ( al.array, al.array+al.n_columns*al.n_rows, v.array );
	std::copy ( b.array, b.array+b.n_columns*al.n_rows, v.array+v.n_rows*al.n_columns );
	return v;
}

/// append the last column of the second array to the entire first array
array_link hstacklastcol ( const array_link &al, const array_link &b )
{
	array_link v ( al.n_rows, al.n_columns+1, array_link::INDEX_NONE );
	std::copy ( al.array, al.array+al.n_columns*al.n_rows, v.array );
	size_t offset = al.n_rows* ( b.n_columns-1 );
	std::copy ( b.array+offset, b.array+offset+al.n_rows, v.array+v.n_rows*al.n_columns );
	return v;
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4; 
