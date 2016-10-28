#include <math.h>
#include <algorithm>

#include "oaoptions.h"
#include "arraytools.h"
#include "lmc.h"
#include "tools.h"
#include "mathtools.h"
#include "extend.h"

#include "nonroot.h"

#ifndef DOOPENMP
// no OpenMP
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

#ifdef DOOPENMP
#include <omp.h>
#endif


using namespace std;

// legacy global static object
LMC_static_struct_t _globalStaticX;
LMC_static_struct_t &getGlobalStaticOne()
{
#ifdef OADEBUG
//  printfd("  return _globalStaticX...\n");
#endif
    return _globalStaticX;
}

LMC_static_struct_t * getGlobalStaticOnePointer()
{
#ifdef OADEBUG
//  printfd("  return globalStaticOne\n");
#endif
    return &_globalStaticX;
}


// worker pool of static objects

static object_pool< LMC_static_struct_t > staticDataPool;


#ifdef OADEBUG
int getGlobalStaticNumber ( LMC_static_struct_t * p )
{
    return p->id; //staticDataPool.id(p);
}
#endif

LMC_static_struct_t * getGlobalStatic()
{
#ifdef NOREUSE
    LMC_static_struct_t *pp = new LMC_static_struct_t();
    pp->id = 100*omp_get_thread_num() + int ( random() % 40 );
    return pp;
#endif

#ifdef OADEBUGMORE
    staticDataPool.verbose=1;
#endif

    LMC_static_struct_t *p = 0;
// printfd("enter critical: %d\n", omp_get_thread_num() );
    #pragma omp critical
    {
        p = staticDataPool.New();
    }
    //printfd("leave critical: %d\n", omp_get_thread_num() );
    return p;

}
void releaseGlobalStatic ( LMC_static_struct_t *p )
{
#ifdef NOREUSE
    delete p;
    return;
#endif

//printfd("releaseGlobalStatic: enter critical: %d\n", omp_get_thread_num() );
    #pragma omp critical
    {
        staticDataPool.Delete ( p );
// printfd("   releaseGlobalStatic: pool size %d\n", staticDataPool.size() );
    }
}

void cleanGlobalStatic()
{
    staticDataPool.reset();
}



// indexed pool of static objects
std::vector<LMC_static_struct_t *> globalStaticPool;
LMC_static_struct_t * getGlobalStaticIndexed ( int n )
{
    const int verbose=0;
    //return &globalStatic;

    if ( verbose )
        printf ( "getGlobalStatic: n %d, pool size %ld\n", n, ( long ) globalStaticPool.size() );
    if ( n>= ( int ( globalStaticPool.size() ) ) ) {
        printf ( "  allocating new element in globalStaticPool: n %d, pool size %ld\n", n, ( long ) globalStaticPool.size() );
        size_t psize = globalStaticPool.size();
        //globalStaticPool.resize(n+1);
        for ( int jj=psize; jj<=n; jj++ ) {
            LMC_static_struct_t *p= new LMC_static_struct_t();
            if ( verbose )
                printf ( "new element jj %d\n", jj );
            globalStaticPool.push_back ( p );
        }
    }

    if ( verbose ) {
        printf ( "   " );
        globalStaticPool[n]->show();
    }
    //printf("  before return...\n");
    return globalStaticPool.at ( n );

}

void cleanGlobalStaticIndexed()
{
    printfd ( "cleanGlobalStaticIndexed: delete %d items\n", globalStaticPool.size() );
    for ( size_t jj=0; jj<globalStaticPool.size(); jj++ ) {
        //printfd("cleanGlobalStatic: delete item %d\n", jj);
        delete globalStaticPool[jj];
    }
    globalStaticPool.resize ( 0 );
}




//

LMC_static_struct_t::~LMC_static_struct_t()
{
#ifdef OADEBUG
//  printf("LMC_static_struct_t::~LMC_static_struct_t: last ref %s\n", ref.c_str() );
#endif

    // call free() function (need to extrect from init code)
    this->freeall();
}

#ifdef OADEBUG
static int LMC_static_count = 0;
#endif

LMC_static_struct_t::LMC_static_struct_t()
{
#ifdef OADEBUG
    //printf("LMC_static_struct_t::LMC_static_struct_t()\n");
//#pragma omp atomic
    //id = LMC_static_count;
    //LMC_static_count++;
#endif
    /* initialize static structure values to zero */
    this->LMC_non_root_init = 0;
    this->LMC_root_init = 0;
    this->LMC_root_rowperms_init = 0;
    this->rootrowperms_full = 0;

    this->ad = 0;
    this->colbuffer = 0;
    this->nrootrowperms = 0;
    this->rootrowperms = 0;
    this->dyndata_p = 0;
    this->current_trans = 0;
    this->colperm_p = 0;
    this->localcolperm_p = 0;
}

/// return name of the algorithm
std::string algnames ( algorithm_t m )
{
    std::string str;
    switch ( m ) {
    case MODE_ORIGINAL:
        str =  stringify ( MODE_ORIGINAL );
        break;
    case MODE_LMC_2LEVEL:
        str =  stringify ( MODE_LMC_DEBUG );
        break;
    case MODE_J4:
        str =  stringify ( MODE_J4 );
        break;
    case MODE_LMC_SYMMETRY:
        str =  stringify ( MODE_LMC_SYMMETRY );
        break;
    case MODE_AUTOSELECT:
        str =  stringify ( MODE_AUTOSELECT );
        break;
    case MODE_J5ORDER:
        str =  stringify ( MODE_J5ORDER );
        break;
    case MODE_J5ORDERX:
        str =  stringify ( MODE_J5ORDERX );
        break;
    case MODE_J5ORDERXFAST:
        str =  stringify ( MODE_J5ORDERXFAST );
        break;
    default
            :
    case MODE_INVALID:
        printfd ( "MODE_INVALID!" );
        str =  stringify ( MODE_INVALID );
        break;
    }
    //printf("%d: %s", m, str.c_str() );
    return str;
}



/**
* @brief Perform a row permutation on an array using a rowsort_t structure
* @param source
* @param target
* @param rs
* @param nrows
* @param ncols
*/
void perform_row_permutation ( const carray_t *source, array_t *target, rowsort_t *rs, rowindex_t nrows, colindex_t ncols )
{
    int t;
    for ( colindex_t i = 0; i < ncols; i++ ) {
        for ( rowindex_t j=0; j<nrows; j++ ) {
            t = rs[j].r; //perm[j];
            target[nrows*i+j]=source[nrows*i+t];
        }
    }
}





/* Debugging functions */
void debug_print_xxx ( lmc_t ret, int cpoffset, int nlevels, const dyndata_t *dyndata, levelperm_t lperm, const LMCreduction_t *reduction, carray_t *original, const arraydata_t *ad )
{
    if ( 0 && ( ret==LMC_LESS ) ) {
        //if (1 && (ret==LMC_LESS || ret==LMC_MORE)) {
        printf ( "### decision: ret %d\n", ret );
        reduction->transformation->show();
        print_array ( "original\n", original, ad->N, ad->ncols );
        // reduction->transformation->apply ( original, reduction->array );
        print_array ( "update:\n", reduction->array, ad->N, ad->ncols );

        /// Compare 2 arrays and return position of first difference
        rowindex_t rpos;
        colindex_t cpos;
        int vv = array_diff ( original+cpoffset,  reduction->array+dyndata->col*ad->N, ad->N, 1, rpos, cpos );
        printf ( "   array_diff result: dyndata->col %d, rpos %d, cpos %d, cpoffset %d, vv %d\n", dyndata->col, rpos, cpos,cpoffset, vv );

        printf ( "lperm: " );
        print_perm ( lperm, nlevels );
        print_perm ( original+cpoffset, ad->N );
        print_perm ( reduction->array+dyndata->col*ad->N, ad->N );
        dyndata->show();
    }

}

void debug_show_reduction ( LMCreduction_t *reduction, carray_t *original, const arraydata_t *ad, const char *str = "##\n" )
{
    printf ( "####3 new reduction:\n" );
    reduction->transformation->show();
    print_array ( "original:\n", original, ad->ncols, ad->N );
    reduction->transformation->apply ( original, reduction->array );
    print_array ( "update:\n", reduction->array, ad->ncols, ad->N );
}

void debug_print_transformations ( LMCreduction_t *reduction, carray_t *original, carray_t *array, const arraydata_t *ad, const dyndata_t *dyndata )
{
    printf ( "### LMCreduce_non_root: col %d\n", dyndata->col );
    printf ( "current array:\n" );

    reduction->transformation->show();
    printf ( "orignal:\n" );
    print_array ( original, ad->ncols, ad->N );
    printf ( "current:\n" );
    reduction->transformation->print_transformed ( original );


}


#ifdef OAOVERFLOW
int check_bounds_rowsort ( arraydata_t *ad )
{
    long long maxval = 1;

    for ( int x=0; x<ad->ncols; x++ ) {
        if ( maxval> ( std::numeric_limits<rowsort_value_t>::max() /ad->s[x] ) ) {
            logstream ( SYSTEM ) << "check_bounds_rowsort: possible overflow";
        }
        maxval *= ad->s[x];
    }
    return 0;
}
#endif





object_pool< larray<rowindex_t> > arraysymmetry::rowpermpool = object_pool< larray<rowindex_t> >(); // rowpermpool;

arraysymmetry::~arraysymmetry()
{
#ifdef USE_ROWPERM_POOL
    arraysymmetry::rowpermpool.Delete ( rowperm );
    //rowperm=0;
#else
    delete rowperm;
    rowperm=0;
    //rowperm = new larray<rowindex_t>;
#endif

    rowperm=0;
}

arraysymmetry::arraysymmetry ( const arraysymmetry &rhs )
{
#ifdef USE_ROWPERM_POOL
    rowperm = arraysymmetry::rowpermpool.New();
#else
    rowperm = new larray<rowindex_t>;
#endif
    *rowperm = * ( rhs.rowperm );
    colperm = rhs.colperm;
}
arraysymmetry & arraysymmetry::operator= ( const arraysymmetry &rhs )
{
    if ( rowperm==0 ) {
#ifdef USE_ROWPERM_POOL
        rowperm = arraysymmetry::rowpermpool.New();
#else
        rowperm = new larray<rowindex_t>;
#endif
    }
    *rowperm = * ( rhs.rowperm );
    colperm = rhs.colperm;
    return *this;
}

arraysymmetry::arraysymmetry ( const dyndata_t *dyndata )
{
//printf("arraysymmetry: constructor from *dyndata\n");
//    rowperm =   dyndata->getRowperm();

#ifdef USE_ROWPERM_POOL
    rowperm =arraysymmetry::rowpermpool.New();
#else
    rowperm = new larray<rowindex_t>;
#endif

// printf( "xxx %d: ,", (*rowperm).n ); // print_perm( rowperm->d, rowperm->n);
    //rowperm = new larray<rowindex_t>();  printf( "xxx %d\n", (*rowperm).n );

    // *rowperm =   dyndata->getRowperm();  //print_perm( rowperm->d, rowperm->n);
    dyndata->getRowperm ( *rowperm ); //print_perm( rowperm->d, rowperm->n);
    dyndata->getColperm ( colperm );
    //printf("---\n");
    // colpermtypelight x = dyndata->getColperm();
//printf("arraysymmetry: done... %d %d\n", rowperm.n, colperm.n);

    // printf("arraysymmetry::arraysymmetry: done\n");
    // this->show();
}

void LMCreduction_t::symm_t::storeSymmetryPermutation ( const dyndata_t *dyndata )
{
    if ( !store )
        return;
    int n = dyndata->col+1;
    // printf("LMCreduction_t::storeSymmetryPermutation: insert symmetry at %d\n", n);
    arraysymmetry as ( dyndata );
    //   printf("LMCreduction_t::storeSymmetryPermutation: calling insertUnique %d %d %d\n", as.rowperm.n, as.colperm.n, as.lperm.n);
    insertUnique ( symmetries[n], as );

    if ( symmetries[n].size() % 500000 == 0 ) {
        int nrows=symmetries[n][0].rowperm->size();

        printf ( "LMCreduction_t::storeSymmetryPermutation: stored permutation %ld for %d cols\n", ( long ) symmetries[n].size(), n );
        long es=0;
        for ( int i=0; i< ( int ) symmetries.size(); i++ ) {
            es+= nrows*sizeof ( rowindex_t ) *symmetries[i].size()  + symmetries[i].size() *sizeof ( colindex_t ) *n;
        }
        printf ( "LMCreduction_t::storeSymmetryPermutation: raw data size %.1f [MB]\n", double ( es ) / ( 1024*1024 ) );

    }
    //   printf("LMCreduction_t::storeSymmetryPermutation: calling insertUnique ... done\n");
}

LMCreduction_t::LMCreduction_t ( const arraydata_t *adp )
{
    //log_print(SYSTEM, "LMCreduction_t::constructor\n");
    mode = OA_TEST;

    transformation = new array_transformation_t ( adp );
    array = create_array ( adp );
    sd=symmdataPointer ( ( symmdata * ) 0 );
    symms.store=0;
    maxdepth = -1;

    staticdata = 0;


    ncols=adp->ncols;
    nrows=adp->N;
    reset();
#ifdef SWIGCODE
    init_state=COPY;
#else
    init_state=INIT_STATE_INVALID;	// remove this later
#endif
}

LMCreduction_t::~LMCreduction_t()
{
    free();

    releaseStatic();
}

void LMCreduction_t::free()
{
    delete transformation;
    transformation = 0;
    destroy_array ( array );
    array = 0;

}

#define MINCOLMAX 1000

LMCreduction_t::LMCreduction_t ( const LMCreduction_t  &at )
{
    mode = at.mode;
    state=at.state;
    init_state=at.init_state;

    maxdepth=at.maxdepth;
    lastcol=at.lastcol;
    nred=at.nred;
    ncols=at.ncols;
    nrows=at.nrows;

    targetcol=at.targetcol;
    mincol=at.mincol;
    staticdata = 0;

    sd=symmdataPointer ( ( symmdata * ) 0 );

    symms=at.symms;
    //symms.store=at.symms.store;
    //colperms=at.colperms;
    //colcombs=at.colcombs;
    //symmetries=at.symmetries;

    transformation = new array_transformation_t ( * ( at.transformation ) );
    array = create_array ( transformation->ad );

}

LMCreduction_t & LMCreduction_t::operator= ( const LMCreduction_t &at )	/// Assignment operator
{
    mode = at.mode;
    state=at.state;
    init_state=at.init_state;
    maxdepth=at.maxdepth;
    lastcol=at.lastcol;
    nred=at.nred;

    targetcol=at.targetcol;
    mincol=at.mincol;
    ncols=at.ncols;
    nrows=at.nrows;

    staticdata = at.staticdata; // ??

    symms=at.symms;
    //store=at.store;
    //colperms=at.colperms;
    //colcombs=at.colcombs;
    //symmetries=at.symmetries;

    // printf("LMCreduction_t::operator=: colperms.size() %zu\n", colperms.size() );
    free();

    transformation = new array_transformation_t ( * ( at.transformation ) );
    array = create_array ( transformation->ad );

    return *this;
}

void LMCreduction_t::reset()
{
    lastcol=-1;
    nred=0;
    state=REDUCTION_INITIAL;
    init_state=COPY; // ?
    transformation->reset();

    targetcol=0;
    mincol=MINCOLMAX;

    //printf("here: ncols %d\n", ncols);
    clearSymmetries();
    //transformation->reset();

}
void LMCreduction_t::clearSymmetries()
{
    symms.store=0;
    symms.ncols=-1;
    // TODO: move into struct
    symms.colcombs.resize ( ncols+1 );
    symms.colperms.resize ( ncols+1 );
    symms.symmetries.resize ( ncols+1 );
    for ( size_t i=0; i<= ( size_t ) ncols; i++ ) {
        symms.colcombs[i].clear();
        symms.colperms[i].clear();
        symms.symmetries[i].clear();
    }
}

void LMCreduction_t::updateTransformation ( const arraydata_t &ad, const dyndata_t &dyndatacpy, levelperm_t *lperm_p, const array_t *original )
{

    nred++;
    this->state = REDUCTION_CHANGED;
    /* copy level permutations, note: expensive operation */
    for ( colindex_t c=0; c<ad.ncols; c++ )
        copy_perm<array_t> ( lperm_p[c], this->transformation->lperms[c], ad.s[c] );

    /* copy row permutation */
    dyndatacpy.getRowperm ( this->transformation->rperm );
//    for ( rowindex_t x=0; x<ad.N; x++ )
    //      this->transformation->rperm[x] = dyndatacpy.rowsort[x].r;
    /* copy column permutation */
    copy_perm ( dyndatacpy.colperm, this->transformation->cperm, ad.ncols );

    // @pte
    //void debug_show_reduction(reduction, original,  ad, "####3 new reduction:\n");



}

void LMCreduction_t::updateFromLoop ( const arraydata_t &ad, const dyndata_t &dyndatacpy, levelperm_t *lperm_p, const array_t *original )
{

    nred++;
    this->state = REDUCTION_CHANGED;
    /* copy level permutations, note: expensive operation */
    for ( colindex_t c=0; c<ad.ncols; c++ )
        copy_perm<array_t> ( lperm_p[c], this->transformation->lperms[c], ad.s[c] );

    /* copy row permutation */
    dyndatacpy.getRowperm ( this->transformation->rperm );

    /* copy column permutation */
    copy_perm ( dyndatacpy.colperm, this->transformation->cperm, ad.ncols );

    // @pte
    //void debug_show_reduction(reduction, original,  ad, "####3 new reduction:\n");

    /* update array used for comparing */

    this->transformation->apply ( original, this->array );

    if ( 0 ) {
        printf ( "LMCreduction_t::updateFromLoop (red %ld): colperm ", nred );
        print_perm ( this->transformation->cperm, this->transformation->ad->ncols );
        printf ( "new array:\n" );
        print_array ( this->array, ad.N, ad.ncols );
    }
    if ( 0 ) {
        array_link tmp ( original, ad.N, ad.ncols );
        array_link tmp2 ( this->array, ad.N, ad.ncols );

        int rr, cc;
        tmp.firstDiff ( tmp2, rr, cc );
        printf ( "updateFromLoop: first diff col %d, row %d\n", cc, rr );
    }
}


void random_transformation ( array_t *array, const arraydata_t *adp )
{
    array_transformation_t *transformation = new array_transformation_t ( adp );
    transformation->randomize();

    array_t *cpy = clone_array ( array, adp->N, adp->ncols );
    transformation->apply ( cpy, array );
}

/// Apply Hadamard transformation to orthogonal array
void apply_hadamard ( array_link &al, colindex_t hcol )
{
    arraydata_t adata = arraylink2arraydata ( al );
    apply_hadamard ( &adata, al.array, hcol );
}


/**
* @brief Apply Hadamard transformation to orthogonal array
* @param source
* @param target
* @param hcol
*/
void apply_hadamard ( const arraydata_t *ad, array_t *array, colindex_t hcol )
{
    if ( ( ad->N-1 ) !=ad->ncols ) {
        printf ( "WARNING: array does not have size N x (N-1), Hadamard transformation makes no sense\n" );
        //return;
    }
    for ( colindex_t z=0; z<ad->ncols; z++ ) {
        if ( ad->s[z]!=2 ) {
            printf ( "error: Hadamard transformation only defined for 2-level arrays\n" );
            return;
        }
    }

    carray_t idperm[2] = {0,1};
    carray_t swapperm[2] = {1,0};
    carray_t *perm = 0;

    for ( rowindex_t r=0; r<ad->N; r++ ) {
        /* detemine action on this row */
        if ( array[r+ad->N*hcol]==0 )
            perm = idperm;
        else
            perm = swapperm;

        for ( colindex_t c=0; c<ad->ncols; c++ ) {
            if ( c==hcol )
                continue;

            array[c*ad->N+r]=perm[array[ad->N*c+r]];
        }
    }
}


void dyndata_t::reset ( )
{
    init_perm ( this->colperm, this->N );

    for ( int i=0; i<N; i++ ) {
        this->rowsort[i].r = i;
        this->rowsort[i].val = 0;
    }
}

/**
* @brief Constructor for the dyndata_t structure
* @param N_
*/
dyndata_t::dyndata_t ( int N_, int col_ )
{
    this->N = N_;
    this->col = col_;
    this->rowsort = ( rowsort_t * ) malloc ( sizeof ( rowsort_t ) *this->N );
    this->colperm = new_perm_init<colindex_t> ( this->N );	// TODO: initialization with N is nonsense?
    this->rowsortl = 0;

    for ( int i=0; i<N; i++ ) {
        this->rowsort[i].r = i;
        this->rowsort[i].val = 0;	// TODO: not always needed
    }
}

/**
* @brief Make a a deep copy of a dyndata_t structure
* @param dd
*/
dyndata_t::dyndata_t ( const dyndata_t &dd )
{
    //	printf("dyndata_t constructor: from other dyndata_t\n");
    this->rowsortl=0;
    this->rowsort=0;
    //printf("dyndata_t constructor: from other * dyndata_t\n");
    initdata ( dd );

}
/**
* @brief Make a a deep copy of a dyndata_t structure
* @param dd
*/
dyndata_t::dyndata_t ( const dyndata_t *dd )
{
    this->rowsortl=0;
    this->rowsort=0;
    //printf("dyndata_t constructor: from other * dyndata_t\n");
    initdata ( *dd );
}

dyndata_t::~dyndata_t()
{
    free ( this->rowsort );
    delete_perm ( colperm );

    deleterowsortl();
//    if(this->rowsortl!=0) {
    //    delete [] this->rowsortl;
    //  this->rowsortl = 0;
    //}

}

void dyndata_t::initdata ( const dyndata_t &dd )
{
    this->col = dd.col;
    this->N = dd.N;
    this->colperm = new_perm<colindex_t> ( this->N );
    copy_perm ( dd.colperm, this->colperm, N );

    if ( dd.rowsort!=0 ) {
        this->rowsort = ( rowsort_t * ) malloc ( sizeof ( rowsort_t ) *this->N );
        memcpy ( this->rowsort, dd.rowsort, N*sizeof ( rowsort_t ) );
    }
    if ( dd.rowsortl!=0 ) {
        this->rowsortl = new_perm<rowindex_t> ( this->N ); // new rowindex_t [N];
        copy_perm ( dd.rowsortl, this->rowsortl, this->N );

    }

}

void dyndata_t::copydata ( const dyndata_t &dd )
{
    if ( this->col==dd.col ) {
        copy_perm ( dd.colperm, this->colperm, N );
    } else {
        this->col = dd.col;
        delete_perm ( colperm );
        this->colperm = new_perm<colindex_t> ( this->N );
        copy_perm ( dd.colperm, this->colperm, N );
    }


    if ( this->N == dd.N ) {
        if ( dd.rowsort!=0 ) {
            memcpy ( this->rowsort, dd.rowsort, N*sizeof ( rowsort_t ) );
        }
        if ( dd.rowsortl!=0 ) {
            copy_perm ( dd.rowsortl, this->rowsortl, this->N );

        }
    } else {
        this->N = dd.N;
        if ( dd.rowsort!=0 ) {
            if ( this->rowsort!=0 ) {
                free ( this->rowsort );
                this->rowsort=0;
            }

            this->rowsort = ( rowsort_t * ) malloc ( sizeof ( rowsort_t ) *this->N );
            memcpy ( this->rowsort, dd.rowsort, N*sizeof ( rowsort_t ) );
        }
        if ( dd.rowsortl!=0 ) {
            deleterowsortl();
            this->rowsortl = new_perm<rowindex_t> ( this->N ); // new rowindex_t [N];
            copy_perm ( dd.rowsortl, this->rowsortl, this->N );

        }
    }
}

/*
void dyndata_t::update( const dyndata_t &dd )
{
    // clear old data
  if(this->rowsort!=0) {
    free ( this->rowsort );
    this->rowsort=0;
  }
    delete_perm ( colperm );

    deleterowsortl();

    //printf("dyndata_t::operator=: from other dyndata_t\n");
    initdata(dd);

 } */

dyndata_t & dyndata_t::operator= ( const dyndata_t &dd )
{
    this->copydata ( dd );
    return *this;
}



void dyndata_t::show() const
{
    printf ( "dyndata_t: N %d, col %d; ", this->N, this->col );
    printf ( "  colperm: " );
    print_perm ( this->colperm, this->col );
    printf ( "  rowsort:\n" );
    print_rowsort ( this->rowsort, this->N );
}

/**
* @brief Convert column permutation of the root to a sorting of the rows (to make the array in root form again)
* @param tmprperm Current rowperm
* @param rperm New permutation
* @param cp Column permutation
*/
inline void cperm_lperms_to_rowsort ( const rowperm_t tmprperm, levelperm_t *lperms, rowperm_t rperm, const colperm_t cp, const arraydata_t *ad )
{
    // assert(ad->s[0]==2); // TODO: check this

    array_t *mults = new array_t[ad->strength+1];

    /// invert column permutation
    colperm_t cpi = invert_permutation ( cp, ad->strength );


    int mult=1; //mults[strength]=mult;
    for ( int l=ad->strength-1; l>=0; l-- ) {
        mults[l]=mult;
        mult*=ad->s[cpi[l]];
    }
    //printf("mults: "); print_perm(mults,ad->strength);


    // loop over all rows
    for ( rowindex_t k=0; k<ad->N; k++ ) {
        //int mult = 1;
        rowindex_t v = k;

        /* calculate the value of the permutation, modulo blocks of oaindex */
        rperm[k] = 0;
        v /= ad->oaindex;
        //printf(" row %d: adding ", k);
        for ( int l=ad->strength-1; l>=0; l-- ) {
            int nl = cp[l];		// new column
            int vmod = v% ( ad->s[nl] );	// value in this row
            v /= ad->s[l];
            rperm[k] += mults[nl]*lperms[l][vmod];
            //	printf(" %d",mults[nl]*lperms[l][vmod]);

        }
        //printf(", total %d ", rperm[k]);
        /* translate value into new position */
        rperm[k] *= ad->oaindex;
        rperm[k] += ( k%ad->oaindex );
        //printf(", fix %d \n", rperm[k]);
    }

    delete_perm ( cpi );
    delete [] mults;
}

/**
* @brief Convert a set of level permutations to the corresponding row permutations (for an array in root form)
* @param rperm
* @param lperms
* @param ad
*/
inline void lperms_to_rowsort ( rowperm_t rperm, levelperm_t *lperms, const arraydata_t *ad )
{
    // loop over all rows
    for ( rowindex_t k=0; k<ad->N; k++ ) {
        int mult = 1;
        rowindex_t v = k;

        /* calculate the value of the permutation, modulo blocks of oaindex */
        rperm[k] = 0;
        v /= ad->oaindex;
        for ( int l=ad->strength-1; l>=0; l-- ) {
            int vmod = v% ( ad->s[l] );
            rperm[k] += mult*lperms[l][vmod];
            v /= ad->s[l];
            mult*=ad->s[l];
        }
        /* translate value into new position */
        rperm[k] *= ad->oaindex;
        rperm[k] += ( k%ad->oaindex );
    }
}

/** @brief Helper function for create_root_permutations_index()
*/
void create_root_permutations_index_helper ( rowperm_t *rperms, levelperm_t *lperms, const arraydata_t *ad, int level, int &permcounter )
{
    // OPTIMIZE: this function can be made much faster?
    if ( level==ad->strength ) {
        /* combine level permutations into root row permutations */

        lperms_to_rowsort ( rperms[permcounter], lperms, ad );

        permcounter++;
    } else {
        int nlperms = init_perm_n<array_t, int> ( lperms[level], ad->s[level] );

        // loop over permutations
        for ( int i=0; i<nlperms; i++ ) {
            create_root_permutations_index_helper ( rperms, lperms, ad, level+1, permcounter );
            next_perm ( lperms[level], ad->s[level] );
        }
    }

    return;
}


/**
* @brief Create in index of all root level permutations
* @param ad Parameters of the array
* @param totalpermsr Returns the total number of permutations found
* @return
*/
rowperm_t* create_root_permutations_index ( const arraydata_t *ad, int &totalpermsr )
{
    /* OPTIMIZE: can this function be made simpler by noting that all permutations of N/oaindex occur ??*/

    totalpermsr = 1;
    for ( int i=0; i<ad->strength; i++ )
        totalpermsr *= factorial<int> ( ad->s[i] );

    log_print ( DEBUG, "create_root_permutations_index: allocating rootrowperms table of size %d*%d*sizeof(rowsort_t)=%ld=%.1f MB (sizeof(rowsort_t)=%d)\n", totalpermsr, ad->N, ( long ) totalpermsr* ( long ) ad->N* ( long ) sizeof ( rowsort_t ), ( long ) totalpermsr* ( long ) ad->N* ( ( double ) sizeof ( rowsort_t ) / ( 1024*1024.0 ) ), sizeof ( rowsort_t ) );

    rowperm_t *rperms = malloc2d<rowindex_t> ( totalpermsr, ad->N );
    array_t **lperms = malloc2d_irr<array_t> ( ad->strength, ( ad->s ) );

    int totalperms2 = 0;
    create_root_permutations_index_helper ( rperms, lperms, ad, 0, totalperms2 );

    free2d ( lperms );
    return rperms;
}

/**
* @brief Create in index of all root level permutations (both level and column permutations)
* @param ad Parameters of the array
* @param totalpermsr Returns the total number of permutations found
* @return
*/
rowperm_t* create_root_permutations_index_full ( const arraydata_t *ad, int &totalpermsr )
{
    if ( ! ( ad->colgroupsize[0]>=ad->strength ) ) {
        // printf("mixed array, not allocating\n");
        if ( log_print ( NORMAL, "" ) ) {
            printf ( "mixed array, not allocating\n" );
            printf ( "ad->strength %d, ad->colgroupsize[0] %d\n", ad->strength, ad->colgroupsize[0] );
            ad->show();
            ad->show_colgroups();
        }
        return 0;
    }
    //   assert ( ad->colgroupsize[0]>=ad->strength );	// make sure root is uniform
    //assert(ad->s[0]==2);	// make sure we are dealing with 2-factor arrays

    totalpermsr = static_cast<int> ( pow ( 2.0, ad->strength ) * factorial<int> ( ad->strength ) );

    log_print ( DEBUG, "create_root_permutations_index_full: allocating rootrowperms table of size %d*%d*sizeof(rowsort_t)=%ld=%.1f MB (sizeof(rowsort_t)=%d)\n", totalpermsr, ad->N, ( long ) totalpermsr* ( long ) ad->N* ( long ) sizeof ( rowsort_t ), ( long ) totalpermsr* ( long ) ad->N* ( ( double ) sizeof ( rowsort_t ) / ( 1024*1024.0 ) ), sizeof ( rowsort_t ) );

    rowperm_t *rperms = malloc2d<rowindex_t> ( totalpermsr, ad->N );
    array_t **lperms = malloc2d_irr<array_t> ( ad->strength, ( ad->s ) );

    rowperm_t tmprperm = new_perm<rowindex_t> ( ad->N );

    colperm_t cp = new_perm_init<colindex_t> ( ad->strength );
    int permcounter=0;
    for ( int i=0; i<factorial<int> ( ad->strength ); i++ ) {
        // loop over all columns permutations
        for ( int j=0; j<pow ( double ( 2 ), ad->strength ); j++ ) {
            // loop over all possible level permutations
            root_row_permutation_from_index ( j,  ad, lperms ); // TODO: this can be done faster

            cperm_lperms_to_rowsort ( tmprperm, lperms, rperms[permcounter],cp, ad );
            permcounter++;

        }
        next_perm ( cp, ad->strength );
    }
    delete_perm ( cp );
    delete_perm ( tmprperm );

    if ( log_print ( DEBUG+1,"" ) ) {
        printf ( "create_root_permutations_index_j4: strength %d\n", ad->strength );
        for ( int k=0; k<permcounter; k++ ) {
            printf ( "perm %3d: ", k );
            print_perm ( rperms[k], ad->N );
        }
    }

    free2d ( lperms );	// we do not need to store the lperms
    return rperms;
}



/**
* @brief From a given permutation index construct the corresponding level permutations
* @param permindex
* @param ad
* @param lperms
*/
void root_row_permutation_from_index ( int permindex,  const arraydata_t *ad, levelperm_t *lperms )
{
    //printf("permindex: %d\n", permindex);
    for ( colindex_t i=0; i<ad->strength; i++ ) {
        const int col = ad->strength-i-1;
        int fac = factorial<int> ( ad->s[col] );

        int idx = ( permindex % fac );
        int rem = ( permindex / fac );

        init_perm ( lperms[col], ad->s[col] );

        permutationLex<array_t> ( idx, lperms[col], ad->s[col] );
        permindex=rem;
    }
}


template <class numtype>
/* create permutation that induces a certain combination selection at a starting point */
void create_perm_from_comb2 ( numtype *perm, const numtype *comb, int k, int ncols )
{
    // initialize
    // init_perm ( perm, ncols );

    int *tmpperm = new_perm_init<int> ( ncols );
    for ( int i=0; i<k; i++ ) {
        perm[i]=tmpperm[comb[i]];
        tmpperm[comb[i]]=-1;
    }
    int j=0;
    for ( int i=k; i<ncols; i++ ) {
        //printf("i %d: j %d, tmpperm[j] %d\n", i, j, tmpperm[j]);
        while ( tmpperm[j]<0 )
            j++;
        perm[i]=tmpperm[j];
        j++;
    }
    delete_perm ( tmpperm );
}

template <class numtype>
/* create permutation that induces a certain combination selection at a starting point. Note we assume the combination is ordered! */
void create_perm_from_comb ( numtype *perm, const numtype *comb, int k, int ncols, numtype startcol )
{
    // initialize
    init_perm ( perm, ncols );

    // this assumes comb is ordered!!??
#ifdef OADEBUG
    for ( int i=0; i< ( k-1 ); i++ )
        if ( comb[i]>=comb[i+1] )
            printf ( "create_perm_from_comb: ERROR!\n" );
#endif
    for ( int i=0; i<k; i++ )
        std::swap ( perm[i+startcol], perm[comb[i]+startcol] );
}



#ifdef OAMEM
inline void static_init_rp ( rowperm_t &rp, rowindex_t N )
{
    static rowindex_t Np;
    static rowperm_t rowp = 0;


    if ( N!=Np ) {
        Np = N;
        if ( rowp!=0 )
            delete_perm ( rowp );
        rowp = new_perm<rowindex_t> ( N );
    }

    rp = rowp;
}
#endif

void LMC_static_struct_t::init_root_stage ( levelperm_t * &lperm_p, colperm_t * &colperm_p, const arraydata_t *adp )
{
    //static_update ( adp );

    /* permutations buffer */
    lperm_p = this->current_trans->lperms;
    colperm_p = this->colperm_p;

    this->LMC_root_init = 1;
}

/**
* @brief Do the static initialization of structures necessary for LMC test and reduction
*
* The argument is used to check whether the static structures need to be updated or not
*
* @param adp
*/
void LMC_static_struct_t::freeall ( )
{
    //arraydata_t *adp = this->ad;
    // log_print ( DEBUG, "LMC_static_struct_t::free: adp->ncols %d, adp->strength %d\n", adp->ncols, adp->strength );


    /* clear old structures */
    if ( this->current_trans!=0 ) {
        //cout << "LMC_static_struct_t::current_trans " << LMC_static_struct_t::current_trans << endl;
        delete this->current_trans;

#ifdef OADEBUG
        this->current_trans=0;
#endif
    }

    if ( this->colperm_p !=0 ) {
        free2d<colindex_t> ( this->colperm_p );
    }
    if ( this->localcolperm_p !=0 ) {
        free2d<colindex_t> ( this->localcolperm_p );
    }

    if ( this->dyndata_p!=0 ) {
        if ( this->ad!=0 ) {

            for ( int z=0; z<this->ad->ncols; z++ )
                delete this->dyndata_p[z];
            delete [] this->dyndata_p;
        }
    }

    if ( this->rootrowperms!=0 )
        free2d ( this->rootrowperms );
    if ( this->rootrowperms_full!=0 )
        free2d ( this->rootrowperms_full );

    if ( this->ad !=0 )
        delete this->ad;

    if ( this->colbuffer!=0 )
        free ( this->colbuffer );

    /* make sure all variables are initialized again */
    this->LMC_non_root_init=0;
    this->LMC_root_init=0;
    this->LMC_root_rowperms_init=0;

}

void LMC_static_struct_t::init ( const arraydata_t *adp )
{
    log_print ( DEBUG, "static_init_LMC: adp->ncols %d, adp->strength %d\n", adp->ncols, adp->strength );

    /* clear old structures */
    this->freeall();


    /* allocate new structures */
    this->ad = new arraydata_t ( *adp );
    this->current_trans = new array_transformation_t ( adp );
    this->colperm_p = malloc2d<colindex_t> ( this->ad->ncols, this->ad->ncols );
    this->localcolperm_p = malloc2d<colindex_t> ( this->ad->ncols, this->ad->ncols );

    /* level permutations of the root structure */
#ifdef OAMEM
    /* for memory efficient version of the algorithm do not create a root_permutations index */
#else
    this->rootrowperms = create_root_permutations_index ( this->ad,this->nrootrowperms );
    this->rootrowperms_full = create_root_permutations_index_full ( this->ad,this->nrootrowperms_full );

#endif

    this->dyndata_p = new dyndata_t * [adp->ncols];
    for ( int z=0; z<adp->ncols; z++ )
        this->dyndata_p[z] = new dyndata_t ( adp->N );

//    if ( this->colbuffer!=0 )
//        free ( this->colbuffer );
    this->colbuffer = ( array_t * ) malloc ( sizeof ( array_t ) *adp->N );

    /* make sure all variables are initialized again */
    this->LMC_non_root_init=0;
    this->LMC_root_init=0;
    this->LMC_root_rowperms_init=0;

}

int LMC_static_struct_t::needUpdate ( const arraydata_t *adp ) const
{
    int update = 0;
    //printf("needUpdate: this->ad %ld\n", long(this->ad) );
    if ( this->ad==0 )
        update = 1;
    else {
        if ( ! ( * ( this->ad ) ==*adp ) ) {
            // 			printf(" LMC_static_struct_t::update: update=1: ad\n");
            // 		this->ad->show();
            // 			printf("  adp\n");
            // 			adp->show();

            update = 1;
        }
    }
    return update;

}

int LMC_static_struct_t::update ( const arraydata_t *adp )
{
    int update = this->needUpdate ( adp );

    //update=1; /* be carefull with updates, other functions might already have initialized the static variables */
    if ( update ) {
        this->init ( adp );
    }
    return update;
}


void LMC_static_struct_t::init_nonroot_stage ( levelperm_t * &lperm_p, colperm_t * &colperm_p, colperm_t * &localcolperm_p, dyndata_t ** &dynd_p, int &dynd_p_nelem,  array_t * &colbuffer, const arraydata_t *adp ) const
{
    int updatevars=1;	// always update variables since there is a second static init function

    // we assume the static initialization has already been done

    //updatevars=1;
    if ( updatevars!=0 ) {
        /* permutations buffer */
        lperm_p = this->current_trans->lperms;
        colperm_p = this->colperm_p;
        localcolperm_p = this->localcolperm_p;

        /* dynamic_data buffer */
        dynd_p = this->dyndata_p;
        dynd_p_nelem = this->ad->ncols;

        colbuffer = this->colbuffer;
    }
}



/**
* @brief Sort an array based on the values of the root
* @param array
* @param ad
* @param dyndata Contains information about current transformation of the array
* @param rowsort
*/
inline void LMC_root_sort ( const carray_t *array, const arraydata_t *ad, const dyndata_t *dyndata )
{
    vindex_t *valueindex = new_valueindex<vindex_t> ( ad->s, ad->strength );	// OPTIMIZE: make static allocation?

    rowsort_t *rowsort = dyndata->rowsort;
    /* perform initial sort, after the initial sort we can sort using pre-calculated structures */
    for ( rowindex_t j = 0; j < ad->N; j++ ) {
        rowsort[j].val = row_rank_partial ( array,	 ad->N, valueindex, ( colindex_t ) 0, ad->strength, dyndata->colperm, j );
    }
    delete [] valueindex;

    // OPTIMIZE: select best sort (stable_sort, sort, oacolSort)
    std::stable_sort ( rowsort, rowsort+ad->N );
}



/**
* @brief Sort an array based on the values of the root
* @param array
* @param ad
* @param rowsort
*/
inline void LMC_root_sort ( carray_t *array, const arraydata_t *ad, rowsort_t *rowsort )
{
    int *valueindex = new_valueindex<int> ( ad->s, ad->strength );	// OPTIMIZE: make static allocation?

    /* perform initial sort, after the initial sort we can sort using pre-calculated structures */
    for ( rowindex_t j = 0; j < ad->N; j++ ) {
        rowsort[j].val = row_rank_partial ( array+rowsort[j].r, 0, ad->strength, ad->N, valueindex );
    }
    delete [] valueindex;

    // OPTIMIZE: select best sort
    std::sort ( rowsort, rowsort+ad->N );
    //bubbleSort(rowsort, ad->N);
}

/** @brief Static initialization of level permutations
 */
inline void static_init_lperms ( const arraydata_t *adp, levelperm_t * &lperm_p, LMC_static_struct_t &tmpStatic )
{
    /* no static update, we assume this has been done already */
    tmpStatic.update ( adp );

    lperm_p =  tmpStatic.current_trans->lperms;
}

/** @brief Static initialization of root row permutations
*/
inline void static_init_rootrowperms ( const arraydata_t *adp, int &totalperms, rowperm_t * &rootrowperms, levelperm_t * &lperm_p,  LMC_static_struct_t &tmpStatic )
{
    /* no static update, we assume this has been done already */
    //static_update ( adp );

    totalperms = tmpStatic.nrootrowperms;
    rootrowperms = tmpStatic.rootrowperms;
    lperm_p =  tmpStatic.current_trans->lperms;

    tmpStatic.LMC_root_rowperms_init = 1;
}

/*
inline void static_init_reduce_rootrowperms ( const arraydata_t *adp, int &totalperms, rowperm_t * &rootrowperms, LMC_static_struct_t &tmpStatic) {
    tmpStatic.update ( adp );

    totalperms = tmpStatic.nrootrowperms;
    rootrowperms = tmpStatic.rootrowperms;

    tmpStatic.LMC_reduce_root_rowperms_init = 1;
}
*/

/// show array with rowsort and colperm
void show_array ( carray_t *array, const int ncols, const int nrows, colperm_t colperm, rowsort_t *rs )
{
    int count;
    for ( int j = 0; j < nrows; j++ ) {
        count = j;
        for ( int k = 0; k < ncols; k++ ) {
            const char *s = ( k<ncols-1 ) ? " ":"\n";
            int val = array[nrows*k + rs[j].r];
            printf ( "%3i%s", static_cast<int> ( val ), s );
            count += nrows;
        }
    }
}

/**
* @brief Perform root level permutations in LMC check (with full root permutation group)
*
* This method assumes J4
*
* @param original
* @param array
* @param ad Static array information
* @param dyndata Dynamic array information
* @return
* @see LMC
*/
lmc_t LMCreduce_root_level_perm_full ( carray_t const *original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic )
{
    //     	if ( array_cmp(original, array, ad->N, ad->ncols) ) {
    // 		printf("LMCreduce_root_level_perm: Huh?\n");
    // 				print_array("original:\n", original, ad->N, ad->ncols);
    // 		print_array("array:\n", array, ad->N, ad->ncols);
    // 		printf("LMCreduce_root_level_perm: Huh?\n");
    //
    // 		exit(0);
    // 	}

    lmc_t ret = LMC_EQUAL;
    rowsort_t *rowsort = dyndata->rowsort;


    /* perform initial sorting of the root */
    LMC_root_sort ( original, ad, dyndata );

    int totalperms = 0;
    rowperm_t* rootrowperms = 0;
    levelperm_t * lperm_p = 0;
    tmpStatic.init_rootrowperms_full ( totalperms, rootrowperms, lperm_p );
    assert ( rootrowperms );

    /* perform fake level permutations */
    dyndata_t dyndatatmp = dyndata_t ( dyndata );


    //printf("   LMCreduce_root_level_perm_full %d transforms\n", totalperms);
    for ( int l=0; l<totalperms; l++ ) {
        // update sort structure
        // OPTIMIZE: only perform partial copy
        // OPTIMIZE: not needed if there is only 1 root stage?
        for ( rowindex_t k=0; k<ad->N; k++ ) {
            //printf(" %d %d %d\n", l, k, rootrowperms[l][k]);
            dyndatatmp.rowsort[rootrowperms[l][k]].r=rowsort[k].r;
        }

        //printf("LMCreduce_root_level_perm: l %d/%d\n", l, totalperms);
        //dyndatatmp.show();

        /* pass to non-root stage */
        ret = LMCreduce_non_root_j4 ( original, ad, &dyndatatmp, reduction, oaextend, tmpStatic );
        //ret = LMCreduce_non_root( original, ad, &dyndatatmp, reduction, tmpStatic );

        if ( ret==LMC_LESS ) {
            break;
        }
    }

    return ret;
}

/**
* @brief Perform root level permutations in LMC check
* @param original
* @param array
* @param ad Static array information
* @param dyndata Dynamic array information
* @return
* @see LMC
*/
lmc_t LMCreduce_root_level_perm ( array_t const *original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic )
{


#ifdef OAANALYZE_DISCR
    analysis_increase_counter ( "root_level_perm" );
#endif

    lmc_t ret = LMC_EQUAL;
    rowsort_t *rowsort = dyndata->rowsort;

    /* perform initial sorting of the root */
    LMC_root_sort ( original, ad, dyndata );
#ifdef OADEBUG2
    printf ( "##### LMCreduce_root_level_perm\n" );
    print_array ( original, ad->N, ad->ncols );
    dyndata->show();
    reduction->show ( 1 );
#endif

    int totalperms = 0;
    rowperm_t* rootrowperms = 0; // pointer to root row permutations
    levelperm_t * lperm_p = 0;
    static_init_rootrowperms ( ad, totalperms, rootrowperms, lperm_p, tmpStatic );

    /* perform fake level permutations */
    dyndata_t dyndatatmp = dyndata_t ( dyndata );

#ifdef OADEBUG
    //oaextend.info();
#endif
    //      printf("#############################\nLMCreduce_root_level_perm: totalperms %d colperm ", totalperms); print_perm(dyndata->colperm, ad->ncols);


    for ( int l=0; l<totalperms; l++ ) {	// loop over root permutations (in levels)
        // update sort structure
        // OPTIMIZE: only perform partial copy

        // TODO: update dyndatatmp rowsort with generic function (e.g. rowsortl if needed)
        for ( rowindex_t k=0; k<ad->N; k++ ) {
            // TODO: is this valid if the root is not properly sorted?
            // 	this assumes that after the LMC_root_sort the root is already in blocks
            dyndatatmp.rowsort[rootrowperms[l][k]].r=rowsort[k].r;

#ifdef OADEBUG
            if ( dyndatatmp.rowsort[rootrowperms[l][k]].r >= ad->N ) {
                printfd ( "error with rowsort values" );
                exit ( 0 );
            }
#endif
        }

        //printfd("l %d/%d\n", l, totalperms); reduction->show();

        // update level permutations structure
        if ( reduction->mode>=OA_REDUCE ) {
            root_row_permutation_from_index ( l, ad, lperm_p );
        }

        /* pass to non-root stage */
        //oaextend.info();
        //printf("oaextend.getAlgorithm() %d, reduction->sd %ld\n", oaextend.getAlgorithm(), (long)reduction->sd.get() );
        if ( oaextend.getAlgorithm() ==MODE_LMC_2LEVEL || ( oaextend.getAlgorithm() ==MODE_J5ORDERX && reduction->sd!=0 ) || ( oaextend.getAlgorithm() ==MODE_J5ORDERXFAST && reduction->sd!=0 ) ) {
//	   myassert(reduction->sd!=0, "reduction_t.sd should be set!" );
            dyndatatmp.initrowsortl();
            ret = LMCreduce_non_root_2level ( original, ad, &dyndatatmp, reduction, oaextend, tmpStatic );	// TODO: this should be a call to j4/tplus code
        } else
            ret = LMCreduce_non_root ( original, ad, &dyndatatmp, reduction, oaextend, tmpStatic );	// TODO: this should be a call to j4/tplus code

        if ( ret==LMC_LESS ) {
            if ( 0 ) {
                printf ( "LMC_non_root: l %d returned LMC_LESS\ntransformed array:\n", l );
            }
            break;
        }
    }

    return ret;
}

#ifdef OAMEM
lmc_t LMCreduce_root_level_perm_ME ( carray_t const *original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction )
{
    lmc_t ret = LMC_EQUAL;
    rowsort_t *rowsort = dyndata->rowsort;

    levelperm_t * lperm_p = 0;
    static_init_lperms ( ad, lperm_p );

    logstream ( DEBUG+1 ) << "LMCreduce_root_level_perm_ME: col " << dyndata->col << endl;
    if ( dyndata->col==0 ) {
        //print_array("LMCreduce_root_level_perm_ME: array\n", array, ad->N, ad->ncols);

        /* perform initial sorting of the root */
        LMC_root_sort ( original, ad, dyndata );

        /* initialize level permutations */
        for ( int x=0; x<ad->strength; x++ ) {
            init_perm<array_t> ( lperm_p[x], ad->s[x] );
        }
    }


    if ( dyndata->col==ad->strength ) {
        /* pass to non-root stage */

        dyndata_t dyndatatmp = dyndata_t ( dyndata ); // OPTIMIZE: make static allocation
        rowperm_t rp;
        static_init_rp ( rp, ad->N );

        lperms_to_rowsort ( rp, lperm_p, ad );
        for ( rowindex_t k=0; k<ad->N; k++ )
            dyndatatmp.rowsort[rp[k]].r=dyndata->rowsort[k].r;

        /* pass to non-root stage */
        //reduction->show();
        ret = LMCreduce_non_root ( original, ad, &dyndatatmp, reduction );

        return ret;
    }

    colindex_t col = dyndata->col;

    /* perform level permutations */
    bool np = true;
    unsigned long nperms = 0;
    while ( np ) {
        dyndata_t *dyndatatmp = new dyndata_t ( dyndata );
        dyndatatmp->col++;
        ret = LMCreduce_root_level_perm_ME ( original, ad, dyndatatmp, reduction );
        delete dyndatatmp;

        if ( ret==LMC_LESS ) {
            if ( 0 && checkloglevel ( NORMAL ) ) {
                printf ( "LMC_non_root: np %d returned LMC_LESS\ntransformed array:\n", np );
                show_array ( original, ad->ncols, ad->N, dyndata->colperm, dyndata->rowsort );
                print_array ( "original:\n" , original, ad->N, ad->ncols );
                printf ( "------\n" );
                print_array ( "LMC_non_root: reduction->array:\n", reduction->array, ad->N, ad->ncols );
            }
            break;
        }
        np = next_permutation ( lperm_p[col], lperm_p[col]+ad->s[col] );

        if ( col==0 && ( nperms%10000 == 0 || nperms<20 ) ) {
            logstream ( DEBUG ) << "LMCreduce_root_level_perm_ME: doing permutation " << nperms << " (of first col)" << endl;
        }
        nperms++;

    }
    return ret;
}

#endif

vector<int> LMCcheckLex ( arraylist_t const &list, arraydata_t const &ad, int verbose )
{
    vector<int> check ( list.size() );

    for ( unsigned int x=0; x<list.size(); x++ ) {
        if ( verbose )
            printf ( "LMCcheck: array %d\n", x );
        check[x] = LMCcheckLex ( list[x], ad );
    }
    return check;
}


static int nfrac1=0;
static int nfrac2=0;
void print_fracs ( int logl )
{
    log_print ( logl, "#### LMCcheckj4: total %d discard %d: gain %.2f\n ", nfrac1, nfrac2, double ( nfrac2 ) /nfrac1 );

}
void clear_fracs()
{
    nfrac1=0;
    nfrac2=0;
}




/** @brief Calculate J-characteristic for a column combination
*
* We assume the array has values 0 and 1
*/
int fastj3 ( carray_t *array, rowindex_t N, const int J, const colindex_t *pp )
{
    //carray_t *x = ar.array;
    int jval=0;

    for ( rowindex_t r=0; r<N; r++ ) {
        array_t tmp=0;
        for ( int i=0; i<J; i++ ) {
            tmp+=array[r+ pp[i]*N];
        }
        tmp %= 2;
        jval += tmp;
    }
    jval = 2*jval-N;
    return ( jval );
}

lmc_t LMCcheckj4 ( array_link const &al, arraydata_t const &adin, LMCreduction_t &reduction, const OAextend &oaextend, int jj )
{
    LMC_static_struct_t &tmpStatic = getGlobalStaticOne();
    //tmpStatic.setRef("LMCcheckj4");

    const int maxjj = 40;
    assert ( jj<maxjj );

    if ( reduction.init_state==INIT_STATE_INVALID ) {
        // TODO: remove this code
        printf ( "LMCreduce: reduction.init_state is INIT_STATE_INVALID, please set it to INIT_STATE::COPY or INIT!\n" );
        fflush ( stdout );
        std::cout.flush();
//	sleep(2);
//	std::cout.flush();
        throw;
    }
    if ( reduction.init_state==COPY ) {
        reduction.setArray ( al );


        // printf("LMCreduce: copy array: %d %d\n", reduction->nrows, reduction->ncols);
        //reduction->getArray().showarray();
    }

    // printf("LMCcheckj4: ad.ncols %d\n", ad.ncols);
    assert ( adin.ncolgroups==1 );
    assert ( adin.s[0]==2 );
    assert ( jj>=adin.strength );
    if ( adin.strength!=3 ) {
        printf ( "LMCcheckj4: error: strength not equal to 3\n" );
        throw;
    }

    // 	printf("entering LMCcheckj4:\n");
    // 	print_array("al.array\n", al.array, adin.N, adin.ncols);
    // 	print_array("reduction.array\n", reduction.array, adin.N, adin.ncols);
    //
    if ( al.n_columns<jj ) {
        // fall back to original method
        dyndata_t dyndata ( adin.N );

        lmc_t ret = LMCreduce ( al.array, al.array, &adin, &dyndata, &reduction, oaextend );
        return ret;
    }

    arraydata_t ad ( adin );
    ad.complete_arraydata_splitn ( jj );	// split column group!

    // TODO: this bites with ad and adin, they are different?!
    tmpStatic.update ( &ad );

    int nc = ncombs ( ad.ncols, jj );
    colindex_t *comb = new_comb_init<colindex_t> ( jj );
    int nrootcombs = ncombs ( jj, ad.strength );
    colindex_t *combroot = new_comb_init<colindex_t> ( ad.strength );	//
    colindex_t combx[maxjj];	// final combination for root selection + extra columns
    colindex_t comby[maxjj];

    carray_t *array = al.array;

    lmc_t ret = LMC_EQUAL;

    /* loop over all possible column combinations for the first jj columns */
    int goodj = abs ( jvaluefast ( array, ad.N, jj, comb ) );

    /* initialize variables outside of loop */
    colperm_t perm = new_perm<colindex_t> ( ad.ncols );
    colperm_t pp = new_perm_init<colindex_t> ( jj );


    for ( int i=0; i<nc; i++ ) {
        nfrac1++;

        int jval=abs ( jvaluefast ( array, ad.N, jj, comb ) );
        //printf("LMCcheckj4: colcomb %d, goodj %d, jval %d\n", i, goodj, jval);
        //int jval=0;
        //    printf("X jval %d, jgood %d\n", jval, goodj);

        if ( jval>goodj ) {
            nfrac2++;
            // this combination will lead to LMC less
            if ( reduction.mode>=OA_REDUCE ) {
                ret=LMC_EQUAL;
                goodj=jval;	// update state

            } else {
                ret=LMC_LESS;
                break;
            }
            //printf("jval %d, jgood %d\n", jval, goodj);
        } else if ( jval<goodj ) {
            nfrac2++;
            // this combination can only lead to LMC more
#ifdef OADEBUG2
            printf ( "jval %d, jgood %d\n", jval, goodj );
#endif
            //printf("checkj4: MORE jval %d, jgood %d\n", jval, goodj);
            next_combination ( comb, jj, ad.ncols );	// increase combination
            continue;
        }
        // do classic check

#ifdef OADEBUG2
        printf ( "#### LMCcheckj4: jval %d, classic checks:\n ", jval );
#endif

        init_comb ( combroot, ad.strength, jj );
        for ( int kk=0; kk<nrootcombs; kk++ ) {
            // prepare for all root combinations

            // convert combination to full permutation
            create_perm_from_comb ( combx, combroot, ad.strength, jj, 0 );
            perform_inv_perm ( comb, comby, jj, combx );

            init_perm ( perm, ad.ncols );
            create_perm_from_comb2<colindex_t> ( perm, comby, jj, ad.ncols );


            if ( 0 ) {
                printf ( "LMCcheckj4: kk %d, jval %d, classic check on perm: ", kk, jval );
                print_perm ( perm, ad.ncols );
                printf ( "combroot: " );
                print_perm ( combroot, ad.strength );
                printf ( "comb: " );
                print_perm ( comb, jj );
                printf ( "combx: " );
                print_perm ( combx, jj );
            }

            dyndata_t dyndata ( ad.N );
            dyndata.col=ad.strength;	// set at col after root
            copy_perm ( perm, dyndata.colperm, ad.ncols );


            ret = LMCreduce_root_level_perm_full ( array, &ad, &dyndata,&reduction, oaextend, tmpStatic );

            if ( reduction.mode>=OA_REDUCE ) {
                if ( ret==LMC_LESS )
                    reduction.state = REDUCTION_CHANGED;
                ret=LMC_EQUAL;

                // do nothing?
                //printf("LMCcheckj4: reduction changed\n");
                //reduction.show();
                //reduction.transformation->show();
            } else {
                if ( ret==LMC_LESS ) {
                    if ( 0 ) {
                        printf ( "LMCcheckj4: reduction changed: " );
                        print_perm ( dyndata.colperm, 4 );
                        reduction.show();
                        reduction.transformation->show();

                        printf ( "xxx %d\n", reduction.getArray() ==al );


                        printf ( "original array:\n" );
                        print_array ( array, ad.N, ad.ncols );
                        printf ( "   full transform:\n" );
                        //show_array ( array, ad.ncols, ad.N, dyndata.colperm, dyndata.rowsort );
                    }
                    break;
                }
            }
            next_combination ( combroot, ad.strength, jj );	// increase combination

        }
        if ( ret==LMC_LESS ) {
            break;
        }

        next_combination ( comb, jj, ad.ncols );	// increase combination
    }
    delete_perm ( perm );
    delete_perm ( pp );
    delete_comb ( combroot );
    delete_comb ( comb );

#ifdef OADEBUG2
    printf ( "returning %d\n", ret );
#endif


    if ( nfrac1%40000==0 ) {
        print_fracs();
    }
    return ret;
}




/// create splits for a symmetry group
std::vector<int> symmetrygroup2splits ( const symmetry_group &sg, int ncols, int verbose=0 )
{

    int jj= sg.n;

    std::vector<int> splits;
    //splits.push_back(0);

    //printf("symmetry indices: "); display_vector(sg.gidx); printf("\n");
    if ( verbose>=2 ) {
        printf ( "symmetry indices: %d: ", sg.ngroups );
        display_vector ( sg.gstart );
        printf ( "\n" );
    }
    for ( int j=0; j<sg.ngroups; j++ ) {
        splits.push_back ( sg.gstart[j] );
    }
    if ( ncols>jj ) {
        splits.push_back ( jj );
    }
    if ( sg.ngroups<2 && 0 ) {
        printf ( "symmetrygroup2splits: ? " );
        display_vector ( splits );
        printf ( "\n" );

        splits.clear();
        splits.push_back ( 0 );
        splits.push_back ( jj );

        printf ( "symmetrygroup2splits: after " );
        display_vector ( splits );
        printf ( "\n" );


    }
    return splits;
}



/// Perform check or reduction using ordering based on delete-one-factor J4 values
int jj45split ( carray_t *array, rowindex_t N, int jj, const colperm_t comb,  const arraydata_t &ad, const OAextend &oaextend, LMC_static_struct_t &tmpStatic, LMCreduction_t &reduction, int verbose=0 )
{
    assert ( jj==5 );
    lmc_t ret;

    //const int verbose=2;
    const int maxjj=40;
    colindex_t comby[maxjj];
    colindex_t perm[maxjj];
    if ( verbose )
        printf ( "-- jj45split --\n" );

    // allocate buffer to hold the values
    std::vector<int> ww ( 5 );

    if ( verbose>=2 ) {
        printf ( "jj45split: init comb " );
        print_perm ( comb, 5 );
    }

    // TODO: replace this block with functions

    /* calculate the J4 values */
    colindex_t lc[4];
    init_perm ( lc, 4 );
    colindex_t lc2[4];
    init_perm ( lc2, 4 );
    for ( size_t i=0; i<5; i++ ) {
        //printf("  comb %d: ", i); print_perm(lc, 4);
        perform_inv_perm ( comb, lc2, 4, lc );

        ww[i]=abs ( jvaluefast ( array, N, 4, lc2 ) );
        if ( verbose>=3 ) {
            printf ( "  comb %d full: val %d: ", ( int ) i, ( int ) ww[i] );
            print_perm ( lc2, 4 );
        }
        next_comb ( lc, 4, 5 );
    }


    //   std::sort ( ww, ww+5 );
    // we reverse the values (to keep ordering matched to the original orderings)
    std::reverse ( ww.begin(), ww.end() );
    if ( verbose>=2 ) {
        printf ( "## jj45split: values (first value is column one removed): " );
        display_vector ( ww );
        printf ( "\n" );
    }

    indexsort s ( 5 );
    s.sort ( ww );
    std::vector<int> wws = s.sorted ( ww );

    // calculate symmetry group of column permutations in the first jj columns
    symmetry_group sg ( wws, true, verbose>=3 );
    if ( verbose>=2 ) {
        // printf("jj45split: values ");  display_vector(ww);   printf("\n");
        printf ( "jj45split: values sorted " );
        display_vector ( wws );
        printf ( "\n  " );
        sg.show();
    }

    // create the sorted column combination
    perform_inv_perm ( comb, comby, jj, s.indices );
    if ( verbose>=2 ) {
        printf ( "initial comb: " );
        print_perm ( comb, jj );
        printf ( " result comb: " );
        print_perm ( comby, jj );
    }
    // make splits
    init_perm ( perm, ad.ncols );
    create_perm_from_comb2<colindex_t> ( perm, comby, jj, ad.ncols );

    if ( verbose>=2 ) {
        printf ( "perm: " );
        print_perm ( perm, ad.ncols );
    }

    dyndata_t dyndata ( ad.N );
    dyndata.col=0; //ad.strength;	// set at col after root
    dyndata.setColperm ( perm, ad.ncols );

    if ( oaextend.getAlgorithm() ==MODE_J5ORDERX || oaextend.getAlgorithm() ==MODE_J5ORDERXFAST ) {
        dyndata.initrowsortl();
    }

    arraydata_t adfix = ad;

    if ( verbose>=4 ) {
        adfix.show_colgroups();
    }

    lmc_t ret1;

    if ( 0 ) {
        // reduce column permutations using the delete-one-factor symmetry group
        adfix.set_colgroups_jj ( sg, 5 );

        if ( verbose>=2 ) {
            adfix.show_colgroups();
        }

        ret = LMCreduce ( array, array, &adfix, &dyndata, &reduction , oaextend );
        ret1=ret;
        if ( verbose )
            printf ( "  ret %d\n", ret );
    }


    std::vector<int> splits = symmetrygroup2splits ( sg, ad.ncols, verbose );

    adfix.set_colgroups ( splits );

    if ( 0 ) {
        printf ( "symmetry indices: %d: ", sg.ngroups );
        display_vector ( sg.gstart );
        printf ( "\n" );
        printf ( "splits: " );
        display_vector ( splits );
        printf ( "\n" );
        adfix.show_colgroups();
    }


    if ( verbose>=2 ) {
        adfix.show_colgroups();
    }

    if ( 0 ) {
        printf ( "entering LMCreduce: " );
        dyndata.show();	// @pte
        print_perm ( dyndata.colperm, adfix.ncols );
        reduction.init_state=COPY;
        printf ( "reduction.init_state %d (COPY %d, INIT %d)\n", reduction.init_state, COPY, INIT );
    }


    ret = LMCreduce ( array, array, &adfix, &dyndata, &reduction , oaextend );

    if ( verbose )
        printf ( "  ret %d\n", ret );

    ret1=ret;
    if ( ret==LMC_LESS && ret1==LMC_MORE ) {
        if ( verbose>=2 ) {
            printf ( "gr\n" );
            //reduction.transformation->show();
        }
    }

    double val=ret1;
    return val;
}

inline double jj452double ( const double *ww )
{
    //	TODO: check multiplication factor
    // the maximum value of the J4 characteristics is N. we create a unique value out of the pair by multiplication
    double val=0;
    for ( size_t i=0; i<6; i++ )
        val = val*64+ww[i];
    return val;

}

/// return value based on J4-J5 ordering
jj45_t jj45val ( carray_t *array, rowindex_t N, int jj, const colperm_t comb, int j5val, int dosort )
{

    double ww[6];

#ifdef OADEBUG
    assert ( N<=MAXROWS );
#endif
    array_t tmpval[MAXROWS];

    std::fill ( tmpval, tmpval+N, 0 );
    fastJupdate ( array, N, jj, comb, tmpval);
    ww[0]=abs ( fastJupdateValue ( N, tmpval ) );

//	colindex_t comb[5];
    for ( size_t i=0; i<5; i++ ) {
//		comb[0]=i;
        int ii = 5-i-1;
        fastJupdate ( array, N, 1, comb+ii, tmpval );
        ww[i+1]=abs ( fastJupdateValue ( N, tmpval ) );
        fastJupdate ( array, N, 1, comb+ii, tmpval ); // OPTIMIZE: eliminate swap back here
    }

    if ( dosort ) {
        std::sort ( ww+1, ww+6,std::greater<int>() );

    }
    //printf("jj45val (dosort %d): j5 %d, j4 seq ", dosort, (int) ww[0]); print_perm(ww+1, 5);

    double val = jj452double ( ww );

    return val;
}

/// return value based on J4-J5 ordering
jj45_t jj45val_orig ( carray_t *array, rowindex_t N, int jj, const colperm_t comb, int j5val=-1, int dosort=1 )
{

    double ww[6];
    if ( j5val==-1 )
        ww[0]=abs ( jvaluefast ( array, N, jj, comb ) );
    else
        ww[0]=j5val;

    colindex_t lc[4];
    init_perm ( lc, 4 );
    colindex_t lc2[4];
    init_perm ( lc2, 4 );
    for ( size_t i=0; i<5; i++ ) {
        //printf("  comb %d: ", i); print_perm(lc, 4);
        perform_inv_perm ( comb, lc2, 4, lc );
        // printf("  comb %zu full: ", i); print_perm(lc2, 4);

        ww[i+1]=abs ( jvaluefast ( array, N, 4, lc2 ) );
        next_comb ( lc, 4, 5 );
    }

    if ( dosort ) {
        std::sort ( ww+1, ww+6,std::greater<int>() );

    }

    double val = jj452double ( ww );

    return val;
}

#ifdef LMCSTATS
// this code is not thread safe!
static long nnn=0;
static long nng=0;
void lmc_stats()
{
    printf ( "nng %ld, nnn %ld, fraction %.3f\n", nng, nnn, double ( nng ) /double ( nnn ) );
}
#endif

lmc_t LMCcheckj5 ( array_link const &al, arraydata_t const &adin, LMCreduction_t &reduction, const OAextend &oaextend, int hack )
{
    const int dverbose=0;
    LMC_static_struct_t &tmpStatic = reduction.getStaticReference();
    //tmpStatic.setRef("LMCcheckj5");

    const int jj=5;
#ifdef OACHECK
    const int maxjj = 10;
    assert ( jj<maxjj );
#endif
    // perform basic integrity checks
    if ( !hack )
        assert ( adin.ncolgroups==1 );

    assert ( adin.s[0]==2 );
    assert ( jj>=adin.strength );
    myassert ( adin.strength>=3, "LMCcheckj5: strength should be >=3!" );

    if ( reduction.mode>=OA_REDUCE && !hack ) {
        printfd ( "error: reduction mode not implemented for J5 ordering!!!\n" );
    }

    if ( reduction.init_state==COPY ) {
        reduction.setArray ( al );
        reduction.init_state=COPY;
    }
    if ( reduction.init_state==INIT_STATE_INVALID ) {
        // TODO: remove this code
        printf ( "LMCcheckj5: reduction.init_state is INVALID\n" );
        throw;
    }


    /*
    log_print ( DEBUG+1, "LMCcheckj5: reduction->mode %d\n", reduction.mode );
    if ( ! ( reduction.mode==LMC_REDUCE_INIT ) && !hack ) {
        log_print ( DEBUG+1, "   copying array for init\n", reduction.mode );
        reduction.setArray(al);
    } */

    // NOTE: the reduction.array is reduced, so reduction.array is a good value to start from.
    // NOTE: for arrays not in root form we set update the reduction.array to make sure that the root is set (this part is not checked!)
    // IMPROVEMENT: only do this if the reduction.array is not in root form (automatic REDUCE_INIT, but this also works for LMC_TEST)
    bool rootform = check_root_form ( reduction.array, adin );

    if ( ! rootform ) {
        log_print ( SYSTEM, "error: LMC test or LMC reduction for arrays not in root form needs special initialization !\n" );
        print_array ( reduction.array, adin.N, adin.ncols );
    }



    if ( al.n_columns<jj ) {
        // fall back to original method
        dyndata_t dyndata ( adin.N );
        lmc_t ret = LMCreduce ( al.array, al.array, &adin, &dyndata, &reduction , oaextend );
        return ret;
    }

    arraydata_t ad ( adin );
    // split column group!
    // the first 5 columns are selected in this function (and if necessary we loop over all permutations). the remaining columns are then processed with the original LMC algorithm
    ad.complete_arraydata_splitn ( jj );

    tmpStatic.update ( &ad );

    colindex_t *firstcolcomb = new_comb_init<colindex_t> ( jj );
    int nrootcombs = ncombs ( jj, ad.strength );
    colindex_t *combroot = new_comb_init<colindex_t> ( ad.strength );	//
    colindex_t combx[maxjj];	// final combination for root selection + extra columns
    colindex_t comby[maxjj];

    carray_t *array = al.array;
    if ( hack )
        array = reduction.array;

    lmc_t ret = LMC_EQUAL;

    /* loop over all possible column combinations for the first jj columns */
    int jbase = abs ( jvaluefast ( array, ad.N, jj, firstcolcomb ) );
    jj45_t wbase = jj45val ( array, ad.N, jj, firstcolcomb, -1, 0 );

    /* initialize variables outside of loop */
    colperm_t perm = new_perm<colindex_t> ( ad.ncols );
    colperm_t pp = new_perm_init<colindex_t> ( jj );

    const int orig = oaextend.j5structure ==J5_ORIGINAL;

    int ncolsfirst=adin.colgroupsize[0]; // = ad.ncols
    int nc = ncombs ( ncolsfirst, jj );	// number of combinations to select the first jj columns

    if ( dverbose ) {
        printf ( "LMCcheckj5: selected ncolsfirst %d, nc %d, jbase %d, wbase %f\n", ncolsfirst, nc, jbase, wbase );
    }

    if ( hack ) {
        if ( hack>=2 ) {
            printf ( "LMCcheckj5: selected ncolsfirst %d, nc %d, jbase %d, wbase %f\n", ncolsfirst, nc, jbase, wbase );
            print_array ( "array:", array,  al.n_rows,5 );
        }
        array=al.array;
        if ( hack>=2 ) {
            printf ( "--\n" );
            print_array ( "al.array:\n", array,  al.n_rows,5 );
        }
//        nc=1; //
    }

    // loop over all possible combinations for the first jj columns
    for ( int i=0; i<nc; i++ ) {

        // determine the J-characteristic for the selected combination
        int j5val=abs ( jvaluefast ( array, ad.N, jj, firstcolcomb ) ); // TODO: this can be made faster by caching sub results...

        if ( dverbose>=2 ) {
            printf ( "   LMCcheckj5 column comb %d/%d (jj=%d): j5val %d, jbase %d\n", i, nc, jj, j5val, jbase );
        }

        // check order based on J5
        if ( j5val ORDER_J5_SMALLER jbase ) {
            ret=LMC_LESS;
            //reduction->lastcol=col;
            reduction.updateLastCol ( 4 );

            if ( dverbose ) {
                printfd ( "LMCcheckj5: ret %d: j5val %d, jbase %d\n", ret, j5val,jbase );
            }
            break;
        } else if ( j5val ORDER_J5_GREATER jbase ) {
            // this combination can only lead to LMC more
            next_combination ( firstcolcomb, jj, ad.ncols );	// increase combination
            continue;
        }

        if ( dverbose>=2 ) {
            printf ( "   LMCcheckj5 column comb %d/%d : going to jjsplit\n", i, nc );
        }

        const int dverbose=0;
        // check order based on J4-J5 value
        if ( 1 ) {
#ifdef LMCSTATS
            nnn++;
#endif
            jj45_t w = jj45val ( array, ad.N, jj, firstcolcomb, j5val );
            if ( w ORDER_J45_SMALLER wbase ) {
                ret=LMC_LESS;
                if ( dverbose )
                    printfd ( "LMCcheckj5: ret %d: w %ld, wbase %ld\n", ( int ) ret, ( long ) w, ( long ) wbase );
                reduction.updateLastCol ( 4 );
                //		if reduction.doBreak(ret)
                break;
            } else if ( w ORDER_J45_GREATER wbase ) {
                // this combination can only lead to LMC more
                next_combination ( firstcolcomb, jj, ad.ncols );	// increase combination
                continue;
            }
#ifdef LMCSTATS
            nng++;
#endif
        }

        //    reduction.storeColumnPermutation(firstcolcomb, jj);

        if ( dverbose>=2 ) {
            printf ( "   LMCcheckj5 column comb %d/%d : going to jjsplit\n", i, nc );
        }


        if ( orig ) {
            // do classic check
            init_comb ( combroot, ad.strength, jj );
            for ( int kk=0; kk<nrootcombs; kk++ ) {
                // prepare for all root combinations

                // convert combination to full permutation
                create_perm_from_comb ( combx, combroot, ad.strength, jj, 0 );
                perform_inv_perm ( firstcolcomb, comby, jj, combx );
                init_perm ( perm, ad.ncols );
                create_perm_from_comb2<colindex_t> ( perm, comby, jj, ad.ncols );

                dyndata_t dyndata ( ad.N );
                dyndata.col=ad.strength;	// set at col after root
                copy_perm ( perm, dyndata.colperm, ad.ncols );
                // copy_array ( array, reduction.array, ad.N, ad.ncols );


                ret = LMCreduce_root_level_perm_full ( array, &ad, &dyndata,&reduction, oaextend, tmpStatic );

                if ( ret==LMC_LESS ) {
                    printfd ( "LMCcheckj5: note: code unchecked...\n" );
                    //    reduction.updateLastCol(4);
                    break;
                }

                next_combination ( combroot, ad.strength, jj );	// increase combination
            }
        }
        int retx =ret;
        if ( !orig ) {

            retx= jj45split ( array, ad.N, jj, firstcolcomb, ad, oaextend, tmpStatic, reduction );

            if ( retx==LMC_LESS ) {
                // NOTE: is this needed?
                //reduction.updateLastCol(4);
            }

        }
        ret= ( lmc_t ) ( int ) retx;

#ifdef OADEBUG
        if ( 0 ) {
            if ( ret!= retx ) {
                printfd ( "##  final ret x %d (new method %d) ##\n", ret, ( int ) retx );
                printf ( "   comb: " );
                print_perm ( firstcolcomb, jj );
                printf ( "   colperm: " );
                print_perm ( perm, ad.ncols );
                setloglevel ( DEBUG );

                double retx = jj45split ( array, ad.N, jj, firstcolcomb, ad, oaextend, tmpStatic, reduction, 2 );
                throw;
            }
        }
#endif

        if ( ret==LMC_LESS ) {
            // printf("LMCcheckj5: found LMC_LESS\n");
            break;
        }

        next_combination ( firstcolcomb, jj, ncolsfirst );	// increase combination
    }
    delete_perm ( perm );
    delete_perm ( pp );
    delete_comb ( combroot );
    delete_comb ( firstcolcomb );

    return ret;
}

lmc_t LMCcheckLex ( array_link const &al, arraydata_t const &ad )
{
    int N = al.n_rows;

    if ( al.n_columns!=ad.ncols ) {
        printf ( "error: number of columns of matrix inconsistent!\n" );
        return LMC_NONSENSE;
    }

    dyndata_t dynd = dyndata_t ( N );
    lmc_t result = LMC_NONSENSE;
    arraydata_t adx = arraydata_t ( &ad, al.n_columns );
    LMCreduction_t *reduction = new LMCreduction_t ( &adx );
    OAextend oaextend;
    result = LMCreduce ( al.array, al.array, &ad, &dynd, reduction , oaextend );

    return result;
}


lmc_t LMCcheckSymmetryMethod ( const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction,  LMCreduction_t &reductionsub, int dverbose )
{

    int special=0;
    double t0=get_time_ms();
    double dt=0;

    if ( dverbose ) {
        printf ( "### LMCcheckSymmetryMethod: ncols %d, dverbose %d\n", ad.ncols, dverbose );
    }

    lmc_t r = LMC_EQUAL;

    // TODO: add assert to type of adata

    std::vector<colindex_t> pinit = permutation<colindex_t> ( ad.ncols );

    int newcol=ad.ncols-1;
    arraydata_t adfix = ad;
    dyndata_t dyndata ( ad.N );

    symmdata sd ( al );
    if ( dverbose>=2 ) {
        al.showarray();
        sd.show();
    }

    LMC_static_struct_t *tmpStatic = & reduction.getStaticReference();

//    LMC_static_struct_t *tmpStatic = getGlobalStatic();
    tmpStatic->update ( &ad );

    reduction.setArray ( al );


    // loop over all symmetries
    for ( int i=ad.strength; i<ad.ncols+1-1; i++ ) {
        // continue;
        if ( dverbose ) {
            dt =get_time_ms()-t0;
            printf ( "LMCcheckSymmetryMethod: symmetryset i %d: %ld, dt %.3f [ms]\n", i, ( long ) reductionsub.symms.symmetries[i].size(), 1e3*dt );
        }
        // split column groups
        std::vector<int> splits; // = symmetrygroup2splits(sg, ad.ncols, verbose);
        for ( int k=0; k<=i+1; k++ )
            splits.push_back ( k );
        adfix.set_colgroups ( splits );	// should not be needed

        if ( ad.ncols==0 )
            printf ( "ad.ncols==0!\n" );
        larray<colindex_t> ww ( ad.ncols );
        larray<colindex_t> wtmp ( ad.ncols );
        for ( int j=0; j< ( int ) reductionsub.symms.symmetries[i].size(); j++ ) {
            cprintf ( dverbose>=2, "LMCcheckSymmetryMethod: testing symmetry %d: newcol %d\n", j, newcol );
//            std::vector<colindex_t> w = reductionsub.symmetries[i][j].colperm;  w.push_back(newcol);
            //larray<colindex_t> w = reductionsub.symmetries[i][j].colperm.addelement(newcol);

            //printf("combadd2perm: newcol %d, ad.cols %d, wtmp.size() %d, colperm ", newcol, ad.ncols, wtmp.size() ); print_perm(reductionsub.symms.symmetries[i][j].colperm);
            combadd2perm ( reductionsub.symms.symmetries[i][j].colperm, newcol, ad.ncols, ww, wtmp ); // expensive
            //larray<colindex_t> ww = comb2perm<colindex_t>(w, ad.ncols);

            // cprintf(dverbose>=2, "newcol %d, i %d. w.size() %d\n", newcol, i, w.size());

            // initialize dyndata with w
            // dyndata.reset();
            dyndata.col=ad.strength;	// set at col after root
            dyndata.col=i; // trick
            //dyndata.col=i+1; //

            dyndata.setColperm ( ww );

            // initialize dyndata with rowperm and init level perms
            dyndata.initsymmetry ( reductionsub.symms.symmetries[i][j], sd, i );
            // note: above call  is expensive, merge with initrowsortl call

            if ( dverbose>=2 ) {
                printf ( "  dyndata perm: " );
                print_perm ( dyndata.colperm, ad.ncols );
            }

            //cprintf(dverbose>=2, "  calling LMCreduce\n");

            if ( dverbose ) {
                //printf("LMCcheckSymmetryMethod: calling LMCreduce_non_root...\n");
            }
            /* pass to non-root stage */
            if ( oaextend.getAlgorithm() ==MODE_LMC_SYMMETRY || oaextend.getAlgorithm() ==MODE_J5ORDERXFAST || oaextend.getAlgorithm() ==MODE_LMC_2LEVEL ) {
                // TODO: make call to LMCreduce_non_root_2level more generic
                dyndata.initrowsortl();
                //dyndata.ro
                r= LMCreduce_non_root_2level ( al.array, &adfix, &dyndata,& reduction, oaextend, *tmpStatic );

            } else {
                printf ( "oaextend problem: " );
                oaextend.info();
                throw;
                r= LMCreduce_non_root ( al.array, &adfix, &dyndata,& reduction, oaextend, *tmpStatic );
            }

            if ( dverbose>=2 ) {
                printf ( "  i %d: lmc_t returned %d, splits: perm ", i, int ( r ) );
                print_perm ( dyndata.colperm, ad.ncols );
                adfix.show_colgroups();
                if ( special ) {
                    // exit(0);
                }
            }

            if ( r==LMC_LESS ) {
                if ( dverbose ) {
                    printf ( "LMCcheckSymmetryMethod: return LMC_LESS: perm " );
                    larray<colindex_t> w = reductionsub.symms.symmetries[i][j].colperm;
                    w=w.addelement ( newcol ); // w.push_back(newcol);
                    print_perm ( w );
                    adfix.show();
                    dyndata.show();
                }
                return r;
            }

        }
    }

    cprintf ( dverbose>=1, "  LMCcheckSymmetryMethod: symmetries done: time: %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );

    // loop over all colcomb pairs
    for ( int i=ad.strength; i<ad.ncols+1-1; i++ ) {
        // break;
        if ( dverbose ) {
            dt =get_time_ms()-t0;
            printf ( "LMCcheckSymmetryMethod: colcombs i %d: %ld, dt %.1f [ms]\n", i, ( long ) reductionsub.symms.colcombs[i].size(), 1e3* dt );
//    printf("original check: %d, time: %.3f [ms]\n", (int)r, 1e3*dt);

        }

        // split column groups
        std::vector<int> splits;
        splits.push_back ( 0 );
        splits.push_back ( 5 );

        reduction.setArray ( al );
        for ( int j=0; j< ( int ) reductionsub.symms.colcombs[i].size(); j++ ) {
            std::vector<int> w = reductionsub.symms.colcombs[i][j];
            w.push_back ( newcol );
            // printf("comb2perm: i %d: al %d %d: reductionsub.symms %d, comb size %zu, ad.ncols %d: ", i, al.n_rows, al.n_columns, reductionsub.symms.ncols, w.size(), ad.ncols ); print_perm(w);
            std::vector<colindex_t> ww = comb2perm<colindex_t> ( w, ad.ncols );

            // initialize dyndata with w
            dyndata.reset();
            dyndata.col=0;

            dyndata.setColperm ( ww );
            if ( dverbose>=2 ) {
                printf ( "dyndata perm: " );
                print_perm ( dyndata.colperm, ad.ncols );
            }

            adfix.set_colgroups ( splits );
            if ( dverbose>=2 ) {
                printf ( "LMCcheckSymmetryMethod: adfix.show_colgroups()\n" );
                adfix.show_colgroups();
            }

            // NOTE: allow non-proper initialization of array
            //reduction.mode=LMC_REDUCE_INIT;
            reduction.mode=OA_TEST;

            array_link alx = al.selectColumns ( ww );

            cprintf ( dverbose>=2, "  LMCcheckSymmetryMethod: colcombs check : starting...\n" );

            reduction.setArray ( al );
            //alx=al;
            r = LMCcheckj5 ( alx, adfix, reduction, oaextend, 1 );


            if ( dverbose>=2 ) {
                printf ( "i %d: lmc_t returned %d, splits: perm ", i, int ( r ) );
                print_perm ( dyndata.colperm, ad.ncols );
                adfix.show_colgroups();

                if ( special ) {
                    // exit(0);
                }
            }

            if ( r==LMC_LESS ) {
                if ( dverbose ) {
                    printf ( "LMCcheckSymmetryMethod: return LMC_LESS: perm " );
                    print_perm ( w );
                    adfix.show();
                    dyndata.show();
                }
                return r;
            }

        }
    }
    cprintf ( dverbose>=1, "  LMCcheckSymmetryMethod: colcombs done: time: %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );

    // loop over all colperm pairs
    for ( int i=ad.strength; i<ad.ncols+1-1; i++ ) {
        if ( dverbose>=2 ) {
            dt =get_time_ms()-t0;
            printf ( "LMCcheckSymmetryMethod: colperms i %d: %ld, dt %.1f [ms]\n", i, ( long ) reductionsub.symms.colperms[i].size(), 1e3*dt );
//    printf("original check: %d, time: %.3f [ms]\n", (int)r, 1e3*dt);

        }
        for ( int j=0; j< ( int ) reductionsub.symms.colperms[i].size(); j++ ) {
            std::vector<int> w = reductionsub.symms.colperms[i][j];
            w.push_back ( newcol );
            //std::vector<colindex_t> ww = perform_perm(pinit, w);
            std::vector<colindex_t> ww = comb2perm<colindex_t> ( w, ad.ncols );

            if ( dverbose>=2 || i==23 ) {
                printf ( "processing perm " );
                print_perm ( w );
                printf ( "	--> " );
                print_perm ( ww );
                printf ( "	--> " );

            }


            // initialize dyndata with w
            dyndata.reset();
            dyndata.col=ad.strength;	// set at col after root
            dyndata.col=0; // ?
            // adfix.strength=1; adfix.calcoaindex(ad.strength); // ad.strength;	// HACK -> should be strength!

            dyndata.setColperm ( ww );
            if ( dverbose>=2 ) {
                printf ( "dyndata perm: " );
                print_perm ( dyndata.colperm, ad.ncols );
            }

            // split column groups
            std::vector<int> splits; // = symmetrygroup2splits(sg, ad.ncols, verbose);
            for ( int k=0; k<=i+1; k++ )
                splits.push_back ( k );
            adfix.set_colgroups ( splits );
            // call original lmc check -> use function ???

            r = LMCreduce ( al.array, al.array, &adfix, &dyndata, &reduction , oaextend );


            if ( dverbose>=2 ) {
                printf ( "i %d: lmc_t returned %d, splits: perm ", i, int ( r ) );
                print_perm ( dyndata.colperm, ad.ncols );
                adfix.show_colgroups();

                if ( special ) {
                    // exit(0);
                }
            }

            if ( r==LMC_LESS ) {
                if ( dverbose ) {
                    printf ( "LMCcheckSymmetryMethod: return LMC_LESS: perm " );
                    print_perm ( w );
                    adfix.show();
                    dyndata.show();
                }
                return r;
            }

        }
    }
    cprintf ( dverbose>=1, "  LMCcheckSymmetryMethod: colperms done: time: %.1f [ms]\n", 1e3* ( get_time_ms()-t0 ) );


    releaseGlobalStatic ( tmpStatic );

    return LMC_MORE;
}


LMCreduction_t calculateSymmetryGroups ( const array_link &al, const arraydata_t &adata, const OAextend &oaextend, int dverbose, int hack )
{

// dverbose=1;
    if ( dverbose ) {
        printf ( "calculateSymmetryGroups: array " );
        al.show();
    }
    // pre-compute
    array_link alsub = al;
    arraydata_t adatasub ( &adata, al.n_columns );
    //array_link alsub = al.deleteColumn(-1);  arraydata_t adatasub(&adata, alsub.n_columns);
    LMCreduction_t reductionsub ( &adatasub );
    reductionsub.init_state=COPY;

    reductionsub.symms.store=2;	// store both symmetries and combintions

    // add elements at strength level
    int n = adatasub.ncols;

    int baselevel = adata.strength;
    if ( oaextend.getAlgorithm() ==MODE_J5ORDERX || oaextend.getAlgorithm() ==MODE_J5ORDERXFAST || oaextend.getAlgorithm() ==MODE_LMC_SYMMETRY ) {
        if ( dverbose )
            printf ( "calculateSymmetryGroups: setting baselevel to that of J45\n" );
        baselevel=4;

        std::vector<colindex_t> combp = permutation<colindex_t> ( baselevel );
        int nc = ncombs ( adatasub.ncols, baselevel );
        int np = factorial<int> ( baselevel );
        np=1;
        if ( dverbose )
            printf ( "calculateSymmetryGroups: colcombs: ncombs %d, np %d\n", nc, np );
        for ( int i=0; i<nc; i++ ) {
            std::vector<colindex_t> px = permutation<colindex_t> ( baselevel );
            for ( int j=0; j<np; j++ ) {

                std::vector<colindex_t>  pv= perform_perm ( combp, px );
                if ( dverbose>=2 )
                    printf ( "calculateSymmetryGroups: colcombs: store combination i %d, j %d\n", i, j );

                reductionsub.symms.storeColumnCombination ( pv );
                next_perm ( px );
            }
            //printf("comb: "); print_perm(combp);
            next_comb ( combp, baselevel, n );
        }


    } else {
        std::vector<colindex_t> combp = permutation<colindex_t> ( baselevel );
        int nc = ncombs ( adatasub.ncols, baselevel );
        int np = factorial<int> ( baselevel );
        //printf("calculateColpermsGroup: ncombs %d, np %d\n", nc, np);
        for ( int i=0; i<nc; i++ ) {
            std::vector<colindex_t> px = permutation<colindex_t> ( baselevel );
            for ( int j=0; j<np; j++ ) {

                // std::vector<colindex_t>  pv= array2vector(combp, adata.strength);
                std::vector<colindex_t>  pv= perform_perm ( combp, px );
                reductionsub.symms.colperms[adata.strength].push_back ( pv );
                next_perm ( px );
            }
            //printf("comb: "); print_perm(combp);
            next_comb ( combp, baselevel, n );
        }
    }


    if ( dverbose>=2 ) {
        reductionsub.symms.showColcombs();
    }

    reductionsub.symms.ncols=al.n_columns;

    if ( alsub.n_columns==adata.strength )
        return reductionsub;

    lmc_t r;
    if ( 0 ) {
        LMCreduction_t reduction ( &adata );
        reductionsub.symms.store=0;
        r = LMCcheck ( al, adata, oaextend, reduction ) ;
        reductionsub.symms.store=0;
        r = LMCcheck ( alsub, adatasub, oaextend, reductionsub ) ;
    }

    //  LMCreduction_t reductionsubG(&adata);
    //r = LMCcheck(al, adata, oaextend, reductionsubG) ;
    //alsub = al.selectFirstColumns(al.n_columns-1);

    //  printf("calculateSymmetryGroups:  reductionsub.symms.store %d\n",  reductionsub.symms.store);

    OAextend x = oaextend;
    // x.setAlgorithm(MODE_J5ORDERX, &adatasub);
    x.setAlgorithm ( MODE_J5ORDERXFAST, &adatasub );
    r = LMCcheck ( alsub, adatasub, x, reductionsub ) ;
    if ( dverbose )
        printf ( "calculateSymmetryGroups: LMCcheck complete...\n" );
    // printf("calculateSymmetryGroups:  reductionsub.symms.store %d\n",  reductionsub.symms.store);
// reductionsub.symms.show();

    //reductionsub=LMCreduction_t(&adatasub);
    //r = LMCcheck(alsub, adatasub, oaextend, reductionsub) ;
    if ( dverbose )
        printf ( "calculateColpermsGroup: sub lmc_t %d\n", ( int ) r );
    if ( dverbose>=2 || 0 ) {
        reductionsub.symms.show ( 2 );
    }

    reductionsub.symms.makeColpermsUnique ( dverbose );


    if ( oaextend.getAlgorithm() ==MODE_J5ORDERX || oaextend.getAlgorithm() ==MODE_J5ORDERXFAST || oaextend.getAlgorithm() ==MODE_LMC_SYMMETRY ) {
        //reductionsub.symms.showColcombs();
        //reductionsub.symms.showSymmetries();


        for ( int i=0; i<5; i++ ) {
            reductionsub.symms.symmetries[i].clear();
        }
    }
    return reductionsub;

}

lmc_t LMCcheckOriginal ( const array_link &al ) {
    assert ( al.is2level() );
    arraydata_t ad = arraylink2arraydata ( al );
    int strength = al.strength();

    OAextend oaextend ( ad );
    LMCreduction_t reduction ( &ad );
	reduction.init_state=INIT;
			copy_array(al.array, reduction.array, ad.N, ad.ncols);

	            int changed = check_root_update ( al.array, ad, reduction.array );

	//            reduction->setArray ( array, ad->N, ad->ncols );

	
    return LMCcheck ( al.array, ad,  oaextend, reduction );
}

lmc_t LMCcheck ( const array_t * array, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction )
{
    lmc_t lmc = LMC_NONSENSE;
    dyndata_t dynd = dyndata_t ( ad.N );

    //printfd("reduction.init_state %d\n", reduction.init_state);

    if ( reduction.init_state==INIT_STATE_INVALID ) {
        printf ( "LMCcheck: reduction.init_state is INVALID, please set it to COPY or INIT!\n" );
        //oaextend.info();
        throw;
    }
    if ( reduction.init_state==SETROOT ) {
        log_print ( DEBUG, "LMCcheck: reduction.init_state is SETROOT\n" );

        array_link alp = rootPlus ( ad );
        reduction.setArray ( alp );
        reduction.init_state=INIT;
        //printf("LMCcheck: updated reduction.array to\n");
        //reduction.getArray().showarray();
    }
    if ( reduction.init_state==COPY ) {
        reduction.setArray ( array, ad.N, ad.ncols );
        //reduction->init_state=COPY;

        //printfd("LMCreduce: copied array N %d, k %d, reduction->init_state %d\n", ad.N, ad.ncols, reduction.init_state);
        //print_array(reduction.array, ad.N, ad.ncols);

    }

    if ( 0 ) {
        printfd ( "LMCcheck %s line %d: switching algorithm: init_state %d\n", __FUNCTION__, __LINE__, reduction.init_state );
        oaextend.info();
    }

    //printfd("oaextend.getAlgorithm() %d, MODE_ORIGINAL %d\n", oaextend.getAlgorithm(),MODE_ORIGINAL );
    switch ( oaextend.getAlgorithm() ) {
    case MODE_ORIGINAL: {
        lmc = LMCreduce ( array, array, &ad, &dynd, &reduction , oaextend );

    }
    break;
    case  MODE_J4: {
        array_link al ( array, ad.N, ad.ncols, -20 );
        // if(reduction.init_state!=COPY) {
        //     printf("LMCcheck: warning: MODE_J4: reduction.init_state!=COPY (reduction.init_state %d)\n", reduction.init_state);
        // }
        // copy_array ( array, reduction.array, ad.N, ad.ncols );
        //reduction.sd = symmdataPointer(new symmdata(al, 1) );
        //printf("debug: reduction->mode %d\n", reduction.mode);
        lmc = LMCcheckj4 ( al, ad, reduction , oaextend );
    }
    break;
    case  MODE_J5ORDER: {
        printf ( "LMCcheck: algorithm MODE_J5ORDER: code path untested\n" );
        lmc = LMCreduce ( array, array, &ad, &dynd, &reduction, oaextend );
    }
    break;
    case MODE_LMC_2LEVEL: {
        array_link al ( array, ad.N, ad.ncols, -20 );
        reduction.updateSDpointer ( al, true );
        lmc = LMCreduce ( array, array, &ad, &dynd, &reduction , oaextend );
    }
    break;
    case MODE_LMC_SYMMETRY: {
        int dverbose=0;

        // testing code!
        array_link al ( array, ad.N, ad.ncols );
        OAextend x = oaextend;
        arraydata_t adx ( ad );
        x.setAlgorithm ( MODE_J5ORDERX, &adx );
        x.setAlgorithm ( MODE_J5ORDERXFAST, &adx );	 // TODO: work this out
        //x.setAlgorithm(MODE_ORIGINAL, &adx);
        // LMCreduction_t reductionsub2 = calculateSymmetryGroups( al, adx,  x, 0, 1);
        if ( dverbose ) {
            printf ( "LMCcheck: MODE_LMC_SYMMETRY: calculateSymmetryGroups " );
        }

        LMCreduction_t reductionsub ( &adx );
        if ( reduction.symms.valid() == 0 || 0 ) {
            printf ( "LMCcheck: MODE_LMC_SYMMETRY calculating new symmetry group (al %d %d)\n", al.n_rows, al.n_columns );
            double t0=get_time_ms();
            reductionsub = calculateSymmetryGroups ( al.deleteColumn ( -1 ), adx,  x, 0 );
            reductionsub.symms.show();
            printf ( " dt symmetry group: %.1f [ms]\n" , 1e3* ( get_time_ms()-t0 ) );
        } else {
            // symmetryPointer;
            // printf("LMCcheck: MODE_LMC_SYMMETRY: skipping calculation of symmetries... (symms.ncols %d, al.ncols %d, symms.symmetries[4].size() %zu, symms.symmetries[5].size() %zu)\n", reduction.symms.ncols, al.n_columns, reduction.symms.symmetries[4].size(), reduction.symms.symmetries[5].size()  );
            //reductionsub.symms.show();
            //vector<symmetryset> vs = reduction.symmetries.get();

            // NOTE: copy symmetries, but extend with new columns as well!
            reductionsub.symms = reduction.symms;
            //symmetries = reduction.symmetries;
        }
        if ( dverbose ) {
            printf ( "LMCcheck: MODE_LMC_SYMMETRY: testing: strength %d, al ", ad.strength );
            al.show();
            reductionsub.symms.showColperms ( 1 );
        }
        // reduction.clearSymmetries();

        OAextend xx =oaextend;
        xx.setAlgorithm ( MODE_J5ORDERXFAST );
        reduction.updateSDpointer ( al );

        lmc = LMCcheckSymmetryMethod ( al, ad, xx, reduction, reductionsub, dverbose ) ;

        if ( 0 ) {
            lmc_t lmcorig=lmc;
            LMCreduction_t reductionx=reduction;
            lmcorig= LMCreduce ( array, array, &ad, &dynd, &reductionx, oaextend );

            if ( lmcorig!=lmc ) {
                printf ( "MODE_LMC_SYMMETRY: difference! lmcorig %d, lmc new %d\n", ( int ) lmcorig, ( int ) lmc );
            }
        }


    }
    break;
    case  MODE_J5ORDERX: {
        array_link al ( array, ad.N, ad.ncols, -20 );
        copy_array ( array, reduction.array, ad.N, ad.ncols );
            lmc = LMCcheckj5 ( al, ad, reduction, oaextend );
    }
    break;
    case  MODE_J5ORDERXFAST: {
        array_link al ( array, ad.N, ad.ncols, -20 );
        copy_array ( array, reduction.array, ad.N, ad.ncols );
        // reduction.sd = symmdataPointer(new symmdata(al, 1) );

        //double t0 = get_time_ms();

        OAextend oaextend2 = oaextend;
        oaextend2.setAlgorithm ( MODE_LMC_2LEVEL );
        reduction.updateSDpointer ( al, false );

        lmc = LMCcheckj5 ( al, ad, reduction, oaextend2 );
    }
    break;

    default: {
        printf ( "file %s: LMCcheck: line %d: error: no such mode %d\n", __FILE__, __LINE__, oaextend.getAlgorithm() );
        printf ( "  algorithm list: %s\n", algorithm_t_list().c_str() );
    }
    break;
    }

    return lmc;
}

lmc_t LMCreduction_train ( const array_link &al, const arraydata_t* ad, LMCreduction_t *reduction, const OAextend &oaextend )
{
    //TODO: pass oaextend options to reduction train
    dyndata_t dyndata ( ad->N, 0 );
    // OAextend oaextend;
    return LMCreduction_train ( al.array, ad, &dyndata, reduction , oaextend ) ;
}

/**
* @brief Perform an LMC reduction in several steps.
* @param original
* @param array
* @param ad
* @param dyndata
* @param reduction
* @return
*/
lmc_t LMCreduction_train ( const array_t* original, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction,const OAextend &oaextend )
{
    // OPTIMIZE: analyze this code for efficiency

    lmc_t ret = LMC_NONSENSE;
    int n;
    colindex_t *cols;

    //printf("ad->ncols %d\n", ad->ncols);
    switch ( ad->ncols ) {
    case 15:
        /* tight */
        n = std::max ( ( int ) floor ( ( double ) ( ad->ncols-7 ) /2 ), 1 );
        cols = new colindex_t[n];
        for ( int z=0; z<n; z++ ) {
            cols[z] = ad->ncols - ( n-z-1 );
        }
        break;
    default
            :
        /* default mode */
        n = std::max ( ( int ) floor ( ( double ) ( ad->ncols-7 ) /2 ), 1 );
        cols = new colindex_t[n];
        for ( int z=0; z<n; z++ ) {
            cols[z] = ad->ncols - ( n-z-1 ) *2 ;
        }
        break;
    }

    logstream ( DEBUG ) << "   cols selected: ";
    if ( checkloglevel ( DEBUG ) )
        print_perm ( cols,n );

    array_t *atmp = create_array ( ad );


    LMCreduction_t reductiontmp ( ad );
    reduction->mode = OA_REDUCE;
    reductiontmp.mode = OA_REDUCE;

    copy_array ( original, reduction->array, ad->N, ad->ncols );
    copy_array ( original, reductiontmp.array, ad->N, ad->ncols );
    copy_array ( original, atmp, ad->N, ad->ncols );
    //*(reduction->transformation) = (*(reductiontmp.transformation))*(*(reduction->transformation));

    const  int verbose=0;

    std::vector<array_link> aldbg ( 100 );

    for ( colindex_t z=0; z<n; z++ ) {
        copy_array ( reductiontmp.array, atmp, ad->N, ad->ncols );
        //printf("z=%d: init array: \n", z); print_array(atmp, ad->N, ad->ncols );

        logstream ( NORMAL ) << printfstring ( "   LMCreduction_train: stage %d/%d", z, n ) << printfstring ( " (current column %d, number of columns %d)", cols[z], cols[n-1] ) << endl;
        //dyndata->show();

        reductiontmp.maxdepth = cols[z];
        reductiontmp.reset();
        reductiontmp.transformation->reset();
        reductiontmp.mode = OA_REDUCE;
        reductiontmp.nred=0;

        //printf("before LMCreduce: "); reductiontmp.show();
        ret = LMCreduce ( atmp, atmp, ad, dyndata, &reductiontmp , oaextend ) ;

        if ( 0 ) {
            // loop check
            array_link altmp ( atmp, ad->N, ad->ncols );
            array_link alreductiontmp ( reductiontmp.array, ad->N, ad->ncols );
            array_link altr = reductiontmp.transformation->apply ( altmp );

            printf ( "#### loop check #### z=%d: ret %d, atmp(warp)==reductiontmp.array %d\n", z, ret, altr==alreductiontmp );
            if ( ! ( altr==alreductiontmp ) ) {
                printf ( "input atmp:\n" );
                altmp.showarray();
                altr.showarray();
                alreductiontmp.showarray();
                reductiontmp.show();
                return LMC_NONSENSE;
            }
        }


        if ( verbose>=3 ) {
            printf ( "### LMCreduction_train: concatenate the transformations (z=%d, column %d)\n", z, cols[z] );
            printf ( "  accum -> \n" );
            reduction->transformation->show();
            printf ( "  tmp -> \n" );
            reductiontmp.transformation->show();
        }

        * ( reduction->transformation ) = ( * ( reductiontmp.transformation ) ) * ( * ( reduction->transformation ) );

    }
    //printf("LMCreduction_train: final iteration\n");
    /* final iteration */
    copy_array ( reductiontmp.array, atmp, ad->N, ad->ncols );
    reduction->maxdepth = -1;
    reductiontmp.transformation->reset();

    //array_link alxx(atmp, ad->N, ad->ncols); printf("  atmp\n"); alxx.showarray();
    ret = LMCreduce ( atmp, atmp, ad, dyndata, &reductiontmp , oaextend ) ;

    if ( verbose>=2 ) {
        printf ( "### LMCreduction_train: concatenate the transformations (final round)\n" );
        printf ( "  accum -> \n" );
        reduction->transformation->show();
        printf ( "  tmp -> \n" );
        reductiontmp.transformation->show();
    }
    * ( reduction->transformation ) = ( * ( reductiontmp.transformation ) ) * ( * ( reduction->transformation ) );
    //	*(reduction->transformation) = (*(reduction->transformation))*(*(reductiontmp.transformation));

    if ( verbose>=2 ) {
        printf ( "  final -> \n" );
        reduction->transformation->show();
    }

    copy_array ( reductiontmp.array, reduction->array, ad->N, ad->ncols );

    destroy_array ( atmp );
    delete [] cols;
    return ret;
}


/// full reduction, no root-trick
lmc_t LMCreduceFull ( carray_t* original, const array_t *array, const arraydata_t* adx, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend, LMC_static_struct_t &tmpStatic )
{
    arraydata_t *ad= new arraydata_t ( *adx ); // NOTE: is this needed?
    ad->oaindex=ad->N;	// NOTE: hack to prevent processing on blocks in LMC_check_col

    if ( dyndata->col==0 ) {
        /* update information of static variables */
        tmpStatic.update ( ad );
    }

    array_t *oarray = clone_array ( original, ad->N, ad->ncols );
    dyndata_t dyndatatmp = dyndata_t ( dyndata );
    if ( log_print ( NORMAL, "" ) ) {
        printf ( "LMCreduceFull:\n" );
        print_array ( oarray, ad->N, ad->ncols );
        dyndata->show();
    }
    copy_array ( original, reduction->array, ad->N, ad->ncols );
    //printf("LMCreduceFull: calling LMCreduce_non_root: col %d\n", dyndatatmp.col);
    lmc_t ret = LMCreduce_non_root ( oarray, ad, &dyndatatmp, reduction, oaextend, tmpStatic );

    destroy_array ( oarray );
    delete ad;
    return ret;
}


/*!
  LMC performs an LMC test or LMC reduction on a given array.

  As input we require the array to be tested to be in root form.

  First we select the necessary column permutations. This is done in steps because we cannot interchange columns
    within a single column group.

  \brief LMC test
  \param original Pointer to the original OA
  \param array Pointer to the array to be checked
  \param ad Pointer to the OA data, these data are static
  \param dyndata Pointer to the dynamic OA data
  \return Value for LMC test
  */
lmc_t LMCreduce ( const array_t* original, const array_t *array, const arraydata_t* ad, const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend &oaextend )
{
//  printf("LMCreduce: ...\n");

    LMC_static_struct_t *tpp=0;
    if ( reduction->staticdata==0 ) {
        tpp = getGlobalStaticOnePointer();
#ifdef OADEBUG
        printfd ( "LMCreduce (debugging): acquired LMC_static_struct_t from single global object\n" );
#endif
    } else {
        tpp = ( reduction->staticdata );
    }

    LMC_static_struct_t &tmpStatic = *tpp;

    //printfd("LMCreduce: reduction->init_state %d\n", reduction->init_state);


    if ( dyndata->col==0 ) {
        /* update information of static variables */
        tmpStatic.update ( ad );

        if ( reduction->init_state==INIT_STATE_INVALID ) {
            // TODO: remove this code
            printf ( "LMCreduce: reduction.init_state is INIT_STATE_INVALID\n" );
            throw;
        }

        if ( reduction->init_state==COPY ) {
            reduction->setArray ( array, ad->N, ad->ncols );
            reduction->init_state=COPY;
        }

        //#ifdef OADEBUG
        bool rootform = check_root_form ( reduction->array, *ad );

        if ( ! rootform ) {
            printfd( "LMCreduce: WARNING: LMC test or LMC reduction for arrays not in root form needs special initialization! reduction->mode %d (OA_TEST %d, OA_REDUCE %d)\n", reduction->mode, OA_TEST, OA_REDUCE );

            int changed = check_root_update ( original, *ad, reduction->array );

            if ( checkloglevel ( DEBUG ) || 0 ) {
                printf ( "original:\n" );
                print_array ( original, ad->N, ad->ncols );
                printf ( "reduction:\n" );
                print_array ( reduction->array, ad->N, ad->ncols );
            }
        }
        //#endif

    }

#ifdef OACHECK
    if ( ad->ncols<=ad->strength ) {
        printf ( "LMC reduction not defined for arrays with number of columns (%d) less or equal to the strength (%d)\n", ad->ncols, ad->strength );
        return LMC_NONSENSE;
    }
#endif


    lmc_t ret = LMC_NONSENSE;

    levelperm_t * lperm_p = 0;
    colperm_t * colperm_p = 0;
    tmpStatic.init_root_stage ( lperm_p, colperm_p, ad );
    log_print ( DEBUG+1, "LMC: entering root LMC phase (level %d) \n", dyndata->col );

    /* check whether the root is full */
    if ( dyndata->col==ad->strength ) {
        if ( log_print ( DEBUG, "" ) ) {
            log_print ( DEBUG, "## LMC: root column selection complete \n" );
            printf ( "   perm: " );
            print_perm ( dyndata->colperm, ad->ncols );

        }

#ifdef OAMEM
        dyndata_t *dd = new dyndata_t ( dyndata );
        dd->col = 0;
        ret = LMCreduce_root_level_perm_ME ( original, ad, dd,reduction );
        delete dd;
#else
        ret = LMCreduce_root_level_perm ( original, ad, dyndata,reduction, oaextend, tmpStatic );
#endif
        //logstream(DEBUG) << printfstring("  return from LMCreduce_root_level_perm: column %d\n", reduction->lastcol);

        return ret;
    }

    int col = dyndata->col;	// current column
    int cg = ad->get_col_group ( col );
    const int ncolsremroot = ad->strength-col;
    const int ncolsremgroup = ad->colgroupsize[cg]+ad->colgroupindex[cg]-col;
    int k = min ( ncolsremroot, ncolsremgroup );	/* number of columns to select in this stage */

    log_print ( DEBUG+1, "LMC: root stage: active column %d, selecting %d columns\n", col, k );
    /* select k columns from the entire group */

    int nc = ncombs ( ncolsremgroup, k );
    colindex_t *comb = new_comb_init<colindex_t> ( k );
    colindex_t *localcolperm = new_perm<colindex_t> ( ad->ncols );

    array_t *cpy = create_array ( ad );

    dyndata_t newdyndata = dyndata_t ( dyndata );
    /* loop over all possible column combinations for the root */

    for ( int i=0; i<nc; i++ ) {
        //if(dyndata->col==0 && (i%10)==0) {
#ifdef OAEXTRA
        logstream ( DEBUG ) << printfstring ( "LMCreduce: root stage (col %d): column comb: %d/%d\n", dyndata->col, i, nc ); // print_perm(comb, k);
#endif
        //}

        create_perm_from_comb<colindex_t> ( localcolperm, comb, k, ad->ncols, col );

        // loop over permutations of selected columns
        for ( int z=0; z<factorial<int> ( k ); z++ ) {

            // set necessary components of dyndata (these might be changed inside LMCreduce)
            newdyndata.copydata ( *dyndata );

            newdyndata.col=newdyndata.col+k;

            perform_inv_perm ( dyndata->colperm, colperm_p[col], ad->ncols, localcolperm );
            copy_perm ( colperm_p[col], newdyndata.colperm, ad->ncols );

            if ( log_print ( DEBUG, "" ) ) {
                printf ( "  colperm root selection " );
                print_perm ( newdyndata.colperm, ad->ncols );
            }

            ret = LMCreduce ( original, array,ad, &newdyndata, reduction, oaextend );

#ifdef OADEBUG2
            printf ( "col %d: colperm ", dyndata->col );
            print_perm ( newdyndata.colperm, ad->ncols ); //printf("\n");
#endif

            if ( ret==LMC_LESS ) {
                log_print ( DEBUG+1, "LMC at level %d: received LMC_less, doing break\n" , dyndata->col );
                break;
            }

            next_perm ( localcolperm+col, k );
        }
        if ( ret==LMC_LESS ) {
            //printf("LMCreduce: found LMC_LESS\n");
            break;
        }
        next_combination ( comb, k, ncolsremgroup );	// increase combination
    }

    destroy_array ( cpy );
    delete_perm ( localcolperm );
    delete_comb ( comb );

    return ret;
}

/**
* @brief Print the contents of a rowsort structure
* @param rowsort Pointer to rowsort structure
* @param N Number of elements
*/
void print_rowsort ( rowsort_t *rowsort, int N )
{
    printf ( "row: " );
    for ( int x=0; x<N; x++ )
        printf ( "%2ld ", static_cast<long> ( rowsort[x].r ) );
    printf ( "\nval: " );
    for ( int x=0; x<N; x++ )
        printf ( "%2ld ", static_cast<long> ( rowsort[x].val ) );
    printf ( "\n" );
}

void print_column_rowsort2 ( const array_t *arraycol, const array_t *arraycol2, rowsort_t *rowsort, levelperm_t lperm, int N )
{
    printf ( "column: \n" );
    for ( int x=0; x<N; x++ )
        printf ( "row %d: %d (orig %d) (value %ld)\n", x, lperm[arraycol[rowsort[x].r]], arraycol2[x], ( long ) rowsort[x].val );
}

void print_column_rowsort ( const array_t *arraycol, rowsort_t *rowsort, int N )
{
    printf ( "column: \n" );
    for ( int x=0; x<N; x++ )
        printf ( "%d: %d (value %ld)\n", x, arraycol[rowsort[x].r], ( long ) rowsort[x].val );
}



/// Print a column of an array with order determined by rowsort structure
void print_col_sorted ( carray_t *array, levelperm_t lperm, rowsort_t *rs, int N )
{
    printf ( "do lperm: " );
    for ( int i=0; i<N; i++ )
        printf ( "%d ", lperm[array[rs[i].r]] );
    printf ( "\n" );
    printf ( "no lperm: " );
    for ( int i=0; i<N; i++ )
        printf ( "%d ", array[rs[i].r] );
    printf ( "\n" );
}

/// Predict J-characteristic of an array based on rowsort of the root
inline int predictJrowsort ( const array_t *array, const int N, const rowsort_t *rs, levelperm_t lperm )
{
    int t = N/4;
    int tt = t/2;

    // x1 is the number of + entries in (s)
    int x1=0;
    for ( int i=0; i<tt; i++ ) {
        int ii = rs[i].r;
        if ( lperm[array[ii]]==0 )
            x1++;
    }
    for ( int i=tt; i<t; i++ ) {
        // TODO: value for second loop can be calculated from result of first loop
        int ii = rs[i].r;
        if ( lperm[array[ii]]==1 )
            x1++;
    }

    return 8*x1-N;
}



/// select the unique arrays in a list, the original list is sorted in place
void selectUniqueArrays ( arraylist_t &xlist, arraylist_t &earrays, int verbose )
{
    sort ( xlist.begin(), xlist.end() ); //solutions.sort();
    if ( xlist.size() >0 ) {
        std::vector<int> vv ( xlist.size() );
        vv[0]=1;

        for ( size_t i=0; i<xlist.size()-1; i++ ) {
            int e = xlist[i]==xlist[i+1];
            if ( !e ) {
                vv[i+1]=1;
            }
            //xlist[i].show();
            if ( verbose>=3 && i<10 ) {
                printf ( "--   array %d==%d: %d\n", ( int ) i, ( int ) i+1, e );
            }
        }

        for ( size_t i=0; i<xlist.size(); i++ ) {
            if ( vv[i] ) {
                if ( ( verbose>=3 && i < 10 ) || verbose>=4 ) {
                    printf ( "  selecting array %d\n", ( int ) i );
                }
                earrays.push_back ( xlist[i] );
            }
        }
    }
}


array_link reduceLMCform ( const array_link &al )
{
    int strength=2; // assume strength is at least 2

    arraydata_t ad = arraylink2arraydata ( al, 0, strength );
    LMCreduction_t *reduction = new LMCreduction_t ( &ad );
    dyndata_t dynd = dyndata_t ( ad.N );
    const array_t *array=al.array;
    reduction->mode = OA_REDUCE;
    reduction->setArray ( al );
    OAextend oaextend;
    lmc_t result = LMCreduce ( array, array, &ad, &dynd, reduction, oaextend );

    array_link nf = reduction->transformation->apply ( al );

    return nf;
}

array_link reduceDOPform ( const array_link &al, int verbose )
{
    int dopruning=0;
    int dolmc=1;
    //assert(strength>=2);
    int strength=2;
    arraylist_t lst;
    lst.push_back ( al );
    arraylist_t earrays;
    if ( verbose>=2 )
        printf ( "reduceDOPform: calling reduceArraysGWLP\n" );
    reduceArraysGWLP ( &lst, earrays,  verbose, dopruning,  strength,  dolmc );
    return earrays[0];
}

lmc_t reduceCanonicalGWLP2 ( const array_link &al, const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t &reduction, array_link &lmc, int verbose )
{

    //printf("reduceCanonicalGWLP: untested\n");

    lmc_t ret = LMC_NONSENSE;
    if ( 0 ) {
        // 1. calculate delete-1-projection values
        std::vector< GWLPvalue > dopgma = projectionGWLPs ( al );
        if ( verbose ) {
            printf ( "LMCcheckGWP: delete-1 GWP:        " );
            display_vector< GWLPvalue > ( dopgma );
            cout << endl;
        }

        // calculate symmetry group
        indexsort is ( dopgma );
        std::vector< GWLPvalue > sgma = is.sorted ( dopgma );
        symmetry_group sg ( sgma, 0 );

        if ( verbose>=1 ) {
            printf ( "LMCcheckGWP inside: indices: " );
            is.show();
            printf ( "\n" );
        }

        // create column sorted sorted array
        array_link alx ( al, is.indices );
        arraydata_t adx ( ad );
        adx.set_colgroups ( sg );

        // do normal reduction
        reduction.mode=OA_REDUCE;

        ret = LMCcheck ( alx, adx, oaextend, reduction );
        //reduction.show();
    } else {
        reduction.mode=OA_REDUCE;
        ret = LMCcheck ( al, ad, oaextend, reduction );

    }
    copy_array ( reduction.array, lmc.array, lmc.n_rows, lmc.n_columns );

    if ( verbose>=2 ) {
        printf ( "LMCcheckGWLP: ret %d\n", ret );
    }
    // fix column sorting
    return ret;
}




/// add mixed values
std::vector< GWLPvalue > mixedProjGWLP ( const std::vector< GWLPvalue > dopgwp, const arraydata_t &ad, int verbose = 0 )
{
    std::vector< GWLPvalue > xx;
    int nc = dopgwp.size();
    for ( int i=0; i<nc; i++ ) {
        GWLPvalue w = dopgwp[i];
        std::vector<double> t;
        t.push_back ( -ad.s[i] );

        t.insert ( t.end(), w.v.begin(), w.v.end() );

        if ( verbose>=2 ) {
            printf ( "reduceArraysGWLP: mixed array: column %d: s[%d]=%d: \n   ", i, i, ad.s[i] );
            display_vector<double> ( t );
            printf ( "\n" );
        }
        xx.push_back ( t );
    }


    if ( verbose>=3 ) {
        printf ( "old indexsort:\n" );
        //is.show();
    }

    return xx;
}

array_transformation_t reductionDOP ( const array_link &al, int verbose )
{
    int strength=al.strength();

    arraydata_t ad =arraylink2arraydata ( al );


    std::vector< GWLPvalue > dopgwp = projectionGWLPs ( al );
    if ( 0 ) { // hack
        printfd ( "hack\n" );
        for ( int i=0; i<5; i++ ) {
            array_link d = al.deleteColumn ( i );
            std::vector<double> gma = GWLP ( d );
            printfd ( " i %d\n", i );
            al.show();
            if ( i<3 ) {
                std::copy ( al.array, al.array+al.n_rows* ( 5-1 ), d.array );
            }
            if ( i==3 && 0 ) {
                std::copy ( al.array, al.array+al.n_rows* ( 5-1 ), d.array );

            }
            gma = GWLP ( d );
            if ( gma[3]==1.1 ) {
                printf ( "huh?\n" );

            }
            for ( int j=0; j<0; j++ )
                gma[4]=i;
            printf ( "gma %d: ", i );
            display_vector ( gma );
            printf ( "\n" );
            dopgwp[i]= gma;

            //dopgwp[i]=i;
        }
    }

    indexsort is ( dopgwp );


    std::vector< GWLPvalue > dofvalues = dopgwp;

    if ( ad.ismixed() ) {
        dofvalues = mixedProjGWLP ( dopgwp,ad, verbose );
    }

    is.init ( dofvalues );


    std::vector<DOFvalue> sdofvalues = is.sorted ( dofvalues );

    if ( verbose>=2 ) {
        printf ( "reductionDOP: delete-1-factor values sorted: " );
        display_vector< DOFvalue > ( sdofvalues );
        std::cout << endl;
    }

    symmetry_group sg ( sdofvalues, 0 );
    if ( verbose>=2 )
        sg.show();

    array_transformation_t at ( &ad );
    at.setcolperm ( is.indices );

    array_link alf ( al, is.indices );

    if ( verbose ) {
        printf ( "is sorting: " );
        display_vector< int > ( is.indices );
        printf ( "\n" );

    }
    // NOTE: since the array was re-ordered, the symmetry group has to be re-ordered
    ad.set_colgroups ( sg );

    OAextend oaextend;
    oaextend.setAlgorithm ( OAextend::getPreferredAlgorithm ( ad ) );
    oaextend.setAlgorithm ( MODE_ORIGINAL );

    LMCreduction_t reduction ( &ad );

    //array_link lm ( al );

    if ( verbose>=3 )
        ad.show ( 2 );


    reduction.mode=OA_REDUCE;
    reduction.init_state=COPY;
    {
        reduction.init_state=INIT;
        //reduction.init_state=COPY;
        reduction.setArray ( alf );
        int changed = check_root_update ( alf.array, ad, reduction.array );
        copy_array ( alf.array, reduction.array, al.n_rows, al.n_columns ); // hack?
    }

    lmc_t ret = LMCcheck ( alf, ad, oaextend, reduction );
    array_transformation_t tmp = ( * ( reduction.transformation ) ) * at ;
    return tmp;
}
/// reduce arrays to canonical form using delete-1-factor ordering
void reduceArraysGWLP ( const arraylist_t *arraylist, arraylist_t &earrays, int verbose, int dopruning, int strength, int dolmc )
{
    OAextend oaextend;
    arraylist_t xlist;
    int npruned=0;

    if ( verbose>=2 )
        printf ( "reduceArraysGWLP: start\n" );

    for ( size_t i=0; i< ( size_t ) arraylist->size(); i++ ) {
        if ( verbose>=2 )
            printf ( "reduceArrays: array %d/%d\n", ( int ) i, ( int ) arraylist->size() );
        const array_link al = arraylist->at ( i );
        // create arraydatya
        int ncols=al.n_columns;

        arraydata_t ad = arraylink2arraydata ( al, 0, strength );

        if ( verbose>=4 ) {
            std::vector<double> gwp = GWLP ( al );
            cout << "GMA: ";
            printf_vector<double> ( gwp, "%.3f " );
            cout << endl;
        }
        std::vector< GWLPvalue > dopgwp = projectionGWLPs ( al );

        if ( verbose>=3 ) {
            printf ( "  delete-1 GWP sequence:        " );
            display_vector< GWLPvalue > ( dopgwp );
            cout << endl;
        }

        if ( dopruning ) {
            //double x = vectormax<double>(dopgwp);
            GWLPvalue x = * ( min_element ( dopgwp.begin(), dopgwp.begin() +ncols-1 ) );
            if ( verbose>=2 ) {
                printf ( "  delete-1 GWP sequence:        " );
                //printf_vector<GWLPvalue> ( dopgwp, "%.3f " );
                display_vector< GWLPvalue > ( dopgwp );

                cout << endl;
                std::cout << "  pruning check: " << dopgwp[ncols-1] << " minmax " << x << std::endl;
            }
            if ( dopgwp[ncols-1]>x ) {
                if ( verbose>=2 ) {
                    printf ( "  reject based on dopgwp ordering\n" );
                }
                npruned++;
                continue;
            }
        }

        if ( verbose>=2 )
            printfd ( "reduceArraysGWLP:\n" );

        std::vector<DOFvalue> dofvalues = dopgwp;
        indexsort is ( dopgwp );

        if ( ad.ismixed() ) {
            //printfd ( "reduceArraysGWLP: warning array is mixed!\n" );

            dofvalues = mixedProjGWLP ( dopgwp,ad, verbose );

            if ( verbose>=3 ) {
                printf ( "old indexsort:\n" );
                is.show();
            }
            //is.init(xx);
        }

        if ( verbose>=2 )
            printfd ( "reduceArraysGWLP:\n" );

        is.init ( dofvalues );
        if ( verbose>=3 ) {
            printf ( "new indexsort:\n" );
            is.show();
            printf ( ", \n" );
        }

        if ( verbose>=3 ) {
            printf ( "  indices: " );
            is.show();
            printf ( "\n" );
        }

        std::vector<DOFvalue> sdofvalues = is.sorted ( dofvalues );

        if ( verbose>=2 ) {
            printf ( "  delete-1-factor values sorted: " );
            display_vector< DOFvalue > ( sdofvalues );

            std::cout << endl;
            //printf("  delete-1 GMA: sorted "); display_vector<double>(sgma);  cout << endl;
        }

        //testTranspose(al);

        symmetry_group sg ( sdofvalues, 0 );
        if ( verbose>=2 )
            sg.show();

        // Done: calculate row symmetry group


        array_link alf ( al, is.indices );

        if ( verbose>=2 )
            printf ( "------\n" );

        // NOTE: since the array was re-ordered, the symmetry group has to be re-ordered
        ad.set_colgroups ( sg );

        OAextend oaextend;
        oaextend.setAlgorithm ( OAextend::getPreferredAlgorithm ( ad ) );
        oaextend.setAlgorithm ( MODE_ORIGINAL );

        LMCreduction_t reduction ( &ad );

        // TODO: adapt extend system

        array_link lm ( al );

        if ( verbose>=3 )
            ad.show ( 2 );


        reduction.mode=OA_REDUCE;
        reduction.init_state=COPY;
        {
            reduction.init_state=INIT;
            reduction.setArray ( alf );

            int changed = check_root_update ( alf.array, ad, reduction.array );
            if ( changed ) {
                //printfd ( "reduceArraysGWLP: array changed %d.\n", changed );
                // reduction.getArray().showarray();
            }
        }

        if ( verbose>=3 ) {
            printfd ( "before LMCcheck\n" );
            oaextend.info();
        }
        lmc_t ret = LMCcheck ( alf, ad, oaextend, reduction );
        copy_array ( reduction.array, lm.array, lm.n_rows, lm.n_columns );
        if ( verbose>=2 )
            printfd ( "reduceArraysGWLP: ret %d (LMC_MORE %d, strength %d)\n", ret, LMC_MORE, ad.strength );

        if ( alf==lm ) {
            if ( verbose>=3 )
                printf ( "  array unchanged (ret %d)\n", ret );
            if ( verbose>=3 ) {
                printf ( "input sort:\n" );
                alf.showarray();
                printf ( "output:\n" );
                lm.showarray();
            }
        } else {
            if ( verbose>=3 ) {
                printf ( "input sort:\n" );
                alf.showarray();
                printf ( "output:\n" );
                lm.showarray();
            }
            if ( verbose>=2 )
                printf ( "  array changed\n" );
        }
        xlist.push_back ( lm );
    }
    if ( verbose ) {
        printf ( "reduceArraysGWLP: input size %d, output size %d, pruned %d\n", ( int ) arraylist->size(), ( int ) xlist.size(), npruned );
    }


    if ( dolmc ) {
        selectUniqueArrays ( xlist, earrays, verbose );

    } else {
        sort ( xlist.begin(), xlist.end() ); //solutions.sort();

        for ( size_t i=0; i<xlist.size(); i++ ) {
            earrays.push_back ( xlist[i] );
        }

    }
    if ( verbose )
        printf ( "  selecting %d/%d arrays\n", ( int ) earrays.size(), ( int ) xlist.size() );

}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
