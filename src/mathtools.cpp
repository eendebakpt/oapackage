#include "mathtools.h"

#include <numeric>

/*
dummy::dummy() { val = 0; }
dummy::dummy(double d) {
    val=d;
}
dummy::~dummy() {}

int dummy::operator< ( const dummy &rhs ) const
    {
      printf("dummy class: operator <: value %f, %f\n", val, rhs.val);

      return val < rhs.val;
    }
*/

template<class Type>
void symmetry_group::init(const std::vector<Type> vals, bool ascendingx, int verbose)
{
    //verbose=2;
  
    if (verbose>=2) {
        printf("symmetry_group::init: %zu elements: ", vals.size() );
        for(size_t i=0; i<vals.size(); i++)
            std::cout << vals[i] << ", ";
        printf("\n");
    }

    n=vals.size();
    ascending= ascendingx;

    // check we are sorted
    if (verbose>=2)
      printf("symmetry_group::init: check sorting\n");
    
    indexsort is(vals.size());

    if (ascending)
        is.sort(vals);
    else
        is.sortdescending(vals);

    //      printf("symmetry_group::init: xxx\n");

    if (verbose>=2 || 0) {
        if (ascending) {
            if (! is.issorted() ) {
                printf("symmetry_group: input group was not sorted!\n");
                is.show();
                printf("\n");
            }
        } else {
            if (! is.issorted() ) {
                printf("symmetry_group: input group was not sorted!\n");
                is.show();
                printf("\n");
            }

        }

    }
    // calc group
    int nsg=0;
    Type prev; //=vals[0];

       //       printf("symmetry_group::init: yy\n");

    prev = std::numeric_limits<Type>::quiet_NaN();
    /* count number of symmetry groups */
    for(int i=0; i<n; i++) {
        if(vals[i]!=prev || i==0) {
            nsg++;
            if (verbose) {
	      printf("  symmetry_group: %d: add ", i);
	      std::cout << prev << "->";
	      std::cout << vals[i] << std::endl;
	    }
        }
        prev=vals[i];
    }

    if (verbose) {
        printf("symmetry_group: %d elements, %d groups\n", n, nsg);

    }

    this->ngroups=nsg;

// initialize structures
    gidx.resize(n);
    gstart.resize(nsg+1);
    gsize.resize(nsg+1);


    /* find starting positions of symmetry group */
    if (verbose>=2)
      printf("symmetry_group::init:  find starting positions\n");
    
    prev =  std::numeric_limits<Type>::quiet_NaN();
    nsg = 0;
    int i;
    for(i=0; i<n; i++) {
        //printf("symm_group_index: nsg: %d i %d\n", nsg, i);
        if(vals[i]!=prev || i==0) {
            gstart[nsg]=i;
            nsg++;
        }
        gidx[i]=nsg-1;
        prev=vals[i];
    }
    gstart[nsg]=n;	/* add dummy element */

    for(int i=0; i<nsg; i++) {
        gsize[i] = gstart[i+1]-gstart[i];
    }
    gsize[nsg]=0;
}

symmetry_group::symmetry_group(const std::vector<float> vals, bool ascending, int verbose)
{
    this->init<float>(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector<double> vals, bool ascending, int verbose)
{
    this->init<double>(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector<short int> vals, bool ascending, int verbose)
{
    this->init(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector<unsigned int> vals, bool ascending, int verbose)
{
    this->init(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector<int> vals, bool ascending, int verbose)
{
    this->init<int>(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector< mvalue_t<double> > vals, bool ascending, int verbose)
{
    this->init(vals, ascending, verbose);
}

symmetry_group::symmetry_group(const std::vector<mvalue_t<int> > vals, bool ascending, int verbose)
{
  //printf("symmetry_group::symmetry_group: type <mvalue_t<int>: %zu, %zu\n", vals.size(), vals[0].size() );
    this->init(vals, ascending, verbose);
}

    symmetry_group::symmetry_group ( const symmetry_group &sgx ) {
      gidx = sgx.gidx;
      gstart = sgx.gstart;
      gsize = sgx.gsize;
      ngroups = sgx.ngroups;
      n = sgx.n;
      ascending = sgx.ascending;
    }
    symmetry_group::symmetry_group (  ) 
    {
      ngroups=0; n=0; ascending=0;
    }
    
    
      void symmetry_group_walker::show(int verbose) const {
   printf("symmetry_group_walker: ");
      if (verbose>=2)
    printf("\n");
   for(size_t i=0; i< (size_t)perms.size(); i++) {
     if (verbose>=2) {
     printf("  block %zu: ", i);
     print_perm(perms[i] );
     }
     else {
     print_perm(perms[i], 100, false ); printf(" ");
     }
   }
   if (verbose==1)
    printf("\n");
  }
  
    bool symmetry_group_walker::nextsub(int g) {
//printf("symmetry_group_walker::nextsub: %d\n", g);
//print_perm(perms[g]);
bool of = next_perm(perms[g]);
//print_perm(perms[g]);
    if (of && g>0)
      of = nextsub(g-1);

    return of;
      
    }

    std::vector<int> symmetry_group_walker::fullperm() const
    {
     std::vector<int> w(sg.n);
     for(int i=0; i<sg.n; i++) w[i]=i;

          std::vector<int> ww(sg.n);
for(size_t j=0; j<perms.size(); j++) {
for(size_t k=0; k<perms[j].size(); k++) {
int offset=sg.gstart[j];
  ww[offset+k] = w[offset+perms[j][k]];
} 
}
     return ww;
    }

// instantiate classes
//template symmetry_group::symmetry_group<int>(std::vector<int>, int);
//template symmetry_group::symmetry_group<double>(std::vector<double>, int);

/* Random number generators */

int g_seed = 123;
void seedfastrand(int s)
{
 g_seed=s; 
}
 int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 
 int fastrandK(int K) {
   return fastrand() % K;
 }

