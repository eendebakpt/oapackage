/*! \file Deff.h
 *  \brief Contains functions to optimize designs
 * 
 */

#include "arraytools.h"
#include "arrayproperties.h"


#include "Deff.h"



double scoreD(const std::vector<double> dd, const std::vector<double> alpha) 
{
  double v=0;
  for(size_t i=0; i<dd.size();  i++)
    v+= dd[i]*alpha[i];
  return v;
  
}

array_link  optimDeff(const array_link &A0,  arraydata_t &arrayclass,  std::vector<double> alpha, int verbose, int method, int niter, int nabort)
{
  
  int N = arrayclass.N;
  int k = arrayclass.ncols;
    //#get factor levels
        std::vector<int> s = arrayclass.getS();

	if(method==1) {
	if ( arrayclass.is2level() )
	  method=11; 
	}
	array_link A = A0;
    symmetry_group sg=symmetry_group( s );    

 std::vector<int> gidx = sg.gidx;
 
 std::vector<double> dd0 = A0.Defficiencies();

     if (verbose) {
        printf("optimDeff: initial D-efficiency %.4f\n",  dd0[0] );
     }
     
     // initialize score
     double d = scoreD(dd0, alpha);
 //   gsize=tuple(sg.gsize)
 //   gstart=tuple(sg.gstart)

     int lc=0;	// index of last change to array
  
  
//#pragma omp for 
     for(int ii=0; ii<niter; ii++) {
       
       // select random row and column
       int r = fastrandK(N);
       int c = fastrandK(k);

       int r2 = fastrandK(N);
       //int c2 = fastrandK(k);
	// make sure column is proper column group
               int c2 = sg.gstart[sg.gidx[c]] + fastrandK(sg.gsize[gidx[c]]); 

               // get values
	       array_t o = A._at(r,c);  array_t o2 = A._at(r2,c2); // no extra error checking
	       
	       // update
	       switch (method) {
		 case 0: // swap
		  A._setvalue(r,c,o2); A._setvalue(r2,c2,o);
		  break;
		 case 1: // random update
		 {
		   int val = fastrandK( s[c]);
		    A._setvalue(r,c,val); 
		   break;
		 }
		 case 11: // flip
		    A._setvalue(r,c,1-o); 
		   break;
	       }
	       // evaluate
	        std::vector<double> dd = A.Defficiencies();
	            double dn = scoreD(dd, alpha);

	       
	       // switch back if necessary
		    
		            if (dn>=d) {
            if (dn>d)  {
                lc=ii;
            if (verbose>=3) 
                printf("optimDeff: ii %d: %.6f -> %.6f\n", ii, d, dn);
            d=dn;
	    }
	  }
        else {
            // restore to original
         switch(method) {
	   case 0:
	        A._setvalue(r,c,o);
	        A._setvalue(r2,c2,o2);
		break;
	   case 11:
	   case 1:
	        A._setvalue(r,c,o);
		break;
	 }
	}
	
	// check progress
	        if ((ii-lc)>nabort ){
            if (verbose)
                printf("optimDeff: early abort ii %d, lc %d\n", ii, lc);
	    //abort=true;
            break;
		}

     }
     
      std::vector<double> dd = A.Defficiencies();
      double dn = scoreD(dd, alpha);

      return A;
//      return std::pair<array_link, std::vector<double> >(A, dd);
}
/*
def optimDeff(A0, arrayclass=None, niter=10000, nabort=2500, verbose=1, alpha=[1,0,0], method=0):
    """ Optimize arrays """
        s=A0.getarray().max(axis=0)+1
    sx=tuple(s.astype(np.int64))
    
    for ii in range(0,niter): 
....


            

    Dfinal=A.Defficiency()
    

    if verbose:
        if Dfinal>Dinitial:
            print('optimDeff: final Deff improved: %.4f -> %.4f' % (Dinitial, A.Defficiency()) )
        else:
            print('optimDeff: final Deff %.4f' % A.Defficiency() )
            
    return d, A

    */


