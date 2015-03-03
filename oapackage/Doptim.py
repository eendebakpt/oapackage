# -*- coding: utf-8 -*-
"""

Collection of helper functions for OA package

Pieter Eendebak <pieter.eendebak@gmail.com>

@author: eendebakpt
"""

from __future__ import print_function

import oalib
import sys, os
import numpy as np
import time
try:
    import matplotlib
    import matplotlib.pyplot as plt
except:
    print('Doptim: matplotlib cannot be found, not all functionality is available')
    pass
import oapackage.markup as markup
import oapackage.oahelper as oahelper

from oapackage.markup import oneliner as e

#%%

def array2Dtable(sols, verbose=1, titlestr=None):
    #page.p()
    na=len(sols)
    page=markup.page()
    page.table( style=' border-collapse: collapse;')
    page.tr(style='font-weight: bold; border-bottom: solid 1px black;')
    page.th( 'Array', style='padding-right:30px; ' )
    page.th( ('D-efficiency', 'Ds-efficiency', 'D1-efficiency'), style='padding-right:14px;'  )
    page.th( ('GWLP'), style='padding-right:14px;'  )
    page.tr.close()
    for ii, al in enumerate(sols):
        aidx=ii
        (D, Ds, D1)=al.Defficiencies()
        gwlp=al.GWLP()
        page.tr(style='font-weight: normal;')
        page.td( '%d' % aidx, style='padding-right:10px;' )
        for v in [D,Ds,D1]:
            page.td( '%.4f' % v, style='padding-right:1px;' )
        gstr=oahelper.gwlp2str(gwlp)
        page.td( e.small(gstr) , style='padding-right:1px;' )
        page.tr.close()
    page.table.close()
    #page.p.close()
    return page

#%%

try:
	import brewer2mpl
except:
	pass


def generateDscatter(dds, si=0, fi=1, lbls=None, nofig=False, fig=20):        
    data=dds.T
    pp = oahelper.createPareto(dds)
    paretoidx=np.array(pp.allindices())

    nn=dds.shape[0]
    area=40*np.ones( nn,) 
    area[np.array(pp.allindices())]=80
    alpha=1.0
    
    if dds.shape[1]>3:
        colors=dds[:,3]
    else:
        colors=np.zeros( (nn,1) )
    
    try:
        mycmap = brewer2mpl.get_map('Set1', 'qualitative', 3).mpl_colors
    except:
        mycmap= [ matplotlib.cm.jet(ii) for ii in range(4)]
        pass
            
    
    idx=np.unique(colors).astype(int)
    

    
    # For remaining spines, thin out their line and change the black to a slightly off-black dark grey
    almost_black = '#202020'
    
    figh=plt.figure(fig) # ,facecolor='red')
    plt.clf()
    figh.set_facecolor('w')
    ax = plt.subplot(111)
    
    nonparetoidx=np.setdiff1d(range(nn), paretoidx)
    ax.scatter(data[fi,nonparetoidx], data[si,nonparetoidx], s=33, c=(.5,.5,.5), linewidths=0, alpha=alpha, label='Non-pareto design')

    for jj, ii in enumerate(idx):
        gidx=(colors==ii).nonzero()[0]
        gp=np.intersect1d(paretoidx, gidx)
        
        color = mycmap[jj]
        cc=[color]*len(gp)
        print('index %d: %d points' % (ii, gidx.size))
        ax.scatter(data[fi,gp], data[si,gp], s=52, c=cc, linewidths=0, alpha=alpha, label=lbls[jj])
        plt.draw()
    

    if data[si,:].std()<1e-3:        
        y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)

    if 0:
        for xi,al in enumerate(sols):
            D, Ds, D1 = al.Defficiencies()
            print('D1 %f Ds %f D %f' % (D1, Ds,D))    
        
            tmp=plt.scatter(Ds, D, s=60, color='r');
            if xi==0:
                tmp=plt.scatter(Ds, D, s=60, color='r', label='Strength 3');
        plt.draw()
        
    plt.xlabel('$D_s$-efficiency', fontsize=16)
    plt.ylabel('D-efficiency', fontsize=16)

    plt.axis('image')    
    plt.legend=ax.legend(loc=3, scatterpoints=1) #, fontcolor=almost_black)
    if not nofig:
        plt.show()
    #time.sleep(0.01)
    plt.grid(b=True, which='both', color='0.85',linestyle='-')
    
    
    if not nofig:
        plt.draw()


#%%

#import researchOA

def generateDpage(outputdir, arrayclass, dds, allarrays, fig=20, optimfunc=[1,0,0],nofig=False,  urlprefix='', makeheader=True, verbose=1, lbls=None):
    
    #%% Prepare data
    pp = oahelper.createPareto(dds)
    paretoidx=np.array(pp.allindices())

    narrays=dds.shape[0]
    npareto=pp.number()
    

    xstr=oahelper.series2htmlstr(arrayclass, case=1)
    xstrplain=oahelper.series2htmlstr(arrayclass, html=0, case=1)

    paretoarrays = oahelper.selectParetoArrays(allarrays, pp)
    at= array2Dtable(paretoarrays, verbose=1)
 
    pfile0='paretoarrays.oa'
    pfile=os.path.join(outputdir, pfile0)
    oalib.writearrayfile(pfile,  paretoarrays) 
            
    istrlnk = markup.oneliner.a('paretoarrays.oa', href=urlprefix+pfile0)

    #lbls= ['Optimization of $D+.5 D_s$', 'Optimization of $D+  0.5 D_s$', 'Optimization of $D+3*Ds$', 'Optimization of $D+3*D_s$']
    if lbls is None:
        lbls= ['Optimization of $D$'] 
    x=generateDscatter(dds, lbls=lbls, fig=fig, nofig=nofig)
 
    scatterfile=os.path.join(outputdir, 'scatterplot.png')
    if verbose:
        print('generateDpage: writen scatterplot to %s' % scatterfile)
    plt.savefig(scatterfile, bbox_inches='tight', pad_inches=0.25)

    #%% Create page

    page = markup.page( )
        
    if makeheader:        
        page.init( title="Class %s" % xstrplain, 
                       css=( '../oastyle.css' ), 
                       lang='en', htmlattrs=dict({'xmlns': 'http://www.w3.org/1999/xhtml', 'xml:lang': 'en'}),
                       header="<!-- Start of page -->", 
                       bodyattrs =dict({'style': 'padding-left: 3px;'}),
                       #doctype=markup.doctype.strict,
                       doctype='<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
                       metainfo=( {'text/html': 'charset=utf-8', 'keywords': 'orthogonal arrays designs', 'robots': 'index, follow', 'description': 'Even-Odd arrays'}),
                       footer="<!-- End of page -->" )
                                      
        
    page.h1('Results for array class %s ' %  xstr)
    

    # mathjax is not updated properly...    
    ss='The Pareto optimaly was calculated according to the statistics \(D\), \(D1\) and \(Ds\).'
    ss='The Pareto optimaly was calculated according to the statistics D, D<sub>1</sub> and D<sub>s</sub>.'
    page.p('Generated %d arrays, %d are Pareto optimal. %s' % (narrays, npareto, ss) )

    if narrays>0:

        scores=calcScore(dds, optimfunc)
        tmp, dd, sols = selectDn(scores, dds, allarrays, nout=1)
        A=sols[0]
        bestdesignfile=os.path.join(outputdir, 'best-design.oa')
        oalib.writearrayfile(bestdesignfile, A)
        
        page.h2('Best design')
        page.p('The best design: %s.' % e.a('best-design.oa', href=os.path.join(urlprefix, 'best-design.oa')))
        
        page.p()
        
        dd=dd[0]
        page.span('D-efficiency: %.4f, ' % A.Defficiency() ); #page.br()
        page.span('D<sub>s</sub>-efficiency: %.4f, ' % dd[1] );# page.br()
        page.span('D<sub>1</sub>-efficiency: %.4f' % dd[2] ); page.br()
        page.span('A-efficiency: %.3f' % A.Aefficiency() ); page.br()
        gwlp=A.GWLP()
        #gwlp=','.join(['%.3f' % xq for xq in gwlp])
        gwlp=oahelper.gwlp2str(gwlp, jstr=', ')
        page.span('Generalized wordlength pattern: %s' % gwlp ); page.br()
        #page.p('D-efficiency: %.3f' % A.Defficiency() )
        pec=oalib.PECsequence(A)
        pec=','.join(['%.3f' % xq for xq in pec])
        page.span('PEC-sequence: %s' % pec )
        page.br()
    
    page.h2('Table of Pareto optimal arrays ')
    
    
    page.span(str(at) )
    page.p('All Pareto optimal arrays: %s' % istrlnk )
    
    
    page.img(src=urlprefix+'scatterplot.png', style="margin: 10px;")
    
    citationstr = markup.oneliner.a('Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays', href='http://dx.doi.org/10.1002/jcd.20236')
    
    
    page.br(clear='both')
    page.p('Citation notice: if you make use of the results on this page, please cite the following paper:')
    page.p('%s, Journal of Combinatorial Designs, Volume 18, Issue 2, pages 123-140, 2010.' % citationstr)
    
    outfile=os.path.join(outputdir ,'Dresults.html' )
    fid=open(outfile, 'wt')
    fid.write(str(page))
    fid.close()
    print('written to file %s' % outfile)

    return outfile

#generateDpage(doptimdir, arrayclass, dds, sols, optimfunc=[1,0,0], makeheader=True)

#%%

#%%
# TODO: do the arrays need to have strength 1?
# TODO: implement in C
# TODO: what is the proper criterium?
# TODO: implement fast update of xf
# TODO: add check to capture mixed arrays

def optimDeffhelper(classdata):
    """ Helper function that is suitable for the multi-processing framework """

    N=classdata[0]
    k=classdata[1]
    alpha=classdata[2]
    method=classdata[3]
    p=classdata[4]
    nabort=p.get('nabort', 2500)    
    niter=p.get('nabort', 12000)
    
    arrayclass=oalib.arraydata_t(2, N, 1, k)
    al=arrayclass.randomarray(1)

    vv= optimDeffPython(al, niter=niter, nabort=nabort, verbose=0, alpha=alpha, method=method)
    return vv[0], vv[1].getarray()
    
def calcScore(dds, optimfunc):
    n=dds.shape[0]
    scores=np.zeros( n, )    
    for ii,d in enumerate(dds):
        alpha=optimfunc        
        scores[ii]= alpha[0]*d[0]+alpha[1]*d[1]+alpha[2]*d[2]
    return scores
        
def optimDeffPython(A0, arrayclass=None, niter=10000, nabort=2500, verbose=1, alpha=[1,0,0], method=0):
    """ Optimize arrays """
    # get factor levels
    if arrayclass is None:
        s=A0.getarray().max(axis=0)+1
    else:
        s=arrayclass.getS()
        arrayclass
    sx=tuple(s.astype(np.int64))
    sg=oalib.symmetry_group( sx )    
    gidx=tuple(sg.gidx)
    gsize=tuple(sg.gsize)
    gstart=tuple(sg.gstart)
    
    nx=0
    
    Dinitial=A0.Defficiency()
    if verbose:
        print('optimDeff: initial Deff %.4f' % Dinitial )
    N=A0.n_rows
    k=A0.n_columns

    # initialize score
    if str(type(alpha))=="<type 'function'>":
        d = alpha(A0)
    else:
        D, Ds, D1=A0.Defficiencies()
        d = alpha[0]*D+alpha[1]*Ds+alpha[2]*D1
    A=A0.clone()
    lc=0
    for ii in range(0,niter): 
        #r=random.randint(0,N-1)
        #c=random.randint(0,k-1)
        r=np.random.randint(N)
        c=np.random.randint(k)
        r2=np.random.randint(N)
        # make sure c2 is in proper column
        #c2=np.random.randint(k)
        c2 = gstart[gidx[c]] + oalib.fastrand() % gsize[gidx[c]]
        #print('  %d,%d <--> %d,%d' % (r,c,r2,c2))


        #o=A[r, c]; o2=A[r2, c2]
        o=A._at(r,c); o2=A._at(r2,c2) # no extra error checking
        # swap
        if method==oalib.DOPTIM_SWAP:
            # swap
            if o==o2:
                continue
            A._setvalue(r,c,o2); A._setvalue(r2,c2,o)   
        elif method==oalib.DOPTIM_NONE:
            tmp=0
        else:
            # flip
            A._setvalue(r,c,1-o); 
            
        if str(type(alpha))=="<type 'function'>":
            dn = alpha(A)
        else:
            D, Ds, D1=A.Defficiencies()
            nx=nx+1
            #print(alpha)
            dn = alpha[0]*D+alpha[1]*Ds+alpha[2]*D1
        if (dn>=d):
            if dn>d:
                lc=ii
                d=dn
            if verbose>=2:
                print('ii %d: %.6f -> %.6f' % (ii, d, dn))
                
        else:
            # restore to original
            if method==oalib.DOPTIM_SWAP:
                A._setvalue(r,c,o)
                A._setvalue(r2,c2,o2)
            elif method==oalib.DOPTIM_NONE:
                tmp=0
            else:
                A._setvalue(r,c,o)
        if (ii-lc)>nabort:
            if verbose:
                print('optimDeff: early abort ii %d, lc %d' % (ii, lc))
            break
    
    Dfinal=A.Defficiency()
    

    if verbose:
        if Dfinal>Dinitial:
            print('optimDeff: final Deff improved: %.4f -> %.4f' % (Dinitial, A.Defficiency()) )
        else:
            print('optimDeff: final Deff %.4f' % A.Defficiency() )

    if verbose>=2:
        print('nx %d' % nx)            
    return d, A



#%%

def filterPareto(scores, dds, sols, verbose=0):
    pp = oahelper.createPareto(dds, verbose=0)   
    paretoidx=np.array(pp.allindices())

    pscores=scores[paretoidx]   
    pdds=dds[paretoidx]       
    psols=[sols[i] for i in paretoidx]

    return pscores, pdds, psols

#%%
   
def selectDn(scores, dds, sols, nout=1):
    idx = np.argsort(-scores.flatten() )        
    #print(idx.shape)
    #print(scores.shape)
    scores=scores[idx]
    dds=dds[idx,:]
    sols=[sols[ii] for ii in idx]
    if not nout is None:
        # sort the arrays
        scores=scores[0:nout]
        dds=dds[range(nout),:]
        sols=sols[0:nout]
    return scores, dds, sols
    
def Doptimize(arrayclass, nrestarts=10, niter=12000, optimfunc=[1,0,0], verbose=1, maxtime=180, selectpareto=True, nout=None, method=0):
    """ Calculate D-optimal designs
    
    
    """
    if verbose:
        print('Doptim: optimization class %s' % arrayclass.idstr() )
    t0=time.time()

    if optimfunc is None:
        optimfunc=[1,2,0]
        
    scores=np.zeros( (0,1))
    dds=np.zeros( (0,3)) 
    sols=[]  #oalib.arraylist_t()
       
    for ii in range(nrestarts):
        if verbose:
            oahelper.tprint('Doptim: iteration %d/%d (time %.1f/%.1f)' % (ii, nrestarts, time.time()-t0, maxtime), dt=4)
        al=arrayclass.randomarray(1)


        
        #array_link  optimDeff(const array_link &A0,  arraydata_t &arrayclass, std::vector<double> alpha, int verbose=1, int method=-1, int niter=10000, int nabort=2500);

        if isinstance(optimfunc, list):
            alpha=optimfunc
            Ax= oalib.optimDeff(al, arrayclass, alpha, verbose>=2, method, niter)
            dd=Ax.Defficiencies()
            score = oalib.scoreD(dd, optimfunc)
        else:
            score, Ax= optimDeffPython(al, verbose=0, niter=niter, alpha=optimfunc, method=method)
            dd=Ax.Defficiencies()
        
        
        if time.time()-t0 > maxtime:
            if verbose:
                print('Doptim: max time exceeded, aborting')
            break

        scores=np.vstack( (scores, [score]) )
        dds=np.vstack( (dds, dd) )
        sols.append(Ax)
        if verbose>=2:
            print('  generated array: %f %f %f' % (dd[0], dd[1], dd[2]))
            
        if selectpareto and ii%502==0:
            scores,dds,sols=filterPareto(scores, dds, sols)

        #print(dds.shape)
        
    if verbose:
        print('Doptim: done' )

    if selectpareto:
        scores,dds,sols=filterPareto(scores, dds, sols)

    #print(dds.shape)
    
    # sort & select
    scores, dds, sols = selectDn(scores, dds, sols, nout=nout)
    
    return scores, dds, sols


#%%
if 0:
    lbls= ['Optimization of $D+.5 D_s$', 'Optimization of $D+  0.5 D_s$', 'Optimization of $D+3*Ds$', 'Optimization of $D+3*D_s$']
    
    x=generateDscatter(dds, si=0, fi=0, lbls=lbls, fig=20)
    
    import pmatlab
    pmatlab.setWindowRectangle(1500,10,w=800,h=600)
    
