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
import matplotlib.pyplot as plt

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

import matplotlib

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
    
    mycmap = brewer2mpl.get_map('Set1', 'qualitative', 3).mpl_colors
    
    
    idx=np.unique(colors).astype(int)
    

    
    # For remaining spines, thin out their line and change the black to a slightly off-black dark grey
    almost_black = '#202020'
    
    figh=plt.figure(fig) # ,facecolor='red')
    plt.clf()
    figh.set_facecolor('w')
    ax = plt.subplot(111)
    
    for jj, ii in enumerate(idx):
        gidx=(colors==ii).nonzero()[0]
        gp=np.intersect1d(paretoidx, gidx)
        
        color = mycmap[jj]
        cc=[color]*len(gp)
        print('index %d: %d points' % (ii, gidx.size))
        ax.scatter(data[fi,gp], data[si,gp], s=52, c=cc, linewidths=0, alpha=alpha, label=lbls[jj])
        plt.draw()
    
    nonparetoidx=np.setdiff1d(range(nn), paretoidx)
    ax.scatter(data[fi,nonparetoidx], data[si,nonparetoidx], s=33, c=(.5,.5,.5), linewidths=0, alpha=alpha, label='Non-pareto design')

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
    legend=ax.legend(loc=3, scatterpoints=1) #, fontcolor=almost_black)
    if not nofig:
        plt.show()
    #time.sleep(0.01)
    plt.grid(b=True, which='both', color='0.85',linestyle='-')
    
    
    if not nofig:
        plt.draw()


#%%

#import researchOA

def generateDpage(outputdir, arrayclass, dds, allarrays, optimfunc=[1,0,0],nofig=False,  urlprefix='', makeheader=True, verbose=1):
    
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
    lbls= ['Optimization of $D$'] 
    x=generateDscatter(dds, lbls=lbls, fig=20, nofig=nofig)
 
    scatterfile=os.path.join(outputdir, 'scatterplot.png')
    if verbose:
        print('generateDpage: writen scatterplot to %s' % scatterfile)
    plt.savefig(scatterfile)

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
    
    
    ss='The Pareto optimaly was calculated according to the statistics \(D\), \(D_1\) and \(D_s\)'
    page.p('Generated %d arrays, %d are Pareto optimal. %s' % (narrays, npareto, ss))
    
    page.h2('Table of Pareto optimal arrays ')
    
    
    page.span(str(at) )
    page.p('All Pareto optimal arrays: %s' % istrlnk )
    
    
    page.img(src=urlprefix+'scatterplot.png')
    
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

def optimDeff(A0, niter=2000, verbose=1, alpha=[1,0,0]):
    Dinitial=A0.Defficiency()
    if verbose:
        print('optimDeff: initial Deff %.4f' % Dinitial )
    N=A0.n_rows
    k=A0.n_columns
    d=A0.Defficiency()
    d=0
    A=A0.clone()
    lc=0
    for ii in range(0,niter): 
        #r=random.randint(0,N-1)
        #c=random.randint(0,k-1)
        r=np.random.randint(N)
        c=np.random.randint(k)
        r2=np.random.randint(N)
        c2=np.random.randint(k)


        #o=A[r, c]; o2=A[r2, c2]
        o=A._at(r,c); o2=A._at(r2,c2) # no-extra error checking
        if o==o2:
            continue
        # swap
        A._setvalue(r,c,o2); A._setvalue(r2,c2,o)   
        if alpha is None:
            dn=A.Defficiency()
        else:
            D, Ds, D1=A.Defficiencies()
            dn = alpha[0]*D+alpha[1]*Ds+alpha[2]*D1
        if (dn>=d):
            if dn>d:
                lc=ii
            d=dn
            if verbose>=2:
                print('ii %d: %.4f' % (ii, dn))
                
        else:
            # restore to original
            A._setvalue(r,c,o)
            A._setvalue(r2,c2,o2)
        if (ii-lc)>1200:
            if verbose:
                print('optimDeff: early abort ii %d, lc %d' % (ii, lc))
            break
    
    Dfinal=A.Defficiency()
    

    if verbose:
        if Dfinal>Dinitial:
            print('optimDeff: final Deff improved: %.4f -> %.4f' % (Dinitial, A.Defficiency()) )
        else:
            print('optimDeff: final Deff %.4f' % A.Defficiency() )
            
    return d, A



#%%

def Doptimize(arrayclass, niter=10, optimfunc=[1,0,0], verbose=1, maxtime=20):
    """ ... """
    if verbose:
        print('Doptim: optimization class %s' % arrayclass.idstr() )
    t0=time.time()

    scores=np.zeros( (0,1))
    dds=np.zeros( (0,3)) 
    sols=oalib.arraylist_t()
    for ii in range(niter):
        if verbose:
            oahelper.tprint('Doptim: iteration %d/%d' % (ii, niter))
        al=arrayclass.randomarray(1)

        score, Ax= optimDeff(al, verbose=0, niter=7000, alpha=optimfunc)
        dd=Ax.Defficiencies()
        if time.time()-t0 > maxtime:
            if verbose:
                print('Doptim: max time exceeded, aborting')
            break

        scores=np.vstack( (scores, [score]) )
        dds=np.vstack( (dds, dd) )
        sols.push_back(Ax)
        if verbose:
            print('  generated array: %f %f %f' % (dd[0], dd[1], dd[2]))
    if verbose:
        print('Doptim: done' )
        
    return scores, dds, sols


#%%
if 0:
    lbls= ['Optimization of $D+.5 D_s$', 'Optimization of $D+  0.5 D_s$', 'Optimization of $D+3*Ds$', 'Optimization of $D+3*D_s$']
    
    x=generateDscatter(dds, si=0, fi=0, lbls=lbls, fig=20)
    
    import pmatlab
    pmatlab.setWindowRectangle(1500,10,w=800,h=600)
    
