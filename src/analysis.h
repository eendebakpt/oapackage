/** \file analysis.h

 \brief Contains code for analysis of the OA algorithm.

 The code can compiled or left out of the program using the OAANALYZE constant.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>
 Copyright: See LICENSE.txt file that comes with this distribution
*/

#ifdef OAANALYZE

/* static struct with all the performance counters */
struct analysisdata {
    static double Tlmc;
    static double Tnonlmc;
    static double Tstatic;
    static double Tgenerate;
    static double Ttotal;

    static long nextchecked ;
    static long nLMC ;
    static long nNonLMC ;

};

inline void resetAnalysisData()
{
    analysisdata::Tlmc = 0;
    analysisdata::Tnonlmc = 0;
    analysisdata::Tstatic = 0;
    analysisdata::Tgenerate = 0;
    analysisdata::Ttotal = 0;

    analysisdata::nextchecked  = 0;
    analysisdata::nLMC  = 0;
    analysisdata::nNonLMC  = 0;
}

inline double addTtotal(double dt = 0)
{
    analysisdata::Ttotal += dt;
    return analysisdata::Ttotal ;
}

inline double addTgenerate(double dt = 0)
{
    analysisdata::Tgenerate += dt;
    return analysisdata::Tgenerate ;
}

inline double addTstatic(double dt = 0)
{
    //cout << "adding to T static: " << setw(4) << dt << endl;
    analysisdata::Tstatic += dt;
    return analysisdata::Tstatic ;
}

inline double addTlmc(double dt = 0)
{
    analysisdata::Tlmc += dt;
    return analysisdata::Tlmc ;
}

inline double addTnonlmc(double dt = 0)
{
    analysisdata::Tnonlmc += dt;
    return analysisdata::Tnonlmc ;
}


#endif
