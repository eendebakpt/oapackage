/** \file mpitools.h

 \brief Legacy functions to work with MPI.

 Author: Pieter Eendebak <pieter.eendebak@gmail.com>, (C) 2008

 Copyright: See COPYING file that comes with this distribution
*/


#ifndef MPITOOLS_H
#define MPITOOLS_H

#include "mpi.h"

#include "extend.h"

#include "tools.h"
#include "arraytools.h"

#ifdef OAEXTEND_SINGLECORE
/* single core implementation */

/** number of master processor */
const int MASTER = 0;

#else
#define OAEXTEND_MPI
/* multi-core implementation using mpi */

/** commands: new array to extend, stop processing, wait for new data */
enum MASTER_COMMANDS {MASTER_EXTENDJOB, MASTER_STOP, MASTER_WAIT, MASTER_LMCJOB, MASTER_ALGORITHM};

/** commands: slave is ready, slave has solutions to be picked up */
enum SLAVE_COMMANDS {SLAVE_READY, SLAVE_SOLUTIONS};

/** number of master processor */
const int MASTER = 0;

/** maximum number of solutions to extend before complete synchronization */
const int nsolsync = 1000;

int stop_slaves();

int collect_solutions_single(arraylist_t &extensions, const arraydata_t *ad);

#endif

#define DEBUGMPI 4	// for log_print

const int SLAVE_PRINT = 0;	/* if set to 1 then slaves will print output to a file */
void slave_print(const int level, const char *message, ...);

inline void send_int(int msg, int slave)
{
    MPI_Send(&(msg), 1, MPI_INT, slave, MASTER_ALGORITHM, MPI_COMM_WORLD);  
}
inline int receive_int_slave()
{
    MPI_Status	status;
  int msg;
  MPI_Recv(&msg, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
return msg;
}

inline void send_array_header(arraydata_t *ad, colindex_t extensioncol, int slave)
{
    colindex_t ncols = extensioncol;
    MPI_Send(&(ad->N), 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
    MPI_Send(&(ncols), 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
    MPI_Send(&(ad->ncols), 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
    MPI_Send(&(ad->strength), 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
    MPI_Send(ad->s, ad->ncols, MPI_ARRAY_T, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
}

void extend_array_mpi(arraylist_t::iterator first, int narrays, arraydata_t *ad, colindex_t extensioncol, int slave); 
void extend_array_mpi(array_t *array, arraydata_t *ad, colindex_t extensioncol, int slave);

/**
 * @brief Collect solutions from all slaves
 *
 * The solutions are collected fromn the firs slave that reports to be done.
 */
int collect_solutions_all(arraylist_t &extensions, const arraydata_t *ad);

/**
 * @brief Main function for the slaves in MPI configuration
 *
 */
int extend_slave_code(const int this_rank, OAextend const &oaextend);

/** @brief Collect all solutions from a single slave
 *
 * The solutions are collected fromn the firs slave that reports to be done.
 */
int collect_solutions_single(arraylist_t &extensions, const arraydata_t *ad);

/**
 * @brief Collect all solutions from the slave specified
 *
 * The solutions are stored in a list that is specified as an argument.
 */
void collect_solutions_slave(arraylist_t &extensions, const arraydata_t *ad, const int slave, int nsols) ;

void mpi_send_solution(array_t *array, rowindex_t nr, colindex_t nc);
void mpi_receive_solution(array_t *array, rowindex_t nr, colindex_t nc, int slave);

#endif

