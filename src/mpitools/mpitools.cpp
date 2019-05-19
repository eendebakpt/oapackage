#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "mpi.h"

#include "tools.h"
#include "mpitools.h"


/**
 * Function to print information from slaves to file (stdout from slaves is not passed to the user).
 * @param level 
 * @param message 
 */
void slave_print(const int level, const char *message, ...)
{
	const int STRBUFSIZE = 4096;
	static FILE **fid = 0;

	if (SLAVE_PRINT) {
		static int 	loglevel;
		const int MAXSLAVES = 140;
		va_list		va;

		if(fid==0) {
			fid = new FILE* [MAXSLAVES];
		}
		int rank =  MPI::COMM_WORLD.Get_rank();
		if(fid[rank]==0) {
			char buf[STRBUFSIZE];
			sprintf(buf, "slavelog-%d.txt", MPI::COMM_WORLD.Get_rank());
			fid[rank] = fopen(buf, "w");
		}

		va_start(va, message);
		char str[STRBUFSIZE];
		vsprintf(str, message, va);

		if(level < 0)			//level < 0 means set level, any message appended is printed as well
		{
			loglevel = -level;
			fprintf(fid[rank], "%s", str);
		}
		else if(level <= loglevel)	//if printing level is high enough, the message is shown
			fprintf(fid[rank], "%s", str);
		va_end(va);

		fflush(fid[rank]);
	}
}


int stop_slaves()
{
    MPI_Status	status;
    int		msg = 0, target, n_procs;

    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    for (int i = 1; i < n_procs; i++)
    {
#ifdef MPIDEBUG
        log_print(QUIET, "Stopping slave %i\n", i);
#endif
        MPI_Recv(&msg, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        //log_print(DEBUG+1, "   slave has status %d, sending MST_STOP\n", status.MPI_SOURCE);
        target = status.MPI_SOURCE;
        MPI_Send(&msg, 1, MPI_INT, status.MPI_SOURCE, MASTER_STOP, MPI_COMM_WORLD);
    }
#ifdef MPIDEBUG
    log_print(QUIET, "M:  Slaves stopped\n");
#endif

    return 0;
}



void extend_array_mpi(arraylist_t::iterator first, int narrays, arraydata_t *ad, colindex_t extensioncol, int slave)
{
    arraylist_t extensions;

#ifdef MPIDEBUG
    log_print(QUIET, "   extend_array_mpi: sending %d arrays (size %d,%d) to slave %d\n", narrays, ad->N, ad->ncols, slave);
#endif

    if (slave<0) {
        log_print(SYSTEM, "no free slaves!!!\n");
        exit(1);
    }
#ifdef MPIDEBUG
    log_print(NORMAL, "extend_array_mpi: using slave %d\n", slave);
    fflush(stdout);
#endif
    MPI_Send(&narrays, 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);

    send_array_header(ad, extensioncol, slave);

    const int base_size = ad->N*(extensioncol);
    arraylist_t::iterator	cur = first;

    for (int j=0; j<narrays; j++) {
        MPI_Send(cur->array, base_size, MPI_ARRAY_T, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
        cur++;
    }
}


void extend_array_mpi(array_t *array, arraydata_t *ad, colindex_t extensioncol, int slave)
{
    int narrays = 1;

    if (slave<0) {
        printf("no free slaves\n");
        exit(1);
    }
#ifdef MPIDEBUG
    log_print(QUIET, "extend_array_mpi: using slave %d\n", slave);
    fflush(stdout);
#endif

    MPI_Send(&narrays, 1, MPI_INT, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);

    send_array_header(ad, extensioncol, slave);

    const int base_size = ad->N*(extensioncol);
    MPI_Send(array, base_size, MPI_ARRAY_T, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD);
}


int extend_slave_code(const int this_rank, OAextend const &oaextend)
{
    MPI_Status	status;
    int msg = 0;

    MPI_Send(&msg, 1, MPI_INT, MASTER, SLAVE_READY, MPI_COMM_WORLD);
    slave_print(-DEBUG, "setting loglevel\n");

#ifdef MPIDEBUG
    slave_print(-QUIET, "setting loglevel\n");
#endif

    do
    {
        slave_print(SYSTEM, "slave %d: ready to receive command\n", this_rank);

        /* receive command from master */
        MPI_Recv(&msg, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        slave_print(SYSTEM, "slave %d: received command\n", this_rank);

        int nextensions, narrays;

        fflush(stdout);
        switch (status.MPI_TAG)
        {
        case MASTER_EXTENDJOB:
        {
            arraylist_t extensions;
            narrays=msg;
            slave_print(DEBUG, "SLAVE %d: received extend job (%d arrays)\n", this_rank, narrays);

            //slave_print(SYSTEM, "slave: receiving\n"); fflush(stdout);
            int ncols;
            int extensioncol;
            int strength;
            int N;
            MPI_Recv(&N, 1, MPI_INT, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);
            MPI_Recv(&extensioncol, 1, MPI_INT, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);
            MPI_Recv(&ncols, 1, MPI_INT, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);
            MPI_Recv(&strength, 1, MPI_INT, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);

            slave_print(DEBUG, "SLAVE %d: received extend job xxx\n", this_rank);

            array_t *s = new array_t [ncols];
            MPI_Recv(s, ncols, MPI_ARRAY_T, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);
            arraydata_t *ad = new arraydata_t(s, N, strength, ncols);
            delete [] s;

            arraylist_t		solutions;
            /* first receive all arrays, then start calculating */
            for (int i=0; i<narrays; i++) {
                slave_print(QUIET, "SLAVE %d: received extend job %d/%d: extensioncol %d ncols %d\n", this_rank, i, narrays, extensioncol, ad->ncols);

                array_t *array = create_array(ad->N, extensioncol+1); /* extended array has one row more!? */
                MPI_Recv(array, ad->N*extensioncol, MPI_ARRAY_T, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);

                array_link	tmp_extension(array, ad->N, extensioncol, -1);
                solutions.push_back(tmp_extension);
            }

            for (arraylist_t::iterator it = solutions.begin(); it != solutions.end(); it++) {
                extend_array(it->array, ad, extensioncol, extensions, oaextend);
            }
            nextensions = extensions.size();
            free_sols(solutions);

            slave_print(NORMAL, "SLAVE %d: sending results: %d arrays\n", this_rank, nextensions);
            MPI_Send(&nextensions, 1, MPI_INT, MASTER, SLAVE_SOLUTIONS, MPI_COMM_WORLD);

            for (arraylist_t::iterator it = extensions.begin(); it != extensions.end(); it++) {
                mpi_send_solution(it->array, ad->N, ad->ncols);
            }
            free_sols(extensions);
            delete ad;
        }
        break;
        case MASTER_LMCJOB:
            slave_print(SYSTEM, "SLAVE %d: received LMC job\n", this_rank);
            slave_print(SYSTEM, "SLAVE %d: LMC test not implemented as a task yet\n");
            break;

        case MASTER_WAIT:		//master is collecting solutions
            slave_print(NORMAL, "S%2i: Processor received wait signal\n", this_rank);
            MPI_Send(&msg, 1, MPI_INT, MASTER, SLAVE_READY, MPI_COMM_WORLD);

            break;
        case MASTER_STOP:
            slave_print(NORMAL, "S%2i: Slave has received STOP signal\n", this_rank);
            MPI_Send(&msg, 1, MPI_INT, MASTER, SLAVE_READY, MPI_COMM_WORLD);
            break;
        default:		//message not understood: stop
            logstream(SYSTEM) << "extend_slave_code: unknown command from MASTER: cmd " << status.MPI_TAG << std::endl;
            break;
        }
    }
    while (status.MPI_TAG != MASTER_STOP);
    return 0;
}

/**
 * @brief Collect solutions from all slaves
 *
 * The solutions are collected fromn the firs slave that reports to be done.
 */
int collect_solutions_all(arraylist_t &extensions, const arraydata_t *ad)
{
    int msg = 0;
    MPI_Status	status;

#ifdef MPIDEBUG
    log_print(QUIET, "M  : collecting all solutions\n");
#endif

    int n_procs = MPI::COMM_WORLD.Get_size();
    for (int i = 1; i < n_procs; i++)
    {
        double Tstart = get_time_ms();
        log_print(DEBUGMPI, "M  : Checking processor %i\n", i);
        MPI_Recv(&msg, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        log_print(DEBUGMPI, "M  : Message received: %d\n", status.MPI_TAG);
        double Tw = get_time_ms();

        if (status.MPI_TAG == SLAVE_SOLUTIONS)
        {
            collect_solutions_slave(extensions, ad, i, msg);
        } else {
            if (status.MPI_TAG==SLAVE_READY) {
                log_print(DEBUGMPI, "slave %d is ready\n", status.MPI_SOURCE);
            }
            else
                log_print(SYSTEM, "slave sent unknown command: %d\n", status.MPI_TAG);
        }
        MPI_Send(&msg, 1, MPI_INT, i, MASTER_WAIT, MPI_COMM_WORLD);
        log_print(NORMAL, "   collect time proc %d: %.2f (wait %.2f)\n", i, (get_time_ms()-Tstart), Tw-Tstart);
    }
    return 0;
}



void mpi_send_solution(array_t *array, rowindex_t nr, colindex_t nc)
{
#ifdef MPIDEBUG
    log_print(QUIET, "sending solution: %d %d: \n", nr, nc);
    fflush(stdout);
#endif
    MPI_Send(array, nr*nc, MPI_ARRAY_T, MASTER, MASTER_EXTENDJOB, MPI_COMM_WORLD);
}

void mpi_receive_solution(array_t *array, rowindex_t nr, colindex_t nc, int slave)
{
    MPI_Status status;
    MPI_Recv(array, nr*nc, MPI_ARRAY_T, slave, MASTER_EXTENDJOB, MPI_COMM_WORLD, &status);
}



/**
 * @brief Collect all solutions from the slave specified
 *
 * The solutions are stored in a list that is specified as an argument.
 */
void collect_solutions_slave(arraylist_t &extensions, const arraydata_t *ad, const int slave, int nsols) {
#ifdef MPIDEBUG
    log_print(QUIET, "M  : Collecting %d solution(s) from slave %i\n", nsols, slave);
#endif


    colindex_t nc = ad->ncols;

    for (int j=0; j<nsols; j++) {
		array_link	*tmp_solution;
		tmp_solution = new array_link(ad->N,nc, 0);

        log_print(DEBUG, "M  : receiving solution %d\n", j);
        mpi_receive_solution(tmp_solution->array, ad->N, nc, slave);
        extensions.push_back(*tmp_solution);

        delete tmp_solution;
    }
}

/** @brief Collect all solutions from a single slave
 *
 * The solutions are collected fromn the firs slave that reports to be done.
 */
int collect_solutions_single(arraylist_t &extensions, const arraydata_t *ad)
{
    int msg = 0;
    MPI_Status	status;

    log_print(DEBUG, "M  : collecting solutions for single slave\n");

    MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int slave = status.MPI_SOURCE;

    if (status.MPI_TAG==SLAVE_SOLUTIONS) {
        /* get solutions from slave */
        collect_solutions_slave(extensions, ad, slave, msg);
    }
    return slave;
}
