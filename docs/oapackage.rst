============================
The Orthogonal Array package
============================

:Author: P.T. Eendebak [1]_
:Date:   July, 2015

.. role:: math(raw)
   :format: html latex
..

.. role:: raw-latex(raw)
   :format: latex
..

|image|

Introduction
============

Orthogonal arrays are an important tool in the design of
experiments :raw-latex:`\cite{Hedayat1999}`. The Orthogonal Array
packagecontains functionality the generate orthogonal arrays and to
analyse their properties. The algorithms and methods in this package
have been described in :raw-latex:`\cite{Eendebak2009}`. A large
collection of arrays can be found on the Orthogonal Array
packagewebsite :raw-latex:`\cite{EendebakOA}` or on the website of Neil
Sloane :raw-latex:`\cite{Sloanewebsite}`.

Example usage
-------------

The Orthogonal Array packagecan be used to work with arrays and
calculate statistical properties. For example to calculate the
:math:`D`-efficiency and rank of a design we can use:

label=Calculate D-efficiency >>> al=oapackage.exampleArray(0) >>>
al.showarray() array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1 >>>
print(’D-efficiency D-efficiency 1.000000, rank 2 >>> print(’Generalized
wordlength pattern: Generalized wordlength pattern: (1.0, 0.0, 0.0)

Interfaces
----------

The Orthogonal Array packagehas several interfaces. First of all there
are command line tools for manipulating arrays and generating. All
functions of the package can be used from either C++ or Python. For a
restricted set of functionality also Matlab or R can be used.

Compilation and installation
----------------------------

The program has been tested using Linux and Windows (XP, Windows 7 and
Windows 8). The Python interface is available from the Python Package
Index https://pypi.python.org/pypi/OApackage/. The package can be
installed from the command line using pip:

> pip install –user OApackage

The R interface to the package is available from CRAN
http://cran.r-project.org/web/packages/oapackage/index.html.

The command line tools uses a cmake build system. From the command line
type:

> mkdir -p build; cd build > cmake .. > make > make install

This creates the command line utilities and a C++ library. To compile
the Python interface using Linux use

> python setup.py build > python setup.py install –user

The Python interface requires Numpy :raw-latex:`\cite{NumPy}`,
Matplotlib :raw-latex:`\cite{Matplotlib}` and Swig. The code has been
tested with Python 2.7, 3.4 and 3.5.

Using Windows start Cygwin or the Visual Studio command prompt. From the
package source directory run:

> python setup.py bdist\_wininst

This creates a binary installer package.

License
-------

The code is available under a BSD style license, see the file LICENSE
for details. If you use this code or any of the results, please cite
this program as follows:

*Complete Enumeration of Pure-Level and Mixed-Level Orthogonal Arrays*,
P.T. Eendebak, E.D. Schoen, M.V.M. Nguyen, Volume 18, Issue 2, pages
123-140, 2010.

Acknowledgements
----------------

The code and ideas for this package have been contributed by Eric
Schoen, Ruben Snepvangers, Vincent Brouerius van Nidek, Alan
Vazquez-Alcocer and Pieter Thijs Eendebak.

The Orthogonal Array package
============================

An orthogonal array (OA) of strength :math:`{t}`, :math:`{N}` runs and
:math:`{n}` factors at :math:`{s}` levels is an :math:`{N}\times {n}`
array of :math:`0,
\ldots,({s}-1)`-valued symbols such that for every :math:`{t}` columns
every :math:`{t}`-tuple occurs equally
often :raw-latex:`\cite{Rao1947}`. The set of all OAs with given
strength, runs and levels is denoted by
:math:`{\operatorname{OA}({N}; {t}; {s}^{n})}`. The OAs are represented
by arrays (in column-major form).

Data structures
---------------

The package contains several data structures. Here we describe the main
structures and their use.

 ``array_link``
    The structure containing an orthogonal array is called the
    \|array\_link\| structure. Lists of arrays are stored in the
    \|arraylist\_t\| object, which as a \|std::deque\| container.

\|arrayfile\_t\|
    This is an objct that allows for reading and writing of arrays to
    disk.

\|arraydata\_t\|
    The structure describing a certain class of orthogonal arrays or
    designs.

\|array\_transformation\_t\|
    This describes a transformation of an array. This includes the row-,
    column- and level-permutations.

Representing arrays
~~~~~~~~~~~~~~~~~~~

The structure containing an orthogonal array is called the
\|array\_link\| structure. It consists of a specified number of rows and
columns, the data (integer values) and an index.

[def:array\_link]

::

    struct array_link

::

    {
        //! Number of rows in array
        rowindex_t n_rows;
        //! Number of columns in array
        colindex_t n_columns;
        //! Index number
        int index;
        //! Pointer to an array data
        array_t* array;

        /// Constructor functions
        array_link();
        array_link(rowindex_t nrows, colindex_t ncols, int index);
        ~array_link();
        array_link(const array_link &);

    public:
        /// print an array to output stream
        friend std::ostream &operator<<(std::ostream &, const array_link &A);

        /// print array to stdout
        void showarray() const;

        // manipulation of arrays
        
        /// return array with selected column removed
        array_link deleteColumn(int index) const;

        /// return array with first n columns selected
        array_link selectFirstColumns(int n) const;

        /// return array with last n columns selected
        array_link selectLastColumns(int n) const;

        /// select columns from an array
        array_link selectColumns(const std::vector<int> c) const;

        /// return transposed array
        array_link transposed() const;

        // statistical properties of the array

        ...

In the Python interface the arraylink object can be indexed just as
normal arrays. It is also possible to return a Numpy array. The
\|array\_link\| object implements to Python array interface, so most
opertations from packages such as Numpy work on the \|array\_link\|
object.

label=Array representation in Python >>> import oapackage >>>
al=oapackage.exampleArray(0) >>> al.showarray() array: 0 0 0 0 0 1 0 1 1
0 1 0 1 1 1 1 >>> al[2,1] 1L >>> X=al.getarray() >>> X array([[0, 0],
[0, 0], [0, 1], [0, 1], [1, 0], [1, 0], [1, 1], [1, 1]], dtype=int32)

Reading and writing arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~

Reading and writing arrays to disk can be done with the \|arrayfile\_t\|
class. For example:

label=Write an array to disk >>> import oapackage >>>
al=oapackage.exampleArray() >>> af=oapackage.arrayfile\_t(’test.oa’,
al.n\_rows, al.n\_columns) >>> af.append\_array(al) >>> print(af) file
test.oa: 8 rows, 2 columns, 1 arrays, mode text, nbits 8 >>>
af.closefile()

The arrays can be written in text or binary format. For more details on
the file format see Section [section:fileformats]. The header of the
\|arrayfile\_t\| class is listed below.

::

    struct arrayfile_t
    {

    public:
        std::string filename;
        int iscompressed;
        int nrows;
        int ncols;

        /// number of bits used when storing an array
        int nbits;

        /// file mode, can be ATEXT or ABINARY
        arrayfilemode_t mode;
        /// file opened for reading or writing
        afilerw_t rwmode;

        int narrays;
        int narraycounter;

    public:

        /// open existing array file
        arrayfile_t(const std::string fname, int verbose = 1);
        /// open new array file for writing
        arrayfile_t(const std::string fname, int nrows, int ncols,
                     int narrays=-1, arrayfilemode_t m = ATEXT, int nb = 8);
        /// destructor function, closes all filehandles
        ~arrayfile_t();

        /// close the array file
        void closefile();
        /// return true if file is open
        int isopen() const;
        /// seek to specified array position
        int seek(int pos);
        /// read array and return index
        int read_array(array_link &a);
        /// return true if the file has binary format
        bool isbinary() const;
        /// append arrays to the file
        int append_arrays(const arraylist_t &arrays, int startidx);
        /// append a single array to the file
        void append_array(const array_link &a, int specialindex=-1);

        ...
        
    }

Array transformations
~~~~~~~~~~~~~~~~~~~~~

Transformations of (orthogonal) arrays consist of row permutations,
level permutations and level transformations. A transformation is
represented by the \|array\_transformation\_t\| object.

For a given transformation the column permutations are applied first,
then the level permutations and finally the row permutations. The level-
and column permutations are not commutative.

[code:arraytransformationt]

::

    class array_transformation_t
    {
    public:
        rowperm_t   rperm;      /// row permutation
        colperm_t   colperm;    /// column permutation
        levelperm_t *lperms;    /// level permutations
        const arraydata_t *ad;  /// type of array

    public:
        array_transformation_t ( const arraydata_t *ad );
        array_transformation_t ( );     /// default constructor
        array_transformation_t ( const array_transformation_t  &at );   
        array_transformation_t & operator= ( const array_transformation_t &at );    
        ~array_transformation_t();  /// destructor

        /// show the array transformation
        void show() const;

        /// return true if the transformation is equal to the identity
        bool isIdentity() const;

        /// return the inverse transformation
        array_transformation_t inverse() const;

        /// return the transformation to the identity transformation
        void reset();

        /// initialize to a random transformation
        void randomize();

        /// initialize with a random column transformation
        void randomizecolperm();

        /// apply transformation to an array_link object
        array_link apply ( const array_link &al ) const;

        /// composition operator. the transformations are applied from the left
        array_transformation_t operator*(const array_transformation_t b);
        
        ...

Classes of arrays
~~~~~~~~~~~~~~~~~

The \|arraydata\_t\| object represents data about a class of orthogonal
arrays, e.g. the class :math:`{\operatorname{OA}(N; t; s^k)}`.

::

    struct arraydata_t
    {
        rowindex_t N;   /** number of runs */
        array_t *s; /** pointer to levels of the array */
        colindex_t ncols; /** total number of columns (factors) in the design */
        colindex_t strength;    /** strength of the design */

        ordering_t  order; /** Ordering used for arrays */

    public:
        /// create new arraydata_t object
        arraydata_t(std::vector<int> s, rowindex_t N_, colindex_t t, colindex_t nc);
        arraydata_t(carray_t *s_, rowindex_t N_, colindex_t t, colindex_t nc);
        arraydata_t(const arraydata_t &adp);
        
        ...
        
        /// return true if the array is of mixed type
        bool ismixed() const;
        /// return true if the array is a 2-level array
        bool is2level() const;
        /// set column group equal to that of a symmetry group
        void set_colgroups(const symmetry_group &sg);
            /// return random array from the class
        array_link randomarray ( int strength = 0, int ncols=-1 ) const;

    }

File formats
------------

The Orthogonal Array packagestored orthogonal arrays in a custom file
format. There is a text format with is easily readable by humans and a
binary format with is faster to process and memory efficient.

Plain text array files
~~~~~~~~~~~~~~~~~~~~~~

Arrays are stored in plain text files with extension .oa. The first line
contains the number of columns, the number of rows and the number of
arrays (or -1 if the number of arrays is not specified). Then for each
array a single line with the index of the array, followed by N lines
containing the array.

A typical example of a text file would be:

[formatcom=,fontsize=,frame=single,framesep=0.8ex,rulecolor=] 5 8 1 1 0
0 0 0 0 0 0 0 1 1 0 1 1 0 0 0 1 1 1 1 1 0 1 0 1 1 0 1 1 0 1 1 0 0 1 1 1
0 1 0 -1

This file contains exactly 1 array with 8 rows and 5 columns.

Binary array files
~~~~~~~~~~~~~~~~~~

Every binary file starts with a header, which has the following format:

[fontsize=] [INT32] 65 (magic identifier) [INT32] b: Format: number of
bits per number. Currently supported are 1 and 8 [INT32] N: number of
rows [INT32] k: kumber of columns [INT32] Number of arrays (can be -1 if
unknown) [INT32] Binary format number: 1001: normal, 1002: binary diff,
1003: binary diff zero [INT32] Reserved integer [INT32] Reserved integer

The normal binary format has the following format. For each array (the
number is specified in the header):

[INT32] Index [Nxk elements] The elements contain b bits

If the number of bits per number is 1 (e.g. a 2-level array) then the
data is padded with zeros to a multiple of 64 bits. The data of the
array is stored in column-major order. The binary file format allows for
random access reading and writing. The \|binary diff\| and \|binary diff
zero\| formats are special formats.

A binary array file can be compressed using gzip. Most tools in the
Orthogonal Array packagecan read these compressed files transparently.
Writing to compressed array files is not supported at the moment.

Data files
~~~~~~~~~~

The analysis tool (\|oaanalyse\|) writes data to disk in binary format.
The format is consists of a binary header:

[FLOAT64] Magic number 30397995; [FLOAT64] Magic number 12224883;
[FLOAT64] nc: Number of rows [FLOAT64] nr: Number of columns

After the header there follow \|nc\*nr [FLOAT64]\| values.

Statistical properties of an array
----------------------------------

Most properties of an array can be calculated using the \|array\_link\|
object. The interface is listed below.

::

    struct array_link
    {
        ... 
        
    public:
        // statistical properties of the array

        /// calculate D-efficiency
        double Defficiency() const;
        /// calculate main effect robustness (or Ds optimality)
        double DsEfficiency(int verbose=0) const;
        /// calculate A-efficiency
        double Aefficiency() const;
        /// calculate E-efficiency
        double Eefficiency() const;
        /// calculate rank of array
        int rank() const;
        /// calculate generalized wordlength pattern
        std::vector<double> GWLP() const;
        /// return true if the array is a foldover array
        bool foldover() const;
        /// calculate centered L2 discrepancy
        double CL2discrepancy() const;
        /// Calculate the projective estimation capacity sequence
        std::vector<double> PECsequence() const;
    }

The :math:`D`-efficiency, :math:`A`-efficiency and :math:`E`-efficiency
are calculated by calculating the SVD of the second order interaction
matrix. The efficiencies can then be calculated using the eigenvalues of
the SVD. For the definition of the :math:`D`-, :math:`A`- and
:math:`E`-efficiency see Definition [definition:Aefficiency]. For the
rank of a matrix the LU decomposition of the matrix is calculated using
the Eigen package :raw-latex:`\cite{eigenweb}`.

[:math:`{{\color{darkblue}D}}`-efficiency and average VIF] Let
:math:`{{\color{darkblue}X}}` be an :math:`N\times k` :math:`2`-factor
array with second order model :math:`{F({{\color{darkblue}X}})}`. Then
we define the :math:`{{\color{darkblue}D}}`-efficiency and the average
variance inflation factor as

.. math::

   \begin{aligned}
   {{\color{darkblue}D}({{\color{darkblue}X}})} = \left( \det {F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}\right)^{1/m} / N , 
   \label{formula:Defficiency} \\
   {\mathrm{VIF}({{\color{darkblue}X}})} = N \operatorname{tr}\left( \frac{1}{ {F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}} \right) /m . \label{formula:VIF} \end{aligned}

 [definition:Aefficiency] [definition:A value] [definition:B value] The
matrix :math:`{F({{\color{darkblue}X}})}^T {F({{\color{darkblue}X}})}`
is called the information matrix. Let
:math:`\lambda_1, \ldots, \lambda_m` be the eigenvalues of the
information matrix. Then the :math:`{{\color{darkblue}E}}`-efficiency of
a matrix is [definition:Eefficiency] defined as

.. math::

   \begin{aligned}
   {{\color{darkblue}E}({{\color{darkblue}X}})} = \min_j \lambda_j .
   \label{formula:E-efficiency}\end{aligned}

 Note that in terms of the eigenvalues we have
:math:`{{\color{darkblue}D}(X)} = (\prod_j \lambda_j)^{1/m} / N` and
:math:`{\mathrm{VIF}(X)} = N (\sum_j \lambda_j^{-1})/m`.

The :math:`D_s`-effiency is the main effect robustness, see the appendix
in :raw-latex:`\cite{Schoen2010}` for more details.

Calculation of D-optimal designs
--------------------------------

D-optimal designs can be calculated with the function \|Doptimize\|.
This function uses a coordinate exchange algorithm to generate designs
with good properties for the :math:`D`-efficiency.

An example script with Python to generate optimal designs with 40 runs
and 7 factors is shown below.

label=Doptimize >>> N=40; s=2; k=7; >>>
arrayclass=oapackage.arraydata\_t(s, N, 0, k) >>> print(’We generate
optimal designs with: We generate optimal designs with: arrayclass: N
40, k 7, strength 0, s 2,2,2,2,2,2,2, order 0. >>> alpha=[1,2,0] >>>
method=oapackage.DOPTIM\_UPDATE >>> scores, dds, designs, ngenerated =
oapackage.Doptimize(arrayclass, nrestarts=40, optimfunc=alpha,
selectpareto=True) Doptim: optimization class 40.2-2-2-2-2-2-2
Doptimize: iteration 0/40 Doptimize: iteration 39/40 Doptim: done (8
arrays, 0.6 [s]) >>> print(’Generated 8 designs, the best D-efficiency
is 0.9098

The parameters of the function are documented in the code.

To calculate properties of designs we can use the following functions.
For :math:`D`-efficiencies we can use

::

    std::vector<double> array_link::Defficiencies ( int verbose ) const;

to calculate the :math:`D`-, :math:`D_s`- and :math:`D_1`-efficiency.
For details see :raw-latex:`\cite{EendebakSO}`.

The projective estimation capacity (PEC) sequence
from :raw-latex:`\cite{loeppky2004ranking}` can be calculated with:

::

    std::vector<double> PECsequence(const array_link &al, int verbose=1);

.. figure:: images/motivating-40-d-2-2-2-2-2-2-2-scatterplot-ndata2.png
   :alt: Scatterplot for the :math:`D`-efficiency and
   :math:`D_s`-efficiency for generated designs in
   :math:`{\operatorname{OA}(40; 2; 2^7)}`. The Pareto optimal designs
   are colored, while the non-Pareto optimal designs are grey. For
   reference the strength-3 orthogonal array with highest D-efficiency
   is also included in the plot.

   Scatterplot for the :math:`D`-efficiency and :math:`D_s`-efficiency
   for generated designs in :math:`{\operatorname{OA}(40; 2; 2^7)}`. The
   Pareto optimal designs are colored, while the non-Pareto optimal
   designs are grey. For reference the strength-3 orthogonal array with
   highest D-efficiency is also included in the plot.

GWLP and J-characteristics
--------------------------

From an \|array\_link\| object we can calculate the generalized
worldlength patterns :raw-latex:`\cite{Xu2001}`, :math:`F`-values and
:math:`J`-characteristics.

label=Calculate GWLP and :math:`F`-values >>>
al=oapackage.exampleArray(1) >>> al.showarray() array: 0 0 0 0 0 0 0 0 0
0 0 0 0 1 1 0 0 1 0 1 0 1 0 1 0 0 1 1 0 0 0 1 1 1 1 0 1 1 1 1 1 0 0 1 1
1 0 1 0 1 1 0 1 1 0 1 0 1 1 0 1 1 0 0 1 1 1 0 0 1 1 1 0 1 0 1 1 1 0 0

>>> g=al.GWLP() >>> print(’GWLP: GWLP: (1.0, 0.0, 0.0, 1.0, 1.0, 0.0)
>>> print(’F3-value: F3-value: (4, 6) >>> print(’F4-value: F4-value: (1,
4) >>> print(’J3-characteristis: J3-characteristis: (8, 8, 0, 0, 0, 8,
0, 8, 0, 0)

Reduction to canonical form
---------------------------

If we introduce an ordering on the set of arrays, then for each
isomorphism class of arrays the minimal element defines a unique
canonical form. The Orthogonal Array packagecontains functions to reduce
any orhogonal array to canonical form with respect to some ordering. The
default ordering for arrays is the lexicographic ordering in
columns :raw-latex:`\cite{Eendebak2009}`. An alternative ordering is the
delete-one-factor projection ordering as described
in :raw-latex:`\cite{EendebakDOF}`.

::

    {}
    /// reduce an array to canonical form using delete-1-factor ordering
    array_link reduceLMCform(const array_link &al);

    /// reduce an array to canonical form using delete-1-factor ordering
    array_link reduceDOPform(const array_link &al);

Another approach to generation of canonical forms for designs is to use
graph-isomophism packages such as
Nauty :raw-latex:`\cite{nautyI, nautyII}` or
Bliss :raw-latex:`\cite{Bliss2007}`. To reduce a general graph to Nauty
canonical form one can use ``reduceNauty``. For orthogonal arrays we can
encode the array structure as a graph. The reduction can then be done
with ``reduceOAnauty``.

label=Reduce a design to normal form using Nauty >>> al =
oapackage.exampleArray(0).randomperm() >>> al.showarray() array: 0 1 1 1
1 0 0 0 0 0 1 1 1 0 0 1 >>> t=oapackage.reduceOAnauty(al, verbose=0) >>>
t.show() array transformation: N 8 column permutation: 0,1 level perms:
0,1 0,1 row permutation: 3,4,0,7,2,6,1,5 >>> alr=t.apply(al) >>>
alr.showarray() array: 0 0 0 0 0 1 0 1 1 0 1 0 1 1 1 1

Generation of arrays
--------------------

A list of arrays in LMC form can be extended to a list of arrays in LMC
form with one additional column. Details for the algorithm are described
in :raw-latex:`\cite{Eendebak2009}`.

The main function for array extension is the following:

::

    {label=C++ interface}
    /// extend a list of arrays
    arraylist_t & extend_arraylist(arraylist_t & alist, arraydata_t &fullad, 
                OAextend const &oaextendoptions);

Here \|fullad\| is the structure describing the type of arrays and
\|oaxextendoptions\| contains various options for the algorithm.

A typical session could be:

label=Extend an array >>> N=8; ncols=3; >>>
arrayclass=oapackage.arraydata\_t(2, N, 2, ncols) >>>
al=arrayclass.create\_root() >>> al.showarray() array: 0 0 0 0 0 1 0 1 1
0 1 0 1 1 1 1 >>> >>> alist=oapackage.extend\_array(al, arrayclass) >>>
for al in alist: ... al.showarray()

array: 0 0 0 0 0 0 0 1 1 0 1 1 1 0 1 1 0 1 1 1 0 1 1 0

array: 0 0 0 0 0 1 0 1 0 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1

Conference designs
------------------

An :math:`n\times k` conference design is an :math:`N\times k` matrix
with entries 0, -1, +1 such that i) in each column the symbol 0 occurs
exactly one time and ii) all columns are orthogonal to each other. A
more detailed description is given
in :raw-latex:`\cite{wiki:ConferenceMatrix}`.

label=Generate conference designs with 8 rows >>> import oapackage >>>
ctype=oapackage.conference\_t(N=8, k=8) >>> >>> al =
ctype.create\_root\_three() >>> al.showarray() array: 0 1 1 1 0 -1 1 1 0
1 1 1 1 1 -1 1 -1 1 1 -1 1 1 -1 -1 >>> l4=oapackage.extend\_conference (
[al], ctype, verbose=0) >>> l5=oapackage.extend\_conference ( l4, ctype,
verbose=0) >>> l6=oapackage.extend\_conference ( l5, ctype, verbose=0)
>>> >>> print(’number of non-isomorphic conference designs: number of
non-isomorphic conference designs: 11

MD5 sums
--------

To check data structures on disk the packages includes functions to
generate MD5 sums. These are:

::

    /// calculate md5 sum of a data block in memory
    std::string md5(void *data, int numbytes);
    /// calculate md5 sum of a file on disk
    std::string md5(const std::string filename);

Command line interface
======================

Included in the packages are several command line tools. For each tool
help can be obtained from the command line by using the switch \|-h\|.
These are:

\|oainfo\|
    This program reads Orthogonal Array packagedata files and reports
    the contents of the files. For example:

     eendebakpt:math:` oainfo result-8.2-2-2-2.oa
    Orthogonal Array package 1.8.7
    oainfo: reading 1 file(s)
    file result-8.2-2-2.oa: 8 rows, 3 columns, 2 arrays, mode text, nbits 0
    ~eendebakpt`

\|oacat\|
    Show the contents of a file with orthogonal arrays for a data file.

\|oacheck\|
    Check or reduce an array to canonical form.

\|oaextendsingle\|
    Extend a set of arrays in LMC form with one or more columns.

\|oacat\|
    Show the contents of an array file or data file.

    Usage: oacat [OPTIONS] [FILES]

\|oajoin\|
    Read one or more files from disk and join all the array files into a
    single list.

    Orthogonal Arrays 1.8.7 For more details see the files README.txt
    and LICENSE.txt

    Orthonal Array Join: join several array files into a single file
    Usage: oajoin [OPTIONS] [FILES]

    -h –help Prints this help -s –sort Sort the arrays -l –latex Output
    with LaTeX format -o [FILE] –output [FILE] Output prefix (default:
    standard output) -f [FORMAT] Output format (TEXT, BINARY (default),
    D (binary difference) )

\|oasplit\|
    Takes a single array file as input and splits the arrays to a
    specified number of output files.

\|oapareto\|
    Calculates the set of Pareto optimal arrays in a file with arrays.

\|oaanalyse\|
    Calculates various statistics of arrays in a file. The statistics
    are described in section [section:properties].

.. [1]
   Corresponding author. E-mail: pieter.eendebak@gmail.com. Address:
   University of Antwerp, Dept. of Mathematics, Statistics, and
   Actuarial Sciences, Prinsstraat 13, 2000 Antwerp, Belgium.

.. |image| image:: images/oaimage-18_2-3-3-3-3-3-n17.png
