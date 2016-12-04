
// File: index.xml

// File: classAnyOption.xml
%feature("docstring") AnyOption "C++ includes: anyoption.h ";

%feature("docstring")  AnyOption::AnyOption "AnyOption::AnyOption()
";

%feature("docstring")  AnyOption::AnyOption "AnyOption::AnyOption(int
maxoptions) ";

%feature("docstring")  AnyOption::AnyOption "AnyOption::AnyOption(int
maxoptions, int maxcharoptions) ";

%feature("docstring")  AnyOption::~AnyOption "AnyOption::~AnyOption()
";

%feature("docstring")  AnyOption::setCommandPrefixChar "void
AnyOption::setCommandPrefixChar(char _prefix) ";

%feature("docstring")  AnyOption::setCommandLongPrefix "void
AnyOption::setCommandLongPrefix(char *_prefix) ";

%feature("docstring")  AnyOption::setFileCommentChar "void
AnyOption::setFileCommentChar(char _comment) ";

%feature("docstring")  AnyOption::setFileDelimiterChar "void
AnyOption::setFileDelimiterChar(char _delimiter) ";

%feature("docstring")  AnyOption::useCommandArgs "void
AnyOption::useCommandArgs(int _argc, char **_argv) ";

%feature("docstring")  AnyOption::useFileName "void
AnyOption::useFileName(const char *_filename) ";

%feature("docstring")  AnyOption::noPOSIX "void AnyOption::noPOSIX()
";

%feature("docstring")  AnyOption::setVerbose "void
AnyOption::setVerbose() ";

%feature("docstring")  AnyOption::setOption "void
AnyOption::setOption(const char *opt_string) ";

%feature("docstring")  AnyOption::setOption "void
AnyOption::setOption(char opt_char) ";

%feature("docstring")  AnyOption::setOption "void
AnyOption::setOption(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::setFlag "void
AnyOption::setFlag(const char *opt_string) ";

%feature("docstring")  AnyOption::setFlag "void
AnyOption::setFlag(char opt_char) ";

%feature("docstring")  AnyOption::setFlag "void
AnyOption::setFlag(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::setCommandOption "void
AnyOption::setCommandOption(const char *opt_string) ";

%feature("docstring")  AnyOption::setCommandOption "void
AnyOption::setCommandOption(char opt_char) ";

%feature("docstring")  AnyOption::setCommandOption "void
AnyOption::setCommandOption(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::setCommandFlag "void
AnyOption::setCommandFlag(const char *opt_string) ";

%feature("docstring")  AnyOption::setCommandFlag "void
AnyOption::setCommandFlag(char opt_char) ";

%feature("docstring")  AnyOption::setCommandFlag "void
AnyOption::setCommandFlag(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::setFileOption "void
AnyOption::setFileOption(const char *opt_string) ";

%feature("docstring")  AnyOption::setFileOption "void
AnyOption::setFileOption(char opt_char) ";

%feature("docstring")  AnyOption::setFileOption "void
AnyOption::setFileOption(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::setFileFlag "void
AnyOption::setFileFlag(const char *opt_string) ";

%feature("docstring")  AnyOption::setFileFlag "void
AnyOption::setFileFlag(char opt_char) ";

%feature("docstring")  AnyOption::setFileFlag "void
AnyOption::setFileFlag(const char *opt_string, char opt_char) ";

%feature("docstring")  AnyOption::processOptions "void
AnyOption::processOptions() ";

%feature("docstring")  AnyOption::processCommandArgs "void
AnyOption::processCommandArgs() ";

%feature("docstring")  AnyOption::processCommandArgs "void
AnyOption::processCommandArgs(int max_args) ";

%feature("docstring")  AnyOption::processFile "bool
AnyOption::processFile() ";

%feature("docstring")  AnyOption::processCommandArgs "void
AnyOption::processCommandArgs(int _argc, char **_argv) ";

%feature("docstring")  AnyOption::processCommandArgs "void
AnyOption::processCommandArgs(int _argc, char **_argv, int max_args)
";

%feature("docstring")  AnyOption::processFile "bool
AnyOption::processFile(const char *_filename) ";

%feature("docstring")  AnyOption::getValue "char*
AnyOption::getValue(const char *_option) ";

%feature("docstring")  AnyOption::getFlag "bool
AnyOption::getFlag(const char *_option) ";

%feature("docstring")  AnyOption::getValue "char*
AnyOption::getValue(char _optchar) ";

%feature("docstring")  AnyOption::getFlag "bool
AnyOption::getFlag(char _optchar) ";

%feature("docstring")  AnyOption::getDoubleValue "double
AnyOption::getDoubleValue(const char *_option, double defaultvalue) ";

%feature("docstring")  AnyOption::getDoubleValue "double
AnyOption::getDoubleValue(char _optchar, double defaultvalue) ";

%feature("docstring")  AnyOption::getIntValue "int
AnyOption::getIntValue(char _optchar, int defaultvalue) ";

%feature("docstring")  AnyOption::getIntValue "int
AnyOption::getIntValue(const char *_option, int defaultvalue) ";

%feature("docstring")  AnyOption::getLongValue "long
AnyOption::getLongValue(const char *_option, long defaultvalue) ";

%feature("docstring")  AnyOption::getStringValue "const char*
AnyOption::getStringValue(const char *_option, const char
*defaultvalue) ";

%feature("docstring")  AnyOption::getStringValue "const char*
AnyOption::getStringValue(char _optchar, const char *defaultvalue) ";

%feature("docstring")  AnyOption::printUsage "void
AnyOption::printUsage() ";

%feature("docstring")  AnyOption::printAutoUsage "void
AnyOption::printAutoUsage() ";

%feature("docstring")  AnyOption::addUsage "void
AnyOption::addUsage(const std::string line) ";

%feature("docstring")  AnyOption::addUsage "void
AnyOption::addUsage(const char *line) ";

%feature("docstring")  AnyOption::printHelp "void
AnyOption::printHelp() ";

%feature("docstring")  AnyOption::autoUsagePrint "void
AnyOption::autoUsagePrint(bool flag) ";

%feature("docstring")  AnyOption::getArgc "int AnyOption::getArgc()
";

%feature("docstring")  AnyOption::getArgv "char*
AnyOption::getArgv(int index) ";

%feature("docstring")  AnyOption::hasOptions "bool
AnyOption::hasOptions() ";


// File: structarray__link.xml
%feature("docstring") array_link "

Wrapper class for an array.

The array_link struct is a struct that represents an arrays. Copying
of array links is done with shallow copy or deep copy depending on
compile time options!

C++ includes: arraytools.h ";

%feature("docstring")  array_link::array_link "array_link::array_link()

Constructor functions. ";

%feature("docstring")  array_link::array_link "array_link::array_link(rowindex_t nrows, colindex_t ncols, int index)
";

%feature("docstring")  array_link::array_link "array_link::array_link(rowindex_t nrows, colindex_t ncols, int index,
carray_t *data) ";

%feature("docstring")  array_link::~array_link "array_link::~array_link() ";

%feature("docstring")  array_link::array_link "array_link::array_link(const array_link &) ";

%feature("docstring")  array_link::array_link "array_link::array_link(Eigen::MatrixXd &m) ";

%feature("docstring")  array_link::clone "array_link
array_link::clone() const ";

%feature("docstring")  array_link::showarray "void
array_link::showarray() const

print array to stdout ";

%feature("docstring")  array_link::showarraycompact "void
array_link::showarraycompact() const

print array to stdout ";

%feature("docstring")  array_link::showproperties "void
array_link::showproperties() const

print array properties to stdout ";

%feature("docstring")  array_link::is2level "bool
array_link::is2level() const

return true if the arra is a 2-level array (e.g. only contains 0 and
1) ";

%feature("docstring")  array_link::is_conference "bool
array_link::is_conference() const

return true if the array is a +1,0, -1 valued array ";

%feature("docstring")  array_link::deleteColumn "array_link
array_link::deleteColumn(int index) const

return array with selected column removed ";

%feature("docstring")  array_link::selectFirstRows "array_link
array_link::selectFirstRows(int n) const

return array with first n rows ";

%feature("docstring")  array_link::selectFirstColumns "array_link
array_link::selectFirstColumns(int n) const

return array with first n columns selected ";

%feature("docstring")  array_link::selectLastColumns "array_link
array_link::selectLastColumns(int n) const

return array with last n columns selected ";

%feature("docstring")  array_link::selectColumns "array_link
array_link::selectColumns(const std::vector< int >c) const

select columns from an array ";

%feature("docstring")  array_link::selectColumns "array_link
array_link::selectColumns(int c) const

select single column from an array ";

%feature("docstring")  array_link::setColumn "void
array_link::setColumn(int c, const std::vector< int >v)

set a column of the array to the given vector ";

%feature("docstring")  array_link::setColumn "void
array_link::setColumn(int c, const std::vector< signed char >v) ";

%feature("docstring")  array_link::transposed "array_link
array_link::transposed() const

return transposed array ";

%feature("docstring")  array_link::Defficiency "double
array_link::Defficiency() const

calculate D-efficiency ";

%feature("docstring")  array_link::DsEfficiency "double
array_link::DsEfficiency(int verbose=0) const

calculate main effect robustness (or Ds-optimality) ";

%feature("docstring")  array_link::Defficiencies "std::vector< double
> array_link::Defficiencies(int verbose=0, int addDs0=0) const

calculate D-efficiency, calculate main effect robustness (or Ds-
optimality) and D1-efficiency ";

%feature("docstring")  array_link::VIFefficiency "double
array_link::VIFefficiency() const ";

%feature("docstring")  array_link::Aefficiency "double
array_link::Aefficiency() const

calculate A-efficiency ";

%feature("docstring")  array_link::Eefficiency "double
array_link::Eefficiency() const

calculate E-efficiency ";

%feature("docstring")  array_link::Fvalues "std::vector< int >
array_link::Fvalues(int jj) const

Calculate F-values of a 2-level matrix. ";

%feature("docstring")  array_link::FvaluesConference "std::vector<
int > array_link::FvaluesConference(int jj) const

Calculate F-values of a conference design. ";

%feature("docstring")  array_link::Jcharacteristics "std::vector< int
> array_link::Jcharacteristics(int jj=4) const

Calculate J-characteristics of matrix (the values are signed) ";

%feature("docstring")  array_link::PECsequence "std::vector< double >
array_link::PECsequence() const

Calculate the projective estimation capacity sequence. ";

%feature("docstring")  array_link::rank "int array_link::rank() const

calculate rank of array ";

%feature("docstring")  array_link::GWLP "std::vector< double >
array_link::GWLP(int truncate=1, int verbose=0) const

calculate generalized wordlength pattern ";

%feature("docstring")  array_link::strength "int
array_link::strength() const

calculate strength of an array ";

%feature("docstring")  array_link::foldover "bool
array_link::foldover() const

return true if the array is a foldover array ";

%feature("docstring")  array_link::min "array_t array_link::min()
const ";

%feature("docstring")  array_link::max "array_t array_link::max()
const ";

%feature("docstring")  array_link::CL2discrepancy "double
array_link::CL2discrepancy() const

calculate centered L2 discrepancy

The method is from \"A connection between uniformity and aberration in
regular fractions of two-level factorials\", Fang and Mukerjee, 2000
";

%feature("docstring")  array_link::randomperm "array_link
array_link::randomperm() const

apply a random permutation of rows, columns and levels ";

%feature("docstring")  array_link::randomcolperm "array_link
array_link::randomcolperm() const

apply a random permutation of columns ";

%feature("docstring")  array_link::randomrowperm "array_link
array_link::randomrowperm() const

apply a random permutation of row ";

%feature("docstring")  array_link::getModelMatrix "MatrixFloat
array_link::getModelMatrix(int order, int intercept=1) const

This function calculates Helmert contrasts for the factors of an input
design. implementation from code written by Eric Schoen, Dept. of
Applied Economics, University of Antwerp, Belgium ";

%feature("docstring")  array_link::deepcopy "array_link &
array_link::deepcopy(const array_link &rhs) ";

%feature("docstring")  array_link::shallowcopy "array_link &
array_link::shallowcopy(const array_link &rhs) ";

%feature("docstring")  array_link::equalsize "int
array_link::equalsize(const array_link &rhs) const ";

%feature("docstring")  array_link::atfast "const array_t&
array_link::atfast(const rowindex_t r, const colindex_t c) const

get element from array, no error checking, inline version ";

%feature("docstring")  array_link::atfast "array_t&
array_link::atfast(const rowindex_t r, const colindex_t c)

get element from array, no error checking, inline version ";

%feature("docstring")  array_link::_at "array_t array_link::_at(const
rowindex_t, const colindex_t) const ";

%feature("docstring")  array_link::_at "array_t array_link::_at(const
int index) const

get element at specified position, no error checking ";

%feature("docstring")  array_link::at "array_t array_link::at(const
rowindex_t, const colindex_t) const

get element at specified position, no error checking ";

%feature("docstring")  array_link::at "array_t array_link::at(const
int index) const

get element at specified position ";

%feature("docstring")  array_link::at "array_t& array_link::at(const
rowindex_t, const colindex_t)

get element at specified position ";

%feature("docstring")  array_link::setconstant "void
array_link::setconstant(array_t val)

get element at specified position

set all elements in the array to a value ";

%feature("docstring")  array_link::setvalue "void
array_link::setvalue(int row, int col, int val)

set value of an array ";

%feature("docstring")  array_link::setvalue "void
array_link::setvalue(int row, int col, double val) ";

%feature("docstring")  array_link::_setvalue "void
array_link::_setvalue(int row, int col, int val)

set value of an array ";

%feature("docstring")  array_link::negateRow "void
array_link::negateRow(rowindex_t r)

set value of an array, no error checking!

multiply a row by -1 ";

%feature("docstring")  array_link::show "void array_link::show()
const

print information about array ";

%feature("docstring")  array_link::showstr "std::string
array_link::showstr() const ";

%feature("docstring")  array_link::md5 "std::string array_link::md5()
const

return md5 sum of array representation (as represented with 32bit int
datatype in memory) ";

%feature("docstring")  array_link::columnEqual "bool
array_link::columnEqual(int rl, const array_link &rhs, int rr) const

return true if two columns are equal ";

%feature("docstring")  array_link::firstColumnDifference "int
array_link::firstColumnDifference(const array_link &A) const

return index of first different column ";

%feature("docstring")  array_link::firstDiff "bool
array_link::firstDiff(const array_link &A, int &r, int &c, int
verbose=1) const

calculate row and column index of first difference between two arrays

The difference is according to the column-major ordering. ";

%feature("docstring")  array_link::create_root "void
array_link::create_root(const arraydata_t &ad)

create root in arraylink ";

%feature("docstring")  array_link::nonzero_fraction "double
array_link::nonzero_fraction() const ";

%feature("docstring")  array_link::clear "void array_link::clear()

fill array with zeros ";

%feature("docstring")  array_link::getarraydata "void
array_link::getarraydata(int *pymat1, int n) ";

%feature("docstring")  array_link::setarraydata "void
array_link::setarraydata(const numtype *tmp, int n)

internal function ";

%feature("docstring")  array_link::setarraydata "void
array_link::setarraydata(std::vector< int >tmp, int n)

special method for SWIG interface ";

%feature("docstring")  array_link::setarraydata "void
array_link::setarraydata(std::vector< numtype > tmp, int n)

internal function ";

%feature("docstring")  array_link::maxelement "array_t
array_link::maxelement() const ";

%feature("docstring")  array_link::setcolumn "void
array_link::setcolumn(int c, const array_link &al, int sc=0)

set column to values ";

%feature("docstring")  array_link::array_link "array_link::array_link(const array_link &, const std::vector< int >
&colperm) ";

%feature("docstring")  array_link::array_link "array_link::array_link(const array_t *array, rowindex_t nrows,
colindex_t ncols, int index=0) ";

%feature("docstring")  array_link::array_link "array_link::array_link(const array_t *array, rowindex_t nrows,
colindex_t ncolsorig, colindex_t ncols, int index) ";

%feature("docstring")  array_link::array_link "array_link::array_link(const std::vector< int > &v, rowindex_t nrows,
colindex_t ncols, int index=0) ";

%feature("docstring")  array_link::init "void
array_link::init(rowindex_t r, colindex_t c) ";

%feature("docstring")  array_link::row_symmetry_group "symmetry_group
array_link::row_symmetry_group() const

return the row_symmetry group of an array ";

%feature("docstring")  array_link::reduceLMC "array_link
array_link::reduceLMC() const

return the LMC form of the array ";

%feature("docstring")  array_link::reduceDOP "array_link
array_link::reduceDOP() const

return the delete-one-factor-projection form of the array ";

%feature("docstring")  array_link::getEigenMatrix "MatrixFloat
array_link::getEigenMatrix() const

return the array as an Eigen matrix ";

%feature("docstring")  array_link::columnGreater "int
array_link::columnGreater(int c1, const array_link &rhs, int c2) const

return true of specified column is smaller than column in another
array ";

%feature("docstring")  array_link::showarrayS "std::string
array_link::showarrayS() const ";


// File: classarray__transformation__t.xml
%feature("docstring") array_transformation_t "

Contains a transformation of an array.

Contains an array transformation. The transformation consists of
column, row and level permutations. The level and column permutations
are not commutative (since the level permutations are tied to a
particular column). We apply the column permutations first.

C++ includes: arraytools.h ";

%feature("docstring")  array_transformation_t::array_transformation_t
"array_transformation_t::array_transformation_t(const arraydata_t
*ad)

type of array ";

%feature("docstring")  array_transformation_t::array_transformation_t
"array_transformation_t::array_transformation_t(const arraydata_t
&ad) ";

%feature("docstring")  array_transformation_t::array_transformation_t
"array_transformation_t::array_transformation_t() ";

%feature("docstring")  array_transformation_t::array_transformation_t
"array_transformation_t::array_transformation_t(const
array_transformation_t &at)

default constructor ";

%feature("docstring")  array_transformation_t::~array_transformation_t
"array_transformation_t::~array_transformation_t()

assignment operator ";

%feature("docstring")  array_transformation_t::show "void
array_transformation_t::show() const

destructor

show the array transformation ";

%feature("docstring")  array_transformation_t::isIdentity "bool
array_transformation_t::isIdentity() const

return true if the transformation is equal to the identity ";

%feature("docstring")  array_transformation_t::inverse "array_transformation_t array_transformation_t::inverse() const

return the inverse transformation ";

%feature("docstring")  array_transformation_t::reset "void
array_transformation_t::reset()

return the transformation to the identity transformation ";

%feature("docstring")  array_transformation_t::randomize "void
array_transformation_t::randomize()

initialize to a random transformation ";

%feature("docstring")  array_transformation_t::randomizecolperm "void
array_transformation_t::randomizecolperm()

initialize with a random column permutation ";

%feature("docstring")  array_transformation_t::randomizerowperm "void
array_transformation_t::randomizerowperm()

initialize with a random row permutation ";

%feature("docstring")  array_transformation_t::apply "array_link
array_transformation_t::apply(const array_link &al) const

apply transformation to an array_link object ";

%feature("docstring")  array_transformation_t::applygeneric "array_link array_transformation_t::applygeneric(const array_link &al)
const

apply transformation to an array_link object ";

%feature("docstring")  array_transformation_t::apply "void
array_transformation_t::apply(array_t *sourcetarget)

apply transformation to an array ";

%feature("docstring")  array_transformation_t::apply "void
array_transformation_t::apply(const array_t *source, array_t *target)
const

apply transformation to an array ";

%feature("docstring")  array_transformation_t::print_transformed "void array_transformation_t::print_transformed(carray_t *source) const

apply transformation and show resulting array ";

%feature("docstring")  array_transformation_t::show "void
array_transformation_t::show(std::ostream &out) const ";

%feature("docstring")  array_transformation_t::rowperm "std::vector<
int > array_transformation_t::rowperm() const ";

%feature("docstring")  array_transformation_t::colperm "std::vector<
int > array_transformation_t::colperm() const

return the row permutation of the transformation ";

%feature("docstring")  array_transformation_t::lvlperm "std::vector<
int > array_transformation_t::lvlperm(int c) const

return the column permutation of the transformation ";

%feature("docstring")  array_transformation_t::setrowperm "void
array_transformation_t::setrowperm(std::vector< int > rp)

return the level permutations of the transformation ";

%feature("docstring")  array_transformation_t::setcolperm "void
array_transformation_t::setcolperm(std::vector< int > colperm) ";

%feature("docstring")  array_transformation_t::setlevelperm "void
array_transformation_t::setlevelperm(int colindex, std::vector< int >
lvlperm) ";


// File: structarraydata__t.xml
%feature("docstring") arraydata_t "

Contains properties of the design (number of rows, columns, levels)

Constructor: arrayclass = arraydata_t(s, N, strength,ncolumns)

C++ includes: arraytools.h ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t(array_t s, rowindex_t N, colindex_t strength,
colindex_t ncols)

create new arraydata_t object ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t(const std::vector< int >s, rowindex_t N,
colindex_t strength, colindex_t ncols) ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t(const array_t *s_, rowindex_t N, colindex_t
strength, colindex_t ncols) ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t(const arraydata_t &adp) ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t(const arraydata_t *adp, colindex_t newncols)

copy constructor ";

%feature("docstring")  arraydata_t::arraydata_t "arraydata_t::arraydata_t()

copy constructor ";

%feature("docstring")  arraydata_t::~arraydata_t "arraydata_t::~arraydata_t()

dummy constructor ";

%feature("docstring")  arraydata_t::ismixed "bool
arraydata_t::ismixed() const

destructor

return true if the array is of mixed type ";

%feature("docstring")  arraydata_t::is2level "bool
arraydata_t::is2level() const

return true if the array is a 2-level array ";

%feature("docstring")  arraydata_t::randomarray "array_link
arraydata_t::randomarray(int strength=0, int ncols=-1) const

return random array from the class. this operation is only valid for
strength 0 or 1 ";

%feature("docstring")  arraydata_t::writeConfigFile "void
arraydata_t::writeConfigFile(const char *filename) const

Write file with design of OA.

Parameters:
-----------

file:

ad:  Arraydata structure to write ";

%feature("docstring")  arraydata_t::idstr "std::string
arraydata_t::idstr() const ";

%feature("docstring")  arraydata_t::idstrseriesfull "std::string
arraydata_t::idstrseriesfull() const ";

%feature("docstring")  arraydata_t::fullidstr "std::string
arraydata_t::fullidstr(int series=0) const ";

%feature("docstring")  arraydata_t::latexstr "std::string
arraydata_t::latexstr(int cmd=0, int series=0) const ";

%feature("docstring")  arraydata_t::reduceColumns "arraydata_t
arraydata_t::reduceColumns(int k) ";

%feature("docstring")  arraydata_t::showstr "std::string
arraydata_t::showstr() const ";

%feature("docstring")  arraydata_t::show "void arraydata_t::show(int
verbose=1) const ";

%feature("docstring")  arraydata_t::complete_arraydata "void
arraydata_t::complete_arraydata() ";

%feature("docstring")  arraydata_t::lmc_overflow_check "void
arraydata_t::lmc_overflow_check() const

check whether the LMC calculation will overflow ";

%feature("docstring")  arraydata_t::complete_arraydata_fixlast "void
arraydata_t::complete_arraydata_fixlast() ";

%feature("docstring")  arraydata_t::complete_arraydata_splitn "void
arraydata_t::complete_arraydata_splitn(int ns) ";

%feature("docstring")  arraydata_t::set_colgroups "void
arraydata_t::set_colgroups(const std::vector< int >splits) ";

%feature("docstring")  arraydata_t::set_colgroups_jj "void
arraydata_t::set_colgroups_jj(const symmetry_group &sg, int jj) ";

%feature("docstring")  arraydata_t::set_colgroups "void
arraydata_t::set_colgroups(const symmetry_group &sg)

set column group equal to that of a symmetry group ";

%feature("docstring")  arraydata_t::show_colgroups "void
arraydata_t::show_colgroups() const ";

%feature("docstring")  arraydata_t::calcoaindex "void
arraydata_t::calcoaindex(colindex_t strength) ";

%feature("docstring")  arraydata_t::create_root "array_link
arraydata_t::create_root() const

return the root array for the class ";

%feature("docstring")  arraydata_t::getfactorlevel "int
arraydata_t::getfactorlevel(int idx) const ";

%feature("docstring")  arraydata_t::getS "std::vector< int >
arraydata_t::getS() const ";

%feature("docstring")  arraydata_t::reset_strength "void
arraydata_t::reset_strength(colindex_t t)

Reset strength of arraydata.

Parameters:
-----------

t:  ";

%feature("docstring")  arraydata_t::get_col_group "colindex_t
arraydata_t::get_col_group(const colindex_t col) const

Return index of the column group for a column. ";


// File: structarrayfile_1_1arrayfile__t.xml
%feature("docstring") arrayfile::arrayfile_t "

Structure for reading or writing a file with arrays.

The format of array files is described in the file FORMAT.txt

C++ includes: arraytools.h ";

%feature("docstring")  arrayfile::arrayfile_t::arrayfile_t "arrayfile::arrayfile_t::arrayfile_t()

default constructor ";

%feature("docstring")  arrayfile::arrayfile_t::arrayfile_t "arrayfile::arrayfile_t::arrayfile_t(const std::string fname, int
verbose=1)

open existing array file ";

%feature("docstring")  arrayfile::arrayfile_t::arrayfile_t "arrayfile::arrayfile_t::arrayfile_t(const std::string fname, int
nrows, int ncols, int narrays=-1, arrayfilemode_t m=ATEXT, int nb=8)

open new array file for writing ";

%feature("docstring")  arrayfile::arrayfile_t::~arrayfile_t "arrayfile::arrayfile_t::~arrayfile_t()

destructor function, closes all filehandles ";

%feature("docstring")  arrayfile::arrayfile_t::createfile "void
arrayfile::arrayfile_t::createfile(const std::string fname, int nrows,
int ncols, int narrays=-1, arrayfilemode_t m=ATEXT, int nb=8)

close current file and open a new file for writing ";

%feature("docstring")  arrayfile::arrayfile_t::closefile "void
arrayfile::arrayfile_t::closefile()

close the array file ";

%feature("docstring")  arrayfile::arrayfile_t::isopen "int
arrayfile::arrayfile_t::isopen() const

return true if file is open ";

%feature("docstring")  arrayfile::arrayfile_t::seek "int
arrayfile::arrayfile_t::seek(int pos)

seek to specified array position ";

%feature("docstring")  arrayfile::arrayfile_t::read_array "int
arrayfile::arrayfile_t::read_array(array_link &a)

read array and return index ";

%feature("docstring")  arrayfile::arrayfile_t::readnext "array_link
arrayfile::arrayfile_t::readnext()

read next array from the file ";

%feature("docstring")  arrayfile::arrayfile_t::readarrays "arraylist_t arrayfile::arrayfile_t::readarrays(int nmax=NARRAYS_MAX,
int verbose=1)

read set of array from the file ";

%feature("docstring")  arrayfile::arrayfile_t::flush "void
arrayfile::arrayfile_t::flush()

flush any open file pointer ";

%feature("docstring")  arrayfile::arrayfile_t::isbinary "bool
arrayfile::arrayfile_t::isbinary() const

return true if the file has binary format ";

%feature("docstring")  arrayfile::arrayfile_t::append_arrays "int
arrayfile::arrayfile_t::append_arrays(const arraylist_t &arrays, int
startidx=-1)

append list of arrays to the file ";

%feature("docstring")  arrayfile::arrayfile_t::append_array "void
arrayfile::arrayfile_t::append_array(const array_link &a, int
specialindex=-1)

append a single array to the file ";

%feature("docstring")  arrayfile::arrayfile_t::swigcheck "int
arrayfile::arrayfile_t::swigcheck() const ";

%feature("docstring")  arrayfile::arrayfile_t::showstr "std::string
arrayfile::arrayfile_t::showstr() const ";

%feature("docstring")  arrayfile::arrayfile_t::pos "size_t
arrayfile::arrayfile_t::pos() const ";

%feature("docstring")  arrayfile::arrayfile_t::hasrandomaccess "bool
arrayfile::arrayfile_t::hasrandomaccess() const

return true of the file format has random access mode ";

%feature("docstring")  arrayfile::arrayfile_t::updatenumbers "void
arrayfile::arrayfile_t::updatenumbers() ";

%feature("docstring")  arrayfile::arrayfile_t::read_array "int
arrayfile::arrayfile_t::read_array(array_t *array, const int nrows,
const int ncols)

read array and return index ";

%feature("docstring")  arrayfile::arrayfile_t::finisharrayfile "void
arrayfile::arrayfile_t::finisharrayfile() ";

%feature("docstring")  arrayfile::arrayfile_t::setVerbose "void
arrayfile::arrayfile_t::setVerbose(int v) ";

%feature("docstring")  arrayfile::arrayfile_t::getnbits "int
arrayfile::arrayfile_t::getnbits() ";


// File: structarraysymmetry.xml
%feature("docstring") arraysymmetry "

class representing array symmetries

C++ includes: lmc.h ";

%feature("docstring")  arraysymmetry::~arraysymmetry "arraysymmetry::~arraysymmetry() ";

%feature("docstring")  arraysymmetry::arraysymmetry "arraysymmetry::arraysymmetry(const dyndata_t *dyndata) ";

%feature("docstring")  arraysymmetry::arraysymmetry "arraysymmetry::arraysymmetry(const arraysymmetry &rhs)

copy constructor ";

%feature("docstring")  arraysymmetry::N "int arraysymmetry::N() const
";

%feature("docstring")  arraysymmetry::show "void
arraysymmetry::show() const ";


// File: structarraywriter__t.xml
%feature("docstring") arraywriter_t "

structure to write arrays to disk, thread safe

C++ includes: arraytools.h ";

%feature("docstring")  arraywriter_t::arraywriter_t "arraywriter_t::arraywriter_t()

verbosity level ";

%feature("docstring")  arraywriter_t::~arraywriter_t "arraywriter_t::~arraywriter_t() ";

%feature("docstring")  arraywriter_t::flush "void
arraywriter_t::flush() ";

%feature("docstring")  arraywriter_t::writeArray "void
arraywriter_t::writeArray(const array_link &A)

write a single array to disk ";

%feature("docstring")  arraywriter_t::writeArray "void
arraywriter_t::writeArray(const arraylist_t &lst) ";

%feature("docstring")  arraywriter_t::initArrayFiles "void
arraywriter_t::initArrayFiles(const arraydata_t &ad, int kstart, const
std::string prefix, arrayfilemode_t mode=ABINARY_DIFF) ";

%feature("docstring")  arraywriter_t::nArraysWritten "int
arraywriter_t::nArraysWritten() const

return the total number arrays written to disk ";

%feature("docstring")  arraywriter_t::closeafiles "void
arraywriter_t::closeafiles() ";


// File: classCandidateGenerator.xml
%feature("docstring") CandidateGenerator "

Class to generate candidate extensions with caching.

C++ includes: conference.h ";

%feature("docstring")  CandidateGenerator::CandidateGenerator "CandidateGenerator::CandidateGenerator(const array_link &al, const
conference_t &ct) ";

%feature("docstring")  CandidateGenerator::generateConfCandidates "const std::vector<cperm>&
CandidateGenerator::generateConfCandidates(const array_link &al, int
kz) const

generate candidates with caching this method uses j2 filtering ";

%feature("docstring")  CandidateGenerator::showCandidates "void
CandidateGenerator::showCandidates() const ";

%feature("docstring")  CandidateGenerator::updateLastValid "void
CandidateGenerator::updateLastValid(int lv) const ";


// File: classCandidateGeneratorDouble.xml
%feature("docstring") CandidateGeneratorDouble "

Class to generate candidate extensions with caching.

C++ includes: conference.h ";

%feature("docstring")
CandidateGeneratorDouble::CandidateGeneratorDouble "CandidateGeneratorDouble::CandidateGeneratorDouble(const array_link
&al, const conference_t &ct) ";

%feature("docstring")
CandidateGeneratorDouble::generateDoubleConfCandidates "const
std::vector<cperm>&
CandidateGeneratorDouble::generateDoubleConfCandidates(const
array_link &al) const

generate candidates with caching this method uses symmetry inflation,
assumes j1=0 and j2=0 ";

%feature("docstring")  CandidateGeneratorDouble::showCandidates "void
CandidateGeneratorDouble::showCandidates() const ";


// File: classgfx_1_1Compare.xml
%feature("docstring") gfx::Compare "C++ includes: timsort.hpp ";

%feature("docstring")  gfx::Compare::Compare "gfx::Compare< Value,
LessFunction >::Compare(LessFunction f) ";

%feature("docstring")  gfx::Compare::Compare "gfx::Compare< Value,
LessFunction >::Compare(const Compare< value_type, func_type > &other)
";

%feature("docstring")  gfx::Compare::lt "bool gfx::Compare< Value,
LessFunction >::lt(value_type x, value_type y) ";

%feature("docstring")  gfx::Compare::le "bool gfx::Compare< Value,
LessFunction >::le(value_type x, value_type y) ";

%feature("docstring")  gfx::Compare::gt "bool gfx::Compare< Value,
LessFunction >::gt(value_type x, value_type y) ";

%feature("docstring")  gfx::Compare::ge "bool gfx::Compare< Value,
LessFunction >::ge(value_type x, value_type y) ";

%feature("docstring")  gfx::Compare::less_function "func_type&
gfx::Compare< Value, LessFunction >::less_function() ";


// File: structconf__candidates__t.xml
%feature("docstring") conf_candidates_t "

structure to cache a list of candidate extensions

C++ includes: conference.h ";

%feature("docstring")  conf_candidates_t::info "void
conf_candidates_t::info(int verbose=1) const ";


// File: structconference__extend__t.xml
%feature("docstring") conference_extend_t "

Helper structure.

C++ includes: conference.h ";

%feature("docstring")  conference_extend_t::combine "cperm
conference_extend_t::combine(int i, int j) const ";

%feature("docstring")  conference_extend_t::nExtensions "size_t
conference_extend_t::nExtensions() const ";

%feature("docstring")  conference_extend_t::getarrays "arraylist_t
conference_extend_t::getarrays(const array_link al) const ";


// File: structconference__options.xml
%feature("docstring") conference_options "C++ includes: conference.h
";

%feature("docstring")  conference_options::conference_options "conference_options::conference_options(int maxpos=-1) ";


// File: classconference__t.xml
%feature("docstring") conference_t "

Structure representing the type of conference designs.

C++ includes: conference.h ";

%feature("docstring")  conference_t::conference_t "conference_t::conference_t(int N, int k, int j1zero)

for the double conference type matrices

create new conference_t object ";

%feature("docstring")  conference_t::conference_t "conference_t::conference_t(const conference_t &rhs) ";

%feature("docstring")  conference_t::create_root "array_link
conference_t::create_root() const

create the unique representative of the 2 column design (for
conference matrices) ";

%feature("docstring")  conference_t::create_root_three "array_link
conference_t::create_root_three() const

create the unique representative of the 3 column design ";

%feature("docstring")  conference_t::createDconferenceRootArrays "arraylist_t conference_t::createDconferenceRootArrays() const ";

%feature("docstring")  conference_t::addRootArrays "void
conference_t::addRootArrays(arraylist_t &lst) const ";

%feature("docstring")  conference_t::__repr__ "std::string
conference_t::__repr__() const

return string representation of the object ";


// File: classconference__transformation__t.xml
%feature("docstring") conference_transformation_t "

Contains a transformation of a conference matrix.

Contains an array transformation. The transformation consists of
column permutations, row permutations and sign switches for both the
rows and columns.

The sign switches and the permutations are not commutative. We apply
the permutations first and then the sign flips.

C++ includes: arraytools.h ";

%feature("docstring")
conference_transformation_t::conference_transformation_t "conference_transformation_t::conference_transformation_t() ";

%feature("docstring")
conference_transformation_t::conference_transformation_t "conference_transformation_t::conference_transformation_t(int nrows,
int ncols)

default constructor ";

%feature("docstring")
conference_transformation_t::conference_transformation_t "conference_transformation_t::conference_transformation_t(const
array_link &al) ";

%feature("docstring")
conference_transformation_t::conference_transformation_t "conference_transformation_t::conference_transformation_t(const
conference_transformation_t &T) ";

%feature("docstring")  conference_transformation_t::show "void
conference_transformation_t::show(int verbose=1) const

show the array transformation ";

%feature("docstring")  conference_transformation_t::isIdentity "bool
conference_transformation_t::isIdentity() const

return true if the transformation is equal to the identity ";

%feature("docstring")  conference_transformation_t::inverse "conference_transformation_t conference_transformation_t::inverse()
const

return the inverse transformation ";

%feature("docstring")  conference_transformation_t::reset "void
conference_transformation_t::reset()

return the transformation to the identity transformation ";

%feature("docstring")  conference_transformation_t::randomize "void
conference_transformation_t::randomize()

initialize to a random transformation ";

%feature("docstring")  conference_transformation_t::randomizecolperm "void conference_transformation_t::randomizecolperm()

initialize with a random column permutation ";

%feature("docstring")  conference_transformation_t::randomizerowperm "void conference_transformation_t::randomizerowperm()

initialize with a random row permutation ";

%feature("docstring")  conference_transformation_t::randomizecolflips
"void conference_transformation_t::randomizecolflips()

initialize with random col switches ";

%feature("docstring")  conference_transformation_t::randomizerowflips
"void conference_transformation_t::randomizerowflips()

initialize with random row switches ";

%feature("docstring")  conference_transformation_t::apply "array_link
conference_transformation_t::apply(const array_link &al) const

apply transformation to an array_link object ";

%feature("docstring")  conference_transformation_t::setrowperm "void
conference_transformation_t::setrowperm(std::vector< int >rp) ";

%feature("docstring")  conference_transformation_t::setcolperm "void
conference_transformation_t::setcolperm(std::vector< int >cp) ";


// File: structcounter__t.xml
%feature("docstring") counter_t "

structure to count and show number of arrays generated, the structure
is thread safe

C++ includes: evenodd.h ";

%feature("docstring")  counter_t::counter_t "counter_t::counter_t(int
n) ";

%feature("docstring")  counter_t::addNfound "void
counter_t::addNfound(int col, int num) ";

%feature("docstring")  counter_t::nArrays "long counter_t::nArrays()
const ";

%feature("docstring")  counter_t::addNumberFound "void
counter_t::addNumberFound(int n, int k) ";

%feature("docstring")  counter_t::clearNumberFound "void
counter_t::clearNumberFound() ";

%feature("docstring")  counter_t::addNumberFound "void
counter_t::addNumberFound(const counter_t &de) ";

%feature("docstring")  counter_t::showcountscompact "void
counter_t::showcountscompact() const

show information about the number of arrays found ";

%feature("docstring")  counter_t::showcounts "void
counter_t::showcounts(const arraydata_t &ad) const

show information about the number of arrays found ";

%feature("docstring")  counter_t::showcounts "void
counter_t::showcounts(const char *str, int first, int last) const

show information about the number of arrays found ";


// File: classDconferenceFilter.xml
%feature("docstring") DconferenceFilter "

class to filter designs

C++ includes: conference.h ";

%feature("docstring")  DconferenceFilter::DconferenceFilter "DconferenceFilter::DconferenceFilter(const array_link &_als, int
filtersymm_, int filterj2_, int filterj3_=1) ";

%feature("docstring")  DconferenceFilter::filterList "std::vector<cperm> DconferenceFilter::filterList(const std::vector<
cperm > &lst, int verbose=0) const

filter a list of cperms using the filter method ";

%feature("docstring")  DconferenceFilter::filterListJ2last "std::vector<cperm> DconferenceFilter::filterListJ2last(const
std::vector< cperm > &lst) const ";

%feature("docstring")  DconferenceFilter::filterListZero "std::vector<cperm> DconferenceFilter::filterListZero(const
std::vector< cperm > &lst) const

filter a list of cperms using the filterZero method ";

%feature("docstring")  DconferenceFilter::filter "bool
DconferenceFilter::filter(const cperm &c) const

return True of the extension satisfies all checks ";

%feature("docstring")  DconferenceFilter::filterJpartial "bool
DconferenceFilter::filterJpartial(const cperm &c, int r) const

filter on partial column (only last col) ";

%feature("docstring")  DconferenceFilter::filterJ "bool
DconferenceFilter::filterJ(const cperm &c, int j2start=0) const

return True of the extension satisfies all J-characteristic checks ";

%feature("docstring")  DconferenceFilter::filterJlast "bool
DconferenceFilter::filterJlast(const cperm &c, int j2start=0) const

return True of the extension satisfies all J-characteristic checks for
the last columns ";

%feature("docstring")  DconferenceFilter::filterReason "bool
DconferenceFilter::filterReason(const cperm &c) const

return True of the extension satisfies all checks ";

%feature("docstring")  DconferenceFilter::filterJ3 "bool
DconferenceFilter::filterJ3(const cperm &c) const

return True of the candidate satisfies the J3 check ";

%feature("docstring")  DconferenceFilter::filterJ3s "bool
DconferenceFilter::filterJ3s(const cperm &c, int idxstart) const

return True of the candidate satisfies the J3 check for specified
pairs ";

%feature("docstring")  DconferenceFilter::filterJ3r "bool
DconferenceFilter::filterJ3r(const cperm &c) const

return True of the candidate satisfies the J3 check ";

%feature("docstring")  DconferenceFilter::filterJ3inline "bool
DconferenceFilter::filterJ3inline(const cperm &c) const

return True of the candidate satisfies the J3 check ";

%feature("docstring")  DconferenceFilter::filterSymmetry "bool
DconferenceFilter::filterSymmetry(const cperm &c) const

return True of the candidate satisfies the symmetry check ";

%feature("docstring")  DconferenceFilter::filterJ2 "bool
DconferenceFilter::filterJ2(const cperm &c) const

return True of the candidate extension satisfies the J2 check ";

%feature("docstring")  DconferenceFilter::filterJ2last "bool
DconferenceFilter::filterJ2last(const cperm &c) const

return True of the candidate extension satisfies the J2 check for the
last column of the array checked against ";

%feature("docstring")  DconferenceFilter::filterZero "bool
DconferenceFilter::filterZero(const cperm &c) const

return True of the candidate extension satisfies the zero ";


// File: structdepth__extend__sub__t.xml
%feature("docstring") depth_extend_sub_t "

Helper structure for dynamic extension

In this structure we keep track of pointers to valid column extensions

C++ includes: evenodd.h ";

%feature("docstring")  depth_extend_sub_t::depth_extend_sub_t "depth_extend_sub_t::depth_extend_sub_t(int nn=0) ";

%feature("docstring")  depth_extend_sub_t::resize "void
depth_extend_sub_t::resize(int nn) ";

%feature("docstring")  depth_extend_sub_t::n "size_t
depth_extend_sub_t::n() const ";

%feature("docstring")  depth_extend_sub_t::updateExtensionPointers "std::vector< int > depth_extend_sub_t::updateExtensionPointers(int
extcol) ";

%feature("docstring")  depth_extend_sub_t::initialize "arraylist_t
depth_extend_sub_t::initialize(const arraylist_t &alist, const
arraydata_t &adf, const OAextend &oaextend)

initialize the new list of extension columns ";

%feature("docstring")  depth_extend_sub_t::selectArraysZ "arraylist_t
depth_extend_sub_t::selectArraysZ(const arraylist_t &alist) const

select the arrays with are LMC and hence need to be written to disk ";

%feature("docstring")  depth_extend_sub_t::selectArraysXX "arraylist_t depth_extend_sub_t::selectArraysXX(const array_link &al,
const arraylist_t &elist) const ";

%feature("docstring")  depth_extend_sub_t::info "void
depth_extend_sub_t::info() const ";


// File: structdepth__extend__t.xml
%feature("docstring") depth_extend_t "

Helper structure for dynamic extension.

This structure allows for writing the generated arrays to disk. It
also contains functions to print progress of the extension.

Multiple copies of this class are made, but they all share the same
counter_t and arraywriter_t object. Also t0 and tp are shared

C++ includes: evenodd.h ";

%feature("docstring")  depth_extend_t::depth_extend_t "depth_extend_t::depth_extend_t(const arraydata_t *ad_, double
_logtime=10000000, int _discardJ5=-1) ";

%feature("docstring")  depth_extend_t::depth_extend_t "depth_extend_t::depth_extend_t(const depth_extend_t &de) ";

%feature("docstring")  depth_extend_t::~depth_extend_t "depth_extend_t::~depth_extend_t() ";

%feature("docstring")  depth_extend_t::show "void
depth_extend_t::show() ";

%feature("docstring")  depth_extend_t::setNarraysMax "void
depth_extend_t::setNarraysMax(long n) ";

%feature("docstring")  depth_extend_t::maxArrayCheck "void
depth_extend_t::maxArrayCheck()

HACK ";

%feature("docstring")  depth_extend_t::showsearchpath "void
depth_extend_t::showsearchpath(int depth) const ";

%feature("docstring")  depth_extend_t::showprogress "bool
depth_extend_t::showprogress(int showtime=1, int depth=0, int
forcelog=0)

show information about the progress of the loop ";

%feature("docstring")  depth_extend_t::info "void
depth_extend_t::info() const ";

%feature("docstring")  depth_extend_t::setposition "void
depth_extend_t::setposition(int k, int c, int m, int extensioncols=-1,
int goodextensioncols=-1)

set the position in the dextend structure ";

%feature("docstring")  depth_extend_t::setpositionGEC "void
depth_extend_t::setpositionGEC(int k, int goodextensioncols)

set the position in the dextend structure ";


// File: structdepth__extensions__storage__t.xml
%feature("docstring") depth_extensions_storage_t "

Helper structure for the even-odd depth extension.

C++ includes: evenodd.h ";

%feature("docstring")  depth_extensions_storage_t::resize "void
depth_extensions_storage_t::resize(size_t s) ";

%feature("docstring")  depth_extensions_storage_t::set "void
depth_extensions_storage_t::set(int ai, const arraylist_t &goodarrays,
const arraylist_t &extension_column_list, depth_alg_t depthalg, const
depth_extend_sub_t &dextendsub) ";


// File: structdepth__path__t.xml
%feature("docstring") depth_path_t "

structure containing current position in search tree

C++ includes: evenodd.h ";

%feature("docstring")  depth_path_t::depth_path_t "depth_path_t::depth_path_t() ";

%feature("docstring")  depth_path_t::updatePositionGEC "void
depth_path_t::updatePositionGEC(int k, int goodextensioncols) ";

%feature("docstring")  depth_path_t::updatePosition "void
depth_path_t::updatePosition(int k, int c, int m, int extensioncols,
int goodextensioncols) ";

%feature("docstring")  depth_path_t::show "void
depth_path_t::show(int depth, int maxentries=8) const ";

%feature("docstring")  depth_path_t::init "void
depth_path_t::init(int ncols, int _depthstart=9) ";


// File: structdextend__t.xml
%feature("docstring") dextend_t "

Helper structure for dynamic extension.

C++ includes: extend.h ";

%feature("docstring")  dextend_t::dextend_t "dextend_t::dextend_t()
";

%feature("docstring")  dextend_t::resize "void dextend_t::resize(int
nn) ";

%feature("docstring")  dextend_t::DefficiencyFilter "void
dextend_t::DefficiencyFilter(double Dfinal, int k, int kfinal, double
Lmax, int verbose=1)

perform filtering using D-efficiency ";

%feature("docstring")  dextend_t::filterArrays "std::vector<int>
dextend_t::filterArrays(const array_link &al, const arraylist_t
&earrays, arraylist_t &earraysout, std::vector< std::vector< double >
> &edata, int verbose=1)

filter the arrays based on values in filter ";


// File: structDoptimReturn.xml
%feature("docstring") DoptimReturn "

Structure containing results of the Doptimize function.

C++ includes: Deff.h ";


// File: structdyndata__t.xml
%feature("docstring") dyndata_t "

Contains dynamic data of an array.

The dynamic data are used in the inner loops of the LMC algorithm. In
particular they keep track of the current row ordering and column
permutation. By not applying these transformations to the array we can
save calculation time.

We try to prevent copying the object, so it is re-used at different
levels in the algorithm. N: static col: changes at each column level

rowsort: changes at each column level, used mainly in non-root stage

colperm: changes at all levels See:   arraydata_t

C++ includes: lmc.h ";

%feature("docstring")  dyndata_t::dyndata_t "dyndata_t::dyndata_t(int
N, int col=0) ";

%feature("docstring")  dyndata_t::dyndata_t "dyndata_t::dyndata_t(const dyndata_t *dd) ";

%feature("docstring")  dyndata_t::dyndata_t "dyndata_t::dyndata_t(const dyndata_t &) ";

%feature("docstring")  dyndata_t::~dyndata_t "dyndata_t::~dyndata_t()
";

%feature("docstring")  dyndata_t::show "void dyndata_t::show() const
";

%feature("docstring")  dyndata_t::reset "void dyndata_t::reset() ";

%feature("docstring")  dyndata_t::setColperm "void
dyndata_t::setColperm(const colperm_t perm, int n) ";

%feature("docstring")  dyndata_t::setColperm "void
dyndata_t::setColperm(const larray< colindex_t > &perm) ";

%feature("docstring")  dyndata_t::setColperm "void
dyndata_t::setColperm(const std::vector< colindex_t > &perm) ";

%feature("docstring")  dyndata_t::initsymmetry "void
dyndata_t::initsymmetry(const arraysymmetry &arraysymm, const symmdata
&sd, int ncols)

initialize the rowsort structure from an arraysymmetry object ";

%feature("docstring")  dyndata_t::getRowperm "void
dyndata_t::getRowperm(rowpermtypelight &rp) const

set lightweight row permutation ";

%feature("docstring")  dyndata_t::getRowperm "void
dyndata_t::getRowperm(rowperm_t &rperm) const

get row permutation ";

%feature("docstring")  dyndata_t::getRowperm "rowpermtypelight
dyndata_t::getRowperm() const

return lightweight row permutation ";

%feature("docstring")  dyndata_t::getColperm "colpermtypelight
dyndata_t::getColperm() const

return column permutation ";

%feature("docstring")  dyndata_t::getColperm "void
dyndata_t::getColperm(colpermtypelight &cp) const

set column permutation ";

%feature("docstring")  dyndata_t::allocrowsortl "void
dyndata_t::allocrowsortl()

allocate lightweight rowsort structure ";

%feature("docstring")  dyndata_t::deleterowsortl "void
dyndata_t::deleterowsortl() ";

%feature("docstring")  dyndata_t::initrowsortl "void
dyndata_t::initrowsortl() ";

%feature("docstring")  dyndata_t::rowsortl2rowsort "void
dyndata_t::rowsortl2rowsort()

helper function ";

%feature("docstring")  dyndata_t::copydata "void
dyndata_t::copydata(const dyndata_t &dd) ";


// File: structextend__data__t.xml
%feature("docstring") extend_data_t "

Contains static data for the extend loop.

C++ includes: strength.h ";

%feature("docstring")  extend_data_t::extend_data_t "extend_data_t::extend_data_t(const arraydata_t *ad, colindex_t extcol)
";

%feature("docstring")  extend_data_t::~extend_data_t "extend_data_t::~extend_data_t() ";


// File: classindexsort.xml
%feature("docstring") indexsort "

Class to sort data without moving the data in memory.

The data is sorted by using a list of indices. A stable sort is being
used.

C++ includes: mathtools.h ";

%feature("docstring")  indexsort::indexsort "indexsort::indexsort(int
nn) ";

%feature("docstring")  indexsort::indexsort "indexsort::indexsort(const std::deque< Type > &vals)

Constructor for deque class. ";

%feature("docstring")  indexsort::indexsort "indexsort::indexsort(const std::vector< Type > &vals)

Constructor for vector class. ";

%feature("docstring")  indexsort::init "void indexsort::init(const
std::deque< Type > &vals)

initialize sorting structure with specified values ";

%feature("docstring")  indexsort::init "void indexsort::init(const
std::vector< Type > &vals)

initialize sorting structure with specified values ";

%feature("docstring")  indexsort::sort "void indexsort::sort(const
Type *vals)

sort values and store the indices ";

%feature("docstring")  indexsort::sort "void indexsort::sort(const
std::vector< Type > &vals)

sort values and store the indices ";

%feature("docstring")  indexsort::sort "void indexsort::sort(const
std::deque< Type > &vals)

sort values and store the indices ";

%feature("docstring")  indexsort::sortdescending "void
indexsort::sortdescending(const std::vector< Type > &vals)

sort values and store the indices ";

%feature("docstring")  indexsort::show "void indexsort::show() const
";

%feature("docstring")  indexsort::sorted "std::vector<Type>
indexsort::sorted(const std::vector< Type > &vals) const

return array sorted using the order from the indexsort structure ";

%feature("docstring")  indexsort::sorted "ContainerType
indexsort::sorted(const ContainerType &vals) const

return array sorted using the order from the indexsort structure ";

%feature("docstring")  indexsort::issorted "bool
indexsort::issorted() const

Returns true of the data is sorted ascending. ";

%feature("docstring")  indexsort::issorteddescending "bool
indexsort::issorteddescending() const

Returns true of the data is sorted descending. ";


// File: classInfInt.xml
%feature("docstring") InfInt "C++ includes: InfInt.h ";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt() ";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(const char *c)
";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(const
std::string &s) ";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(int l) ";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(long l) ";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(unsigned int l)
";

%feature("docstring")  InfInt::InfInt "InfInt::InfInt(unsigned long
l) ";

%feature("docstring")  InfInt::intSqrt "InfInt InfInt::intSqrt()
const ";

%feature("docstring")  InfInt::digitAt "char InfInt::digitAt(size_t
i) const ";

%feature("docstring")  InfInt::numberOfDigits "size_t
InfInt::numberOfDigits() const ";

%feature("docstring")  InfInt::size "size_t InfInt::size() const ";

%feature("docstring")  InfInt::toString "std::string
InfInt::toString() const ";

%feature("docstring")  InfInt::toInt "int InfInt::toInt() const ";

%feature("docstring")  InfInt::toLong "long InfInt::toLong() const ";

%feature("docstring")  InfInt::toUnsignedInt "unsigned int
InfInt::toUnsignedInt() const ";


// File: classJcounter.xml
%feature("docstring") Jcounter "

object to hold counts of maximum J_k-values

C++ includes: evenodd.h ";

%feature("docstring")  Jcounter::Jcounter "Jcounter::Jcounter()

time needed for calculation ";

%feature("docstring")  Jcounter::Jcounter "Jcounter::Jcounter(int N,
int jj=5, int k=-1) ";

%feature("docstring")  Jcounter::validData "bool
Jcounter::validData() ";

%feature("docstring")  Jcounter::hasColumn "bool
Jcounter::hasColumn(int col) const

return true if specified column is in the data ";

%feature("docstring")  Jcounter::isOpen "bool Jcounter::isOpen()
const ";

%feature("docstring")  Jcounter::showPerformance "void
Jcounter::showPerformance() const ";

%feature("docstring")  Jcounter::narrays "long Jcounter::narrays()
const ";

%feature("docstring")  Jcounter::show "void Jcounter::show() const

show statistics of the object ";

%feature("docstring")  Jcounter::maxCols "int Jcounter::maxCols()
const ";

%feature("docstring")  Jcounter::getCount "long
Jcounter::getCount(int k, int j) const ";

%feature("docstring")  Jcounter::getTotalsJvalue "std::vector<long>
Jcounter::getTotalsJvalue(int jval) const ";

%feature("docstring")  Jcounter::getTotals "std::vector<long>
Jcounter::getTotals() const ";

%feature("docstring")  Jcounter::showcompact "void
Jcounter::showcompact() const

show statistics of the object ";

%feature("docstring")  Jcounter::addArrays "void
Jcounter::addArrays(const arraylist_t &arraylist, int verbose=0)

add list of arrays to object ";

%feature("docstring")  Jcounter::addArray "void
Jcounter::addArray(const array_link &al, int verbose=0)

add single array to statistics object ";


// File: structjindex__t.xml
%feature("docstring") jindex_t "

helper class for indexing statistics of designs

The index consists of the number of columns and the value for the
J-characteristic

C++ includes: evenodd.h ";

%feature("docstring")  jindex_t::jindex_t "jindex_t::jindex_t(int
colindex, int jvalue) ";

%feature("docstring")  jindex_t::toString "std::string
jindex_t::toString() const ";


// File: classjstruct__t.xml
%feature("docstring") jstruct_t "

struct to hold data of an array, e.g. J-characteristic, rank

See papers: Minimum G2-aberration properties of two-level foldover
designs, Butler, 2004 Design Selection and Classification for Hadamard
Matrices Using Generalized Minimum Aberration Criteria, Deng and Tang

C++ includes: arraytools.h ";

%feature("docstring")  jstruct_t::jstruct_t "jstruct_t::jstruct_t()
";

%feature("docstring")  jstruct_t::jstruct_t "jstruct_t::jstruct_t(const int N, const int K, const int jj=4) ";

%feature("docstring")  jstruct_t::jstruct_t "jstruct_t::jstruct_t(const jstruct_t &js) ";

%feature("docstring")  jstruct_t::jstruct_t "jstruct_t::jstruct_t(const array_link &al, int jj=4) ";

%feature("docstring")  jstruct_t::~jstruct_t "jstruct_t::~jstruct_t()
";

%feature("docstring")  jstruct_t::maxJ "int jstruct_t::maxJ() const

calculate maximum J value ";

%feature("docstring")  jstruct_t::Fval "std::vector< int >
jstruct_t::Fval(int strength=3) const

calculate possible values in F vector ";

%feature("docstring")  jstruct_t::calculateF "std::vector< int >
jstruct_t::calculateF(int strength=3) const

calculate histogram of J values for a 2-level array ";

%feature("docstring")  jstruct_t::calculateAberration "void
jstruct_t::calculateAberration() ";

%feature("docstring")  jstruct_t::show "void jstruct_t::show()

Show contents of structure. ";

%feature("docstring")  jstruct_t::showdata "void
jstruct_t::showdata() ";

%feature("docstring")  jstruct_t::showstr "std::string
jstruct_t::showstr() ";

%feature("docstring")  jstruct_t::allzero "int jstruct_t::allzero()

return 1 if all vals are zero ";


// File: classjstructbase__t.xml
%feature("docstring") jstructbase_t "

struct to hold data of an array, e.g. J-characteristic. Abstract base
class

C++ includes: arraytools.h ";

%feature("docstring")  jstructbase_t::maxJ "int jstructbase_t::maxJ()
const

calculate maximum J value ";

%feature("docstring")  jstructbase_t::Jvalues "std::vector< int >
jstructbase_t::Jvalues() const

calculate possible values in F vector ";

%feature("docstring")  jstructbase_t::calculateF "std::vector< int >
jstructbase_t::calculateF() const

calculate histogram of J values ";

%feature("docstring")  jstructbase_t::calc "virtual void
jstructbase_t::calc(const array_link &al)=0 ";

%feature("docstring")  jstructbase_t::show "void
jstructbase_t::show()

Show contents of structure. ";

%feature("docstring")  jstructbase_t::showdata "void
jstructbase_t::showdata(int verbose=1) ";

%feature("docstring")  jstructbase_t::showstr "std::string
jstructbase_t::showstr() ";

%feature("docstring")  jstructbase_t::allzero "int
jstructbase_t::allzero()

return 1 if all vals are zero ";


// File: classjstructconference__t.xml
%feature("docstring") jstructconference_t "C++ includes: arraytools.h
";

%feature("docstring")  jstructconference_t::jstructconference_t "jstructconference_t::jstructconference_t(int N, int jj=4) ";

%feature("docstring")  jstructconference_t::jstructconference_t "jstructconference_t::jstructconference_t(const array_link &al, int
jj=4) ";


// File: classlarray.xml
%feature("docstring") larray "

lightweight array class

C++ includes: mathtools.h ";

%feature("docstring")  larray::larray "larray< numtype >::larray() ";

%feature("docstring")  larray::larray "larray< numtype
>::larray(const numtype *data, int nn) ";

%feature("docstring")  larray::larray "larray< numtype >::larray(int
nn) ";

%feature("docstring")  larray::begin "numtype* larray< numtype
>::begin() ";

%feature("docstring")  larray::larray "larray< numtype
>::larray(const larray &rhs) ";

%feature("docstring")  larray::~larray "larray< numtype >::~larray()
";

%feature("docstring")  larray::resize "void larray< numtype
>::resize(size_t n) ";

%feature("docstring")  larray::size "size_t larray< numtype >::size()
const ";

%feature("docstring")  larray::at "numtype& larray< numtype
>::at(size_t i) ";

%feature("docstring")  larray::addelement "larray larray< numtype
>::addelement(numtype v) const ";


// File: structLMC__static__struct__t.xml
%feature("docstring") LMC_static_struct_t "

Contains initialization data for static allocations.

Part of the allocations is for structures that are constant and are
re-used each time an LMC calculation is performed. Some other
structures are temporary buffers that are written to all the time.

C++ includes: lmc.h ";

%feature("docstring")  LMC_static_struct_t::LMC_static_struct_t "LMC_static_struct_t::LMC_static_struct_t() ";

%feature("docstring")  LMC_static_struct_t::~LMC_static_struct_t "LMC_static_struct_t::~LMC_static_struct_t() ";

%feature("docstring")  LMC_static_struct_t::show "void
LMC_static_struct_t::show(int verbose=1) const ";

%feature("docstring")  LMC_static_struct_t::init "void
LMC_static_struct_t::init(const arraydata_t *adp) ";

%feature("docstring")  LMC_static_struct_t::freeall "void
LMC_static_struct_t::freeall() ";

%feature("docstring")  LMC_static_struct_t::update "int
LMC_static_struct_t::update(const arraydata_t *adp)

update structure with new design specification ";

%feature("docstring")  LMC_static_struct_t::needUpdate "int
LMC_static_struct_t::needUpdate(const arraydata_t *adp) const ";

%feature("docstring")  LMC_static_struct_t::init_root_stage "void
LMC_static_struct_t::init_root_stage(levelperm_t *&lperm_p, colperm_t
*&colperm_p, const arraydata_t *adp) ";

%feature("docstring")  LMC_static_struct_t::init_nonroot_stage "void
LMC_static_struct_t::init_nonroot_stage(levelperm_t *&lperm_p,
colperm_t *&colperm_p, colperm_t *&localcolperm_p, dyndata_t
**&dynd_p, int &dynd_p_nelem, array_t *&colbuffer, const arraydata_t
*adp) const ";

%feature("docstring")  LMC_static_struct_t::init_rootrowperms "void
LMC_static_struct_t::init_rootrowperms(int &totalperms, rowperm_t
*&rootrowperms, levelperm_t *&lperm_p)

Static initialization of root row permutations. ";

%feature("docstring")  LMC_static_struct_t::init_rootrowperms_full "void LMC_static_struct_t::init_rootrowperms_full(int &totalperms,
rowperm_t *&rootrowperms, levelperm_t *&lperm_p)

Static initialization of root row permutations (full group) ";


// File: structLMCreduction__t.xml
%feature("docstring") LMCreduction_t "

Class to describe an LMC reduction.

The most important variable is the transformation itself, contained in
transformation. The state contains information about how the reduction
was performed.

C++ includes: lmc.h ";

%feature("docstring")  LMCreduction_t::LMCreduction_t "LMCreduction_t::LMCreduction_t(const LMCreduction_t &at) ";

%feature("docstring")  LMCreduction_t::LMCreduction_t "LMCreduction_t::LMCreduction_t(const arraydata_t *ad)

copy constructor ";

%feature("docstring")  LMCreduction_t::~LMCreduction_t "LMCreduction_t::~LMCreduction_t() ";

%feature("docstring")  LMCreduction_t::getArray "array_link
LMCreduction_t::getArray() const

Assignment operator. ";

%feature("docstring")  LMCreduction_t::setArray "void
LMCreduction_t::setArray(const array_link al) ";

%feature("docstring")  LMCreduction_t::setArray "void
LMCreduction_t::setArray(const array_t *array, int nrows, int ncols)
";

%feature("docstring")  LMCreduction_t::updateSDpointer "void
LMCreduction_t::updateSDpointer(const array_link al, bool cache=false)
";

%feature("docstring")  LMCreduction_t::clearSymmetries "void
LMCreduction_t::clearSymmetries() ";

%feature("docstring")  LMCreduction_t::releaseStatic "void
LMCreduction_t::releaseStatic() ";

%feature("docstring")  LMCreduction_t::initStatic "void
LMCreduction_t::initStatic()

acquire a reference to a LMC_static_struct_t object ";

%feature("docstring")  LMCreduction_t::getStaticReference "LMC_static_struct_t& LMCreduction_t::getStaticReference()

return a reference to a LMC_static_struct_t object ";

%feature("docstring")  LMCreduction_t::reset "void
LMCreduction_t::reset()

reset the reduction: clears the symmetries and sets the transformation
to zero ";

%feature("docstring")  LMCreduction_t::show "void
LMCreduction_t::show(int verbose=2) const ";

%feature("docstring")  LMCreduction_t::__repr__ "std::string
LMCreduction_t::__repr__() const ";

%feature("docstring")  LMCreduction_t::updateFromLoop "void
LMCreduction_t::updateFromLoop(const arraydata_t &ad, const dyndata_t
&dynd, levelperm_t *lperms, const array_t *original)

called whenever we find a reduction ";

%feature("docstring")  LMCreduction_t::updateTransformation "void
LMCreduction_t::updateTransformation(const arraydata_t &ad, const
dyndata_t &dynd, levelperm_t *lperms, const array_t *original) ";

%feature("docstring")  LMCreduction_t::doBreak "bool
LMCreduction_t::doBreak(lmc_t ret) ";

%feature("docstring")  LMCreduction_t::updateLastCol "void
LMCreduction_t::updateLastCol(int col) ";


// File: structmvalue__t.xml
%feature("docstring") mvalue_t "

Multi-value type.

This object represents a multi-valued object. The objects are ordered
using lexicographic ordering.

C++ includes: mathtools.h ";

%feature("docstring")  mvalue_t::mvalue_t "mvalue_t< NumericType
>::mvalue_t() ";

%feature("docstring")  mvalue_t::~mvalue_t "mvalue_t< NumericType
>::~mvalue_t() ";

%feature("docstring")  mvalue_t::mvalue_t "mvalue_t< NumericType
>::mvalue_t(NumericType m, direction_t dd=HIGH) ";

%feature("docstring")  mvalue_t::mvalue_t "mvalue_t< NumericType
>::mvalue_t(std::vector< NumericType > vv, direction_t dd=HIGH) ";

%feature("docstring")  mvalue_t::mvalue_t "mvalue_t< NumericType
>::mvalue_t(std::vector< T > vv, direction_t dd=HIGH) ";

%feature("docstring")  mvalue_t::size "size_t mvalue_t< NumericType
>::size() const ";

%feature("docstring")  mvalue_t::show_integer "void mvalue_t<
NumericType >::show_integer() const ";


// File: classOAextend.xml
%feature("docstring") OAextend "

Options for the extend code.

class containing parameters of the extension and LMC algorithm

C++ includes: extend.h ";

%feature("docstring")  OAextend::OAextend "OAextend::OAextend() ";

%feature("docstring")  OAextend::OAextend "OAextend::OAextend(const
OAextend &o) ";

%feature("docstring")  OAextend::OAextend "OAextend::OAextend(arraydata_t &ad) ";

%feature("docstring")  OAextend::setAlgorithm "void
OAextend::setAlgorithm(algorithm_t algorithm, arraydata_t *ad=0)

Set the algorithm to use for LMC checks. ";

%feature("docstring")  OAextend::setAlgorithmAuto "void
OAextend::setAlgorithmAuto(arraydata_t *ad=0)

Set the algorithm automatically. ";

%feature("docstring")  OAextend::getAlgorithm "algorithm_t
OAextend::getAlgorithm() const

Return algorithm used. ";

%feature("docstring")  OAextend::getAlgorithmName "std::string
OAextend::getAlgorithmName() const

Return algorithm used (as string) ";

%feature("docstring")  OAextend::updateArraydata "void
OAextend::updateArraydata(arraydata_t *ad=0) const ";

%feature("docstring")  OAextend::info "void OAextend::info(int vb=1)
const

print configuration to stdout ";

%feature("docstring")  OAextend::__repr__ "std::string
OAextend::__repr__() const ";


// File: classobject__pool.xml
%feature("docstring") object_pool "

Class to make a pool of objects that can be re-used.

C++ includes: mathtools.h ";

%feature("docstring")  object_pool::object_pool "object_pool< TYPE
>::object_pool()

constructor ";

%feature("docstring")  object_pool::reset "void object_pool< TYPE
>::reset()

assume the pool is filled with pointers: remove them all and call the
destructor ";

%feature("docstring")  object_pool::begin "iterator object_pool< TYPE
>::begin() ";

%feature("docstring")  object_pool::end "iterator object_pool< TYPE
>::end() ";

%feature("docstring")  object_pool::begin "const_iterator
object_pool< TYPE >::begin() const ";

%feature("docstring")  object_pool::end "const_iterator object_pool<
TYPE >::end() const ";

%feature("docstring")  object_pool::at "TYPE& object_pool< TYPE
>::at(std::size_t index) ";

%feature("docstring")  object_pool::at "TYPE const& object_pool< TYPE
>::at(std::size_t index) const ";

%feature("docstring")  object_pool::size "size_t object_pool< TYPE
>::size() const ";

%feature("docstring")  object_pool::New "TYPE* object_pool< TYPE
>::New() ";

%feature("docstring")  object_pool::Delete "void object_pool< TYPE
>::Delete(TYPE *t) ";


// File: classPareto.xml
%feature("docstring") Pareto "

Class to the calculate Pareto optimal elements.

The class is templated by the type of values to be compared and an
index type. The index type is used to index the elements.

C++ includes: pareto.h ";

%feature("docstring")  Pareto::Pareto "Pareto< ValueType, IndexType
>::Pareto()

constructor ";

%feature("docstring")  Pareto::~Pareto "Pareto< ValueType, IndexType
>::~Pareto() ";

%feature("docstring")  Pareto::number "int Pareto< ValueType,
IndexType >::number() const

return the total number of Pareto optimal values ";

%feature("docstring")  Pareto::numberindices "int Pareto< ValueType,
IndexType >::numberindices() const

return the toal number Pareto optimal objects ";

%feature("docstring")  Pareto::__repr__ "std::string Pareto<
ValueType, IndexType >::__repr__() const ";

%feature("docstring")  Pareto::show "void Pareto< ValueType,
IndexType >::show(int verbose=1)

show the current set of Pareto optimal elements ";

%feature("docstring")  Pareto::allindicesdeque "std::deque<IndexType>
Pareto< ValueType, IndexType >::allindicesdeque() const

return all indices of the Pareto optimal elements as a std::deque ";

%feature("docstring")  Pareto::allindices "std::vector<IndexType>
Pareto< ValueType, IndexType >::allindices() const

return all indices of the Pareto optimal elements ";

%feature("docstring")  Pareto::allvalues "std::vector<pValue> Pareto<
ValueType, IndexType >::allvalues() const

return all Paretop optimal elements ";

%feature("docstring")  Pareto::addvalue "bool Pareto< ValueType,
IndexType >::addvalue(const pValue val, const IndexType idx)

add a new element ";


// File: structpareto__element.xml
%feature("docstring") pareto_element "

helper class for the Pareto class

C++ includes: pareto.h ";

%feature("docstring")  pareto_element::dominates "bool
pareto_element< ValueType, IndexType >::dominates(pValue v)

return true of the argument element dominates this value ";

%feature("docstring")  pareto_element::isdominated "bool
pareto_element< ValueType, IndexType >::isdominated(pValue v) ";

%feature("docstring")  pareto_element::equal "bool pareto_element<
ValueType, IndexType >::equal(pValue v)

return true of the argument element is equal to this element ";


// File: classrankStructure.xml
%feature("docstring") rankStructure "

Structure to efficiently calculate the rank of the second order
interaction matrix of many arrays sharing a common subarray

The input arrays are assumed to be of the form A_i = [A_0 X_i]

C++ includes: arrayproperties.h ";

%feature("docstring")  rankStructure::rankStructure "rankStructure::rankStructure(const array_link &al, int nsub=3, int
verbose=0)

constructor ";

%feature("docstring")  rankStructure::rankStructure "rankStructure::rankStructure(int nsub=3, int id=-1)

constructor ";

%feature("docstring")  rankStructure::info "void
rankStructure::info() const ";

%feature("docstring")  rankStructure::updateStructure "void
rankStructure::updateStructure(const array_link &al)

update the structure cache with a new array ";

%feature("docstring")  rankStructure::matrixP "EigenDecomp::PermutationType rankStructure::matrixP() const

helper function ";

%feature("docstring")  rankStructure::rankdirect "int
rankStructure::rankdirect(const Eigen::MatrixXd &A) const

calculate the rank of an array directly, uses special threshold ";

%feature("docstring")  rankStructure::rankxfdirect "int
rankStructure::rankxfdirect(const array_link &al) const

calculate the rank of the second order interaction matrix of an array
directly ";

%feature("docstring")  rankStructure::rankxf "int
rankStructure::rankxf(const array_link &al)

calculate the rank of the second order interaction matrix of an array
using the cache system ";


// File: structgfx_1_1TimSort_1_1run.xml


// File: classsort__indices.xml
%feature("docstring") sort_indices "

Helper class.

C++ includes: mathtools.h ";

%feature("docstring")  sort_indices::sort_indices "sort_indices< Type
>::sort_indices(Type *parr) ";


// File: classsort__indices__container.xml
%feature("docstring") sort_indices_container "

Helper class.

C++ includes: mathtools.h ";

%feature("docstring")  sort_indices_container::sort_indices_container
"sort_indices_container< ContainerType
>::sort_indices_container(const ContainerType *parr, bool
ascendingx=true) ";


// File: classsort__indices__deque.xml
%feature("docstring") sort_indices_deque "

Helper class.

C++ includes: mathtools.h ";

%feature("docstring")  sort_indices_deque::sort_indices_deque "sort_indices_deque< Type >::sort_indices_deque(const std::deque< Type
> *parr, bool ascendingx=true) ";


// File: classsort__indices__vector.xml
%feature("docstring") sort_indices_vector "

Helper class.

C++ includes: mathtools.h ";

%feature("docstring")  sort_indices_vector::sort_indices_vector "sort_indices_vector< Type >::sort_indices_vector(const std::vector<
Type > *parr, bool ascendingx=true) ";


// File: classLMCreduction__t_1_1symm__t.xml
%feature("docstring") LMCreduction_t::symm_t "

store column permutations from array symmetry group

C++ includes: lmc.h ";

%feature("docstring")  LMCreduction_t::symm_t::symm_t "LMCreduction_t::symm_t::symm_t() ";

%feature("docstring")  LMCreduction_t::symm_t::show "void
LMCreduction_t::symm_t::show(int verbose=1) const ";

%feature("docstring")  LMCreduction_t::symm_t::valid "bool
LMCreduction_t::symm_t::valid() const ";

%feature("docstring")  LMCreduction_t::symm_t::makeColpermsUnique "void LMCreduction_t::symm_t::makeColpermsUnique(int dverbose=0) ";

%feature("docstring")  LMCreduction_t::symm_t::storeColumnCombination
"void LMCreduction_t::symm_t::storeColumnCombination(colpermtype cpv)
";

%feature("docstring")  LMCreduction_t::symm_t::storeColumnCombination
"void LMCreduction_t::symm_t::storeColumnCombination(const colperm_t
cp, int n) ";

%feature("docstring")  LMCreduction_t::symm_t::storeColumnPermutation
"void LMCreduction_t::symm_t::storeColumnPermutation(const colperm_t
cp, int n) ";

%feature("docstring")  LMCreduction_t::symm_t::showColperms "void
LMCreduction_t::symm_t::showColperms(int verbose=1) const ";

%feature("docstring")  LMCreduction_t::symm_t::showColcombs "void
LMCreduction_t::symm_t::showColcombs(int verbose=1) const ";

%feature("docstring")  LMCreduction_t::symm_t::showSymmetries "void
LMCreduction_t::symm_t::showSymmetries(int verbose=1) const ";

%feature("docstring")
LMCreduction_t::symm_t::storeSymmetryPermutation "void
LMCreduction_t::symm_t::storeSymmetryPermutation(const dyndata_t
*dyndata) ";


// File: structsymmdata.xml
%feature("docstring") symmdata "

structure containing data related to symmetries of arrays

C++ includes: arraytools.h ";

%feature("docstring")  symmdata::symmdata "symmdata::symmdata(const
array_link &al, int minlen=1) ";

%feature("docstring")  symmdata::show "void symmdata::show(int
verbose=1) const ";

%feature("docstring")  symmdata::checkIdx "std::vector<int>
symmdata::checkIdx(int col=-1) const

list with indices set to check for symmetry reductions ";


// File: classsymmetry__group.xml
%feature("docstring") symmetry_group "

Class to describe the symmetry group of a list of elements.

The class assumes the list is sorted. The symmetry group is then a
direct product of full permutation groups.

We do not implement this using templates because we want to export to
Python.

C++ includes: mathtools.h ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< int > &vals, bool
ascending=true, int verbose=0)

ordering of elements ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< double > &vals, bool
ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< float > &vals, bool
ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< short int > &vals,
bool ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< unsigned int >
&vals, bool ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< mvalue_t< double > >
&vals, bool ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const std::vector< mvalue_t< int > >
&vals, bool ascending=true, int verbose=0) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group(const symmetry_group &sgx) ";

%feature("docstring")  symmetry_group::symmetry_group "symmetry_group::symmetry_group()

copy constructor ";

%feature("docstring")  symmetry_group::permsize "perm_return_type
symmetry_group::permsize() const

return size of the group of all permutations respecting the symmetry
";

%feature("docstring")  symmetry_group::permsize_large "InfInt
symmetry_group::permsize_large() const

return size of the group of all permutations respecting the symmetry
";

%feature("docstring")  symmetry_group::checkIndices "std::vector<int>
symmetry_group::checkIndices() const

list with indices set to check for symmetry reductions ";

%feature("docstring")  symmetry_group::__repr__ "std::string
symmetry_group::__repr__() const

representation function (for python interface) ";

%feature("docstring")  symmetry_group::show "void
symmetry_group::show(int verbose=1) const

show the symmetry group ";


// File: classsymmetry__group__walker.xml
%feature("docstring") symmetry_group_walker "

class to walk over all element of a symmetry group

C++ includes: mathtools.h ";

%feature("docstring")  symmetry_group_walker::symmetry_group_walker "symmetry_group_walker::symmetry_group_walker(symmetry_group sg) ";

%feature("docstring")  symmetry_group_walker::show "void
symmetry_group_walker::show(int verbose=1) const ";

%feature("docstring")  symmetry_group_walker::init "void
symmetry_group_walker::init() ";

%feature("docstring")  symmetry_group_walker::next "bool
symmetry_group_walker::next() ";

%feature("docstring")  symmetry_group_walker::fullperm "std::vector<int> symmetry_group_walker::fullperm() const ";


// File: classgfx_1_1TimSort.xml
%feature("docstring") gfx::TimSort "C++ includes: timsort.hpp ";


// File: namespacearrayfile.xml


// File: namespaceEigen.xml


// File: namespacegfx.xml
%feature("docstring")  gfx::timsort "void
gfx::timsort(RandomAccessIterator const first, RandomAccessIterator
const last)

Same as std::stable_sort(first, last). ";

%feature("docstring")  gfx::timsort "void
gfx::timsort(RandomAccessIterator const first, RandomAccessIterator
const last, LessFunction compare)

Same as std::stable_sort(first, last, c). ";


// File: namespacenauty.xml
%feature("docstring")  nauty::reduceNauty "std::vector<int>
nauty::reduceNauty(const array_link &G, std::vector< int > colors, int
verbose=0)

reduce a colored graph to Nauty minimal form

The transformation returned is from the normal form to the specified
graph. ";


// File: namespacestd.xml


// File: analysis_8h.xml


// File: anyoption_8h.xml


// File: arrayproperties_8h.xml
%feature("docstring")  detXtX "double detXtX(const Eigen::MatrixXd
&mymatrix, int verbose=1)

calculate determinant of X^T X by using the SVD ";

%feature("docstring")  detXtXfloat "double detXtXfloat(const
Eigen::MatrixXf &mymatrix, int verbose=1) ";

%feature("docstring")  DAEefficiecyWithSVD "void
DAEefficiecyWithSVD(const Eigen::MatrixXd &x, double &Deff, double
&vif, double &Eeff, int &rank, int verbose)

Calculate D-efficiency and VIF-efficiency and E-efficiency values
using SVD. ";

%feature("docstring")  array_rank_D_B "int array_rank_D_B(const
array_link &al, std::vector< double > *ret=0, int verbose=0)

Calculate the rank of the second order interaction matrix of an
orthogonal array, the rank, D-efficiency, VIF-efficiency and
E-efficiency are appended to the second argument. ";

%feature("docstring")  Defficiency "double Defficiency(const
array_link &al, int verbose=0)

Calculate D-efficiency for a 2-level array using symmetric eigenvalue
decomposition. ";

%feature("docstring")  Defficiencies "std::vector< double >
Defficiencies(const array_link &al, const arraydata_t &arrayclass, int
verbose=0, int addDs0=0) ";

%feature("docstring")  VIFefficiency "double VIFefficiency(const
array_link &al, int verbose=0)

Calculate VIF-efficiency of matrix. ";

%feature("docstring")  Aefficiency "double Aefficiency(const
array_link &al, int verbose=0)

Calculate A-efficiency of matrix. ";

%feature("docstring")  Eefficiency "double Eefficiency(const
array_link &al, int verbose=0)

Calculate E-efficiency of matrix (1 over the VIF-efficiency) ";

%feature("docstring")  Aefficiencies "std::vector< double >
Aefficiencies(const array_link &al, int verbose=0)

calculate various A-efficiencies ";

%feature("docstring")  projDeff "std::vector< double > projDeff(const
array_link &al, int kp, int verbose)

Return the D-efficiencies for the projection designs. ";

%feature("docstring")  PECsequence "std::vector< double >
PECsequence(const array_link &al, int verbose=0)

Return the projection estimation capacity sequence of a design. ";

%feature("docstring")  distance_distribution "std::vector< double >
distance_distribution(const array_link &al)

Return the distance distribution of a design. ";

%feature("docstring")  Jcharacteristics "std::vector< int >
Jcharacteristics(const array_link &al, int jj=4, int verbose=0)

Calculate J-characteristics of matrix (the values are signed) ";

%feature("docstring")  GWLP "std::vector< double > GWLP(const
array_link &al, int verbose=0, int truncate=1)

calculate GWLP (generalized wordlength pattern)

The method used for calculation is from Xu and Wu (2001),
\"Generalized minimum aberration for asymmetrical fractional factorial
desings\" The non-symmetric arrays see \"Algorithmic Construction of
Efficient Fractional Factorial Designs With Large Run Sizes\", Xu ";

%feature("docstring")  GWLPmixed "std::vector< double >
GWLPmixed(const array_link &al, int verbose=0, int truncate=1) ";

%feature("docstring")  projectionGWLPs "std::vector< GWLPvalue >
projectionGWLPs(const array_link &al)

calculate delete-one-factor GWLP (generalized wordlength pattern)
projections ";

%feature("docstring")  sortGWLP "std::vector< GWLPvalue >
sortGWLP(std::vector< GWLPvalue >) ";

%feature("docstring")  projectionGWLPvalues "std::vector< double >
projectionGWLPvalues(const array_link &al)

calculate delete-one-factor GWLP (generalized wordlength pattern)
projection values ";

%feature("docstring")  CL2discrepancy "double CL2discrepancy(const
array_link &al)

calculate centered L2-discrepancy

The method is from \"A connection between uniformity and aberration in
regular fractions of two-level factorials\", Fang and Mukerjee, 2000
";

%feature("docstring")  array2secondorder "array_link
array2secondorder(const array_link &al)

calculate second order interaction matrix for 2-level array (does not
include intercept and matrix itself) ";

%feature("docstring")  array2xf "array_link array2xf(const array_link
&al)

add intercept and second order interactions to a 2-level array ";

%feature("docstring")  array2xfeigen "Eigen::MatrixXd
array2xfeigen(const array_link &al)

add intercept and second order interactions to an array ";

%feature("docstring")  arrayrankFullPivQR "int
arrayrankFullPivQR(const array_link &al, double threshold=-1)

return rank of an array based on FullPivHouseholderQR (defined at
compile time) ";

%feature("docstring")  arrayrankColPivQR "int arrayrankColPivQR(const
array_link &al, double threshold=-1)

return rank of an array based on ColPivHouseholderQR (defined at
compile time) ";

%feature("docstring")  arrayrankFullPivLU "int
arrayrankFullPivLU(const array_link &al, double threshold=-1) ";

%feature("docstring")  arrayrankSVD "int arrayrankSVD(const
array_link &al, double threshold=-1) ";

%feature("docstring")  arrayrank "int arrayrank(const array_link &al)

calculate the rank of an array ";

%feature("docstring")  arrayrankInfo "int arrayrankInfo(const
Eigen::MatrixXd &, int verbose=1)

return rank of an array based on Eigen::FullPivLU ";

%feature("docstring")  arrayrankInfo "int arrayrankInfo(const
array_link &al, int verbose=1)

print information related to rank calculations ";

%feature("docstring")  arraylink2eigen "Eigen::MatrixXd
arraylink2eigen(const array_link &al)

convert array_link to Eigen matrix ";

%feature("docstring")  conditionNumber "double conditionNumber(const
array_link &M)

return the condition number of a matrix ";

%feature("docstring")  calculateParetoEvenOdd "void
calculateParetoEvenOdd(const std::vector< std::string > infiles, const
char *outfile, int verbose=1, arrayfilemode_t afmode=ABINARY, int
nrows=-1, int ncols=-1, paretomethod_t
paretomethod=PARETOFUNCTION_DEFAULT)

Calculate the Pareto optimal arrays from a list of array files

Pareto optimality is calculated according to (rank; A3,A4; F4) ";

%feature("docstring")  parsePareto "Pareto< mvalue_t < long >, long >
parsePareto(const arraylist_t &arraylist, int verbose, paretomethod_t
paretomethod=PARETOFUNCTION_DEFAULT) ";

%feature("docstring")  A3A4 "mvalue_t< long > A3A4(const array_link
&al)

calculate A3 and A4 value for array ";

%feature("docstring")  F4 "mvalue_t< long > F4(const array_link &al,
int verbose=1)

calculate F4 value for 2-level array ";

%feature("docstring")  calculateArrayParetoRankFA "Pareto< mvalue_t <
long >, IndexType >::pValue calculateArrayParetoRankFA(const
array_link &al, int verbose)

Add array to list of Pareto optimal arrays

The values to be optimized are:

1) Rank (higher is better) 2) A3, A4 (lower is better) 3) F4 (?? is
better, sum of elements is constant)

Valid for 2-level arrays of strength at least 3 ";

%feature("docstring")  addJmax "void addJmax(const array_link &al,
typename Pareto< mvalue_t< long >, IndexType >::pValue &p, int
verbose=1)

add Jmax criterium to Pareto set ";

%feature("docstring")  calculateArrayParetoJ5 "Pareto< mvalue_t <
long >, IndexType >::pValue calculateArrayParetoJ5(const array_link
&al, int verbose) ";

%feature("docstring")  parseArrayPareto "void parseArrayPareto(const
array_link &al, IndexType i, Pareto< mvalue_t< long >, IndexType >
&pset, int verbose)

Add array to list of Pareto optimal arrays

The values to be optimized are:

1) Rank (higher is better) 2) A3, A4 (lower is better) 3) F4 (?? is
better, sum of elements is constant) ";

%feature("docstring")  Cvalue2Dvalue "double Cvalue2Dvalue(double C,
int ka)

convert C value to D-efficiency value ";

%feature("docstring")  Dvalue2Cvalue "double Dvalue2Cvalue(double A,
int ka)

convert D-efficiency value to C value ";

%feature("docstring")  get_oaindex "int get_oaindex(const array_t *s,
const colindex_t strength, const colindex_t N)

Return index of an array. ";


// File: arraytools_8h.xml
%feature("docstring")  Eigen::eigenInfo "void eigenInfo(const
MatrixFloat m, const char *str=\"eigen\", int verbose=1)

show information about Eigen matrix ";

%feature("docstring")  Eigen::eigen2numpyHelper "void
eigen2numpyHelper(double *pymat1, int n, const MatrixFloat &m) ";

%feature("docstring")  Eigen::dummy "Eigen::VectorXd dummy() ";

%feature("docstring")  Eigen::dummy2 "Eigen::MatrixXd dummy2() ";

%feature("docstring")  Eigen::Fval "std::vector< int > Fval(int N,
int strength)

possible values for J-values of 2-level design ";

%feature("docstring")  Eigen::file_exists "bool file_exists(const
std::string filename)

return true if the specified file exists ";

%feature("docstring")  Eigen::file_exists "bool file_exists(const
char *filename)

return true if the specified file exists ";

%feature("docstring")  Eigen::oa_file_exists "bool
oa_file_exists(const char *filename)

return true if the specified oa file exists ";

%feature("docstring")  Eigen::oa_file_exists "bool
oa_file_exists(const std::string filename)

return true if the specified oa file exists ";

%feature("docstring")  Eigen::createJ2tableConference "array_link
createJ2tableConference(const array_link &confmatrix)

create J2 table as intermediate result for J-characteristic
calculations for conference matrices ";

%feature("docstring")  Eigen::createJdtable "array_link
createJdtable(const array_link &al)

create J2 table as intermediate result for J-characteristic
calculations ";

%feature("docstring")  Eigen::readConfigFile "arraydata_t*
readConfigFile(const char *file)

Read array configuration from file. ";

%feature("docstring")  Eigen::printfstring "std::string
printfstring(const char *message,...)

Function similar to printf returning C++ style string.

printf-style function that returns std::string

Parameters:
-----------

message:  ";

%feature("docstring")  Eigen::copy_array "void copy_array(const
array_t *src, array_t *const dst, const int nrows, const int ncols)

Make a copy of an array. ";

%feature("docstring")  Eigen::destroy_array "int
destroy_array(array_t *array)

Delete an array.

Parameters:
-----------

array:  ";

%feature("docstring")  Eigen::create_array "static array_t*
create_array(const int nrows, const int ncols)

Create an array.

Parameters:
-----------

nrows:  Number of rows

ncols:  Number of columns ";

%feature("docstring")  Eigen::create_array "array_t*
create_array(const arraydata_t *ad)

Create an array from an arraydata_t structure. ";

%feature("docstring")  Eigen::equal_array_cols "int
equal_array_cols(carray_t *A, colindex_t col, colindex_t col2,
rowindex_t nrows, rowindex_t rstart, rowindex_t rend)

Compare 2 columns of an array.

Parameters:
-----------

A:

col:

col2:

nrows:

rstart:

rend:  ";

%feature("docstring")  Eigen::clone_array "array_t* clone_array(const
array_t *const array, const rowindex_t nrows, const colindex_t ncols)

Clone an array. ";

%feature("docstring")  Eigen::perform_inv_column_permutation "void
perform_inv_column_permutation(const array_t *source, array_t *target,
colperm_t perm, int nrows, int ncols)

Perform inverse column permutation on an array.

Parameters:
-----------

source:

target:

perm:

nrows:

ncols:  ";

%feature("docstring")  Eigen::perform_column_permutation "void
perform_column_permutation(carray_t *source, array_t *target,
colperm_t perm, int nrows, int ncols) ";

%feature("docstring")  Eigen::perform_row_permutation "void
perform_row_permutation(const array_t *source, array_t *target,
rowperm_t perm, int nrows, int ncols)

Perform a row permutation.

Parameters:
-----------

source:  Source array

target:  Target array

perm:  Permutation to perform

nrows:  Number of rows

ncols:  Numer of columns ";

%feature("docstring")  Eigen::perform_inv_row_permutation "void
perform_inv_row_permutation(const array_t *source, array_t *target,
rowperm_t perm, int nrows, int ncols)

apply inverse row permutation ";

%feature("docstring")  Eigen::exampleArray "array_link
exampleArray(int idx=0, int verbose=0)

Return example array ";

%feature("docstring")  Eigen::Jcharacteristics_conference "std::vector<int> Jcharacteristics_conference(const array_link &al, int
jj, int verbose=0)

calculate J-characteristics for a conference design ";

%feature("docstring")  Eigen::hstack "array_link hstack(const
array_link &al, const array_link &b) ";

%feature("docstring")  Eigen::hstack "array_link hstack(const
array_link &al, const cperm &b) ";

%feature("docstring")  Eigen::hstacklastcol "array_link
hstacklastcol(const array_link &A, const array_link &B) ";

%feature("docstring")  Eigen::vstack "cperm vstack(const cperm &A,
const cperm &B) ";

%feature("docstring")  Eigen::perform_column_permutation "void
perform_column_permutation(const array_link source, array_link
&target, const std::vector< int >perm)

perform column permutation for an array ";

%feature("docstring")  Eigen::perform_row_permutation "void
perform_row_permutation(const array_link source, array_link &target,
const std::vector< int >perm)

perform row permutation for an array ";

%feature("docstring")  Eigen::arraylink2arraydata "arraydata_t
arraylink2arraydata(const array_link &al, int extracols=0, int
strength=2)

create arraydata_t structure from array ";

%feature("docstring")  Eigen::addConstant "arraylist_t
addConstant(const arraylist_t &lst, int v)

add a constant value to all arrays in a list ";

%feature("docstring")  Eigen::getJcounts "std::vector< int >
getJcounts(arraylist_t *arraylist, int N, int k, int verbose=1)

Return number of arrays with j_{2n+1}=0 for n<m ";

%feature("docstring")  Eigen::predictJ "int predictJ(const array_t
*array, const int N, const int k)

Predict j4(1,2,3,k) using the theorem from Deng This works only for
2-level arrays. The 0 corresponds to a +. ";

%feature("docstring")  Eigen::create_root "void create_root(array_t
*array, const arraydata_t *ad)

set first columns of an array to root form ";

%feature("docstring")  Eigen::create_root "void create_root(const
arraydata_t *ad, arraylist_t &solutions)

Creates the root of an OA. The root is appended to the current list of
arrays. ";

%feature("docstring")  Eigen::array_diff "int array_diff(carray_p A,
carray_p B, const rowindex_t r, const colindex_t c, rowindex_t &rpos,
colindex_t &cpos)

Compare 2 arrays and return position of first difference. ";

%feature("docstring")  Eigen::array_cmp "int array_cmp(carray_p A,
carray_p B, const rowindex_t r, const colindex_t c)

Compare 2 arrays and return 0 if equal. ";

%feature("docstring")  Eigen::fastJupdateValue "int
fastJupdateValue(rowindex_t N, carray_t *tmpval)

helper function to calculate J-values ";

%feature("docstring")  Eigen::fastJupdate "void fastJupdate(const
array_t *array, rowindex_t N, const int J, const colindex_t *pp,
array_t *tmp)

helper function to calculate J-values ";

%feature("docstring")  Eigen::jvalue "int jvalue(const array_link
&ar, const int J, const int *pp)

Calculate J-value for a 2-level array ";

%feature("docstring")  Eigen::jvaluefast "int jvaluefast(const
array_t *array, rowindex_t N, const int J, const colindex_t *pp)

calculate J-value for a 2-level array ";

%feature("docstring")  Eigen::analyseArrays "std::vector< jstruct_t >
analyseArrays(const arraylist_t &arraylist, const int verbose, const
int jj=4)

Analyse a list of arrays. ";

%feature("docstring")  Eigen::showArrayList "void showArrayList(const
arraylist_t &lst)

functions for working with array files ";

%feature("docstring")  Eigen::nArrays "long nArrays(const char
*fname)

return number of arrays in an array file ";

%feature("docstring")  Eigen::arrayfileinfo "void arrayfileinfo(const
char *fname, int &n, int &nr, int &nc)

return number of arrays in an array file ";

%feature("docstring")  Eigen::readarrayfile "int readarrayfile(const
char *fname, arraylist_t *arraylist, int verbose=1, int *setcols=0,
rowindex_t *setrows=0, int *setbits=0)

read list of arrays from file and append to list ";

%feature("docstring")  Eigen::readarrayfile "arraylist_t
readarrayfile(const char *fname, int verbose=1, int *setcols=0)

read list of arrays from file ";

%feature("docstring")  Eigen::writearrayfile "int
writearrayfile(const char *fname, const arraylist_t *arraylist,
arrayfile::arrayfilemode_t mode=arrayfile::ATEXT, int nrows=NRAUTO,
int ncols=NRAUTO)

write a list of arrays to file on disk ";

%feature("docstring")  Eigen::writearrayfile "int
writearrayfile(const char *fname, const arraylist_t arraylist,
arrayfile::arrayfilemode_t mode=arrayfile::ATEXT, int nrows=NRAUTO,
int ncols=NRAUTO)

write a list of arrays to file on disk ";

%feature("docstring")  Eigen::writearrayfile "int
writearrayfile(const char *fname, const array_link &al,
arrayfile::arrayfilemode_t mode=arrayfile::ATEXT)

write a single array to file ";

%feature("docstring")  Eigen::appendarrayfile "int
appendarrayfile(const char *fname, const array_link al)

append a single array to an array file. creates a new file if no file
exists ";

%feature("docstring")  Eigen::selectArrays "void selectArrays(const
std::string filename, std::vector< int > &idx, arraylist_t &fl, int
verbose=0)

Make a selection of arrays from binary array file, append to list. ";

%feature("docstring")  Eigen::selectArrays "array_link
selectArrays(std::string filename, int ii)

Select a single array from a file. ";

%feature("docstring")  Eigen::create_arrayfile "arrayfile_t*
create_arrayfile(const char *fname, int rows, int cols, int narrays,
arrayfile::arrayfilemode_t mode=arrayfile::ATEXT, int nbits=8) ";

%feature("docstring")  Eigen::save_arrays "int
save_arrays(arraylist_t &solutions, const arraydata_t *ad, const int
n_arrays, const int n_procs, const char *resultprefix,
arrayfile::arrayfilemode_t mode=ATEXT) ";

%feature("docstring")  Eigen::write_array_format "void
write_array_format(std::ostream &ss, const atype *array, const int
nrows, const int ncols, int width=3)

write array to output stream ";

%feature("docstring")  Eigen::selectArrays "arraylist_t
selectArrays(const arraylist_t &al, std::vector< int > &idx)

Make a selection of arrays. ";

%feature("docstring")  Eigen::selectArrays "arraylist_t
selectArrays(const arraylist_t &al, std::vector< long > &idx)

Make a selection of arrays. ";

%feature("docstring")  Eigen::selectArrays "void selectArrays(const
arraylist_t &al, std::vector< int > &idx, arraylist_t &fl)

Make a selection of arrays, append to list. ";

%feature("docstring")  Eigen::selectArrays "void selectArrays(const
arraylist_t &al, std::vector< long > &idx, arraylist_t &fl) ";

%feature("docstring")  Eigen::keepElements "void
keepElements(Container &al, std::vector< IntType > &idx)

Make a selection of arrays, keep. ";

%feature("docstring")  Eigen::removeElements "void
removeElements(Container &al, std::vector< IntType > &idx)

Make a selection of arrays, remove. ";

%feature("docstring")  Eigen::selectArraysMask "void
selectArraysMask(const arraylist_t &al, std::vector< MType > &mask,
arraylist_t &rl)

Make a selection of arrays from a list, append to list. ";

%feature("docstring")  Eigen::appendArrays "void appendArrays(const
arraylist_t &al, const typename std::vector< IndexType > &idx,
arraylist_t &lst)

Append selection of arrays to existing list. ";

%feature("docstring")  Eigen::appendArrays "void appendArrays(const
arraylist_t &al, arraylist_t &dst)

Append selection of arrays to existing list. ";

%feature("docstring")  Eigen::write_array_format "void
write_array_format(const atype *array, const int nrows, const int
ncols, int width=3)

Write a formatted array ";

%feature("docstring")  Eigen::write_array_format "void
write_array_format(FILE *fid, const atype *array, const int nrows,
const int ncols)

Write a formatted array. ";

%feature("docstring")  Eigen::write_array_latex "void
write_array_latex(std::ostream &ss, const atype *array, const int
nrows, const int ncols)

write an array in latex style ";

%feature("docstring")  Eigen::readbinheader "bool readbinheader(FILE
*fid, int &nr, int &nc)

Read header for binary data file. Return true if valid header file. ";

%feature("docstring")  Eigen::writebinheader "void
writebinheader(FILE *fid, int nr, int nc)

Write header for binary data file. ";

%feature("docstring")  Eigen::doublevector2binfile "void
doublevector2binfile(const std::string fname, std::vector< Type >
vals, int writeheader=1)

Write a vector of integer elements to file. ";

%feature("docstring")  Eigen::vectorvector2binfile "void
vectorvector2binfile(const std::string fname, const std::vector<
std::vector< double > >vals, int writeheader, int na)

Write a vector of vector elements to binary file. ";

%feature("docstring")  Eigen::array2eigenX2 "MatrixFloat
array2eigenX2(const array_link &al)

convert 2-level array to second order interaction matrix in Eigen
format ";

%feature("docstring")  Eigen::array2eigenX1 "MatrixFloat
array2eigenX1(const array_link &al, int intercept=1) ";

%feature("docstring")  Eigen::array2eigenModelMatrix "MatrixFloat
array2eigenModelMatrix(const array_link &al)

convert 2-level array to second order model matrix (intercept, X1, X2)
";

%feature("docstring")  Eigen::array2eigenME "MatrixFloat
array2eigenME(const array_link &al, int verbose=1)

convert array to model matrix in Eigen format ";

%feature("docstring")  Eigen::array2eigenModelMatrixMixed "std::pair<
MatrixFloat, MatrixFloat > array2eigenModelMatrixMixed(const
array_link &al, int verbose=1)

create first and second order model matrix for mixed-level array ";

%feature("docstring")  Eigen::arrayInFile "int arrayInFile(const
array_link &al, const char *afile, int verbose=1)

return index of specified array in a file. returns -1 if array is not
found ";

%feature("docstring")  Eigen::arrayInList "int arrayInList(const
array_link &al, const arraylist_t &ll, int verbose=1)

return index of specified array in a list. returns -1 if array is not
found ";


// File: conference_8h.xml
%feature("docstring")  print_cperm "void print_cperm(const cperm &c)

print a candidate extension ";

%feature("docstring")  partial_inner_product "int
partial_inner_product(const cperm &a, const array_link &al, int col,
int rmax)

partial inner product ";

%feature("docstring")  showCandidates "void showCandidates(const
std::vector< cperm > &cc)

show a list of candidate extensions ";

%feature("docstring")  reduceConference "array_link
reduceConference(const array_link &, int verbose=0)

reduce conference matrix to normal form ";

%feature("docstring")  reduceConferenceTransformation "conference_transformation_t reduceConferenceTransformation(const
array_link &al, int verbose)

reduce conference matrix to normal form ";

%feature("docstring")  extend_conference_matrix "conference_extend_t
extend_conference_matrix(const array_link &al, const conference_t &ct,
int extcol, int verbose=1, int maxzpos=-1)

Extend a single conference design with candidate columns ";

%feature("docstring")  extend_conference_matrix_generator "conference_extend_t extend_conference_matrix_generator(const
array_link &al, const conference_t &ct, int extcol, int verbose, int
maxzpos, const CandidateGenerator &cgenerator)

helper function ";

%feature("docstring")  extend_conference "arraylist_t
extend_conference(const arraylist_t &lst, const conference_t ctype,
int verbose, int select_isomorphism_classes=0)

Extend a list of conference designs with a single column. ";

%feature("docstring")  extend_conference_plain "arraylist_t
extend_conference_plain(const arraylist_t &lst, const conference_t
ctype, int verbose, int select_isomorphism_classes=0)

plain version without caching ";

%feature("docstring")  extend_conference_restricted "arraylist_t
extend_conference_restricted(const arraylist_t &lst, const
conference_t ctype, int verbose)

Extend a list of conference designs with a single column ";

%feature("docstring")  extend_double_conference "arraylist_t
extend_double_conference(const arraylist_t &lst, const conference_t
ctype, int verbose) ";

%feature("docstring")  selectConferenceIsomorpismClasses "arraylist_t
selectConferenceIsomorpismClasses(const arraylist_t &list, int
verbose, matrix_isomorphism_t itype=CONFERENCE_ISOMORPHISM)

select representatives for the isomorphism classes of a list of
conference arrays ";

%feature("docstring")  selectConferenceIsomorpismIndices "std::vector<int> selectConferenceIsomorpismIndices(const arraylist_t
&lst, int verbose, matrix_isomorphism_t itype=CONFERENCE_ISOMORPHISM)

select representatives for the isomorphism classes of a list of
conference arrays, return indices of classes ";

%feature("docstring")  selectLMC0 "arraylist_t selectLMC0(const
arraylist_t &list, int verbose, const conference_t &ctype) ";

%feature("docstring")  generateConferenceExtensions "std::vector<cperm> generateConferenceExtensions(const array_link &al,
const conference_t &ct, int kz, int verbose=1, int filtersymm=1, int
filterip=1)

Generate candidate extensions

Parameters:
-----------

al:  design to be extended

kz:  index of zero in candidate column ";

%feature("docstring")  generateConferenceRestrictedExtensions "std::vector<cperm> generateConferenceRestrictedExtensions(const
array_link &al, const conference_t &ct, int kz, int verbose=1, int
filtersymm=1, int filterip=1)

Generate candidate extensions for restricted isomorphism classes ";

%feature("docstring")  generateDoubleConferenceExtensions "std::vector<cperm> generateDoubleConferenceExtensions(const array_link
&al, const conference_t &ct, int verbose=1, int filtersymm=1, int
filterip=1, int filterJ3=0, int filtersymminline=1) ";

%feature("docstring")  maxz "int maxz(const array_link &al, int k=-1)

return max position of zero in array, returns -1 if no zero is found

The parameter k specifies the column to search in. For k=-1 all
columns are searched. ";

%feature("docstring")  compareLMC0 "bool compareLMC0(const array_link
&alL, const array_link &alR)

Return true of the array is smaller in LMC-0 ordering ";

%feature("docstring")  sortLMC0 "arraylist_t sortLMC0(const
arraylist_t &lst)

sort list of arrays according to LMC-0 ordering ";

%feature("docstring")  LMC0check "lmc_t LMC0check(const array_link
&al, int verbose=0)

check if array is in LM0 form ";

%feature("docstring")  isConferenceFoldover "bool
isConferenceFoldover(const array_link &al, int verbose=0)

return true if the design is a foldover array ";

%feature("docstring")  satisfy_symm "int satisfy_symm(const cperm &c,
const symmdata &sd, int rowstart=2)

helper function, return true if a candidate extensions satisfies the
symmetry test ";

%feature("docstring")  satisfy_symm "int satisfy_symm(const cperm &c,
const std::vector< int > &check_indices, int rowstart=2)

helper function, return true if a candidate extensions satisfies the
symmetry test ";

%feature("docstring")  satisfy_symm "int satisfy_symm(const cperm &c,
const std::vector< int > &check_indices, int rowstart, int rowend)

helper function, return true if a candidate extensions satisfies the
symmetry test ";

%feature("docstring")  ipcheck "int ipcheck(const cperm &col, const
array_link &al, int cstart=2, int verbose=0) ";

%feature("docstring")  minz "int minz(const array_link &al, int k)

return minimal position of zero in design ";

%feature("docstring")  generateDoubleConferenceExtensionsInflate "std::vector<cperm> generateDoubleConferenceExtensionsInflate(const
array_link &al, const conference_t &ct, int verbose, int filterj2, int
filterj3, int kstart=2) ";

%feature("docstring")  inflateCandidateExtension "std::vector<cperm>
inflateCandidateExtension(const cperm &basecandidate, const array_link
&als, const symmetry_group &alsg, const std::vector< int >
&check_indices, const conference_t &ct, int verbose, const
DconferenceFilter &filter) ";


// File: Deff_8h.xml
%feature("docstring")  scoreD "double scoreD(const std::vector<
double > dd, const std::vector< double > alpha)

calculate score from from set of efficiencies ";

%feature("docstring")  optimDeff "array_link optimDeff(const
array_link &A0, const arraydata_t &arrayclass, std::vector< double >
alpha, int verbose=1, int optimmethod=DOPTIM_AUTOMATIC, int
niter=100000, int nabort=0)

Optimize a design according to the optimization function specified.

Arguments: arrayclass: structure describing the design class alpha:
(3x1 array) verbose: output level ";

%feature("docstring")  optimDeff2level "array_link
optimDeff2level(const array_link &A0, const arraydata_t &arrayclass,
std::vector< double > alpha, int verbose=1, int
optimmethod=DOPTIM_AUTOMATIC, int niter=100000, int nabort=0)

debugging function ";

%feature("docstring")  Doptimize "DoptimReturn Doptimize(const
arraydata_t &arrayclass, int nrestarts, const std::vector< double >
alpha, int verbose, int method=DOPTIM_AUTOMATIC, int niter=300000,
double maxtime=100000, int nabort=5000)

function to generate optimal designs ";

%feature("docstring")  DoptimizeDebug "DoptimReturn
DoptimizeDebug(const arraydata_t &arrayclass, int nrestarts, const
std::vector< double > alpha, int verbose, int method=DOPTIM_AUTOMATIC,
int niter=300000, double maxtime=100000, int nabort=5000)

helper function ";

%feature("docstring")  DoptimizeMixed "DoptimReturn
DoptimizeMixed(const arraylist_t &sols, const arraydata_t &arrayclass,
const std::vector< double > alpha, int verbose=1, int nabort=-1) ";


// File: evenodd_8h.xml
%feature("docstring")  processDepth "void processDepth(const
arraylist_t &goodarrays, depth_alg_t depthalg, depth_extend_t
&dextend, depth_extend_sub_t &dextendsublight, int extensioncol, int
verbose=0)

Extend arrays using a depth-first or breadth-first approach

Parameters:
-----------

goodarrays:  List of arrays to extend

depthalg:  Extend using depth-first or breadth-first

dextend:  Option structure for the extension

dextendsublight:  Data structure for the extensions

extensioncol:  Column to extend

verbose:  Verbosity level ";

%feature("docstring")  depth_extend_hybrid "void
depth_extend_hybrid(const arraylist_t &alist, depth_extend_t &dextend,
int extcol, const OAextend &oaextendx, int verbose)

depth-first extension of arrays. depending on the symmetry group of
the array to be extended a direct method is used or a method with
caching of candidate columns ";

%feature("docstring")  depth_extend_direct "void
depth_extend_direct(const arraylist_t &alist, depth_extend_t &dextend,
int extcol, const OAextend &oaextendx, int verbose)

variation of depth_extend for arrays with large symmetry groups ";

%feature("docstring")  depth_extend_array "void
depth_extend_array(const array_link &al, depth_extend_t &dextend,
const arraydata_t &adfull, int verbose, depth_extensions_storage_t
*ds=0, int=0)

perform depth-first extension

The arrays generated are pruned by keeping a list of possible
extension valuesdepth extend a single array ";

%feature("docstring")  calculateArrayParetoJ5Cache "Pareto< mvalue_t
< long >, IndexType >::pValue calculateArrayParetoJ5Cache(const
array_link &al, int verbose, rankStructure &rs) ";

%feature("docstring")  addArraysToPareto "void
addArraysToPareto(Pareto< mvalue_t< long >, array_link > &pset,
pareto_cb paretofunction, const arraylist_t &arraylist, int jj, int
verbose)

add arrays to set of Pareto results ";

%feature("docstring")  addArraysToPareto "void
addArraysToPareto(Pareto< mvalue_t< long >, array_link > &pset,
pareto_cb_cache paretofunction, const arraylist_t &arraylist, int jj,
int verbose)

add arrays to set of Pareto results ";

%feature("docstring")  readStatisticsFile "Jcounter
readStatisticsFile(const char *numbersfile, int verbose)

read statistics object from disk ";

%feature("docstring")  writeStatisticsFile "void
writeStatisticsFile(const char *numbersfile, const Jcounter &jc, int
verbose)

write statistics object to disk ";

%feature("docstring")  calculateJstatistics "Jcounter
calculateJstatistics(const char *afile, int jj=5, int verbose=1)

calculate J-value statistics ";


// File: extend_8h.xml
%feature("docstring")  progress_column "double
progress_column(array_t *column, extendpos *p) ";

%feature("docstring")  extend_arraylist "int extend_arraylist(const
arraylist_t &alist, arraydata_t &fullad, OAextend const &oaextend,
colindex_t extensioncol, arraylist_t &extensions)

extend a list of arrays ";

%feature("docstring")  extend_arraylist "arraylist_t
extend_arraylist(const arraylist_t &alist, arraydata_t &fullad,
OAextend const &oaextend)

extend a list of arrays ";

%feature("docstring")  extend_arraylist "arraylist_t
extend_arraylist(const arraylist_t &alist, const arraydata_t &fullad)

extend a list of arrays with default options ";

%feature("docstring")  extend_array "arraylist_t extend_array(const
array_link &al, arraydata_t &fullad, OAextend const &oaextend)

extend a single array ";

%feature("docstring")  extend_array "arraylist_t extend_array(const
array_link &al, arraydata_t &arrayclass)

extend a single array with the default LMC algorithm ";

%feature("docstring")  extend_array "int extend_array(carray_t
*array, const arraydata_t *, const colindex_t extensioncol,
arraylist_t &solutions, OAextend const &oaextend)

extend an array with a single column ";

%feature("docstring")  runExtendRoot "arraylist_t
runExtendRoot(arraydata_t adata, int nmax, int verbose=0)

simple wrapper function ";


// File: graphtools_8h.xml
%feature("docstring")  nauty::transformGraph "array_link
transformGraph(const array_link &G, const std::vector< int > tr, int
verbose=1)

apply a vertex permutation to a graph ";

%feature("docstring")  nauty::reduceOAnauty "array_transformation_t
reduceOAnauty(const array_link &al, int verbose=0)

reduce an orthogonal array to Nauty minimal form. the array
transformation is returned ";

%feature("docstring")  nauty::reduceOAnauty "array_transformation_t
reduceOAnauty(const array_link &al, int verbose, const arraydata_t
&ad) ";

%feature("docstring")  nauty::array2graph "std::pair<array_link,
std::vector<int> > array2graph(const array_link &al, int verbose=1)

Convert orthogonal array to graph representation

The conversion method is as in Ryan and Bulutoglu. The resulting graph
is bi-partite. The graph representation can be used for isomorphism
testing. ";

%feature("docstring")  nauty::array2graph "std::pair<array_link,
std::vector<int> > array2graph(const array_link &al, int verbose,
const arraydata_t &ad)

Convert orthogonal array to graph representation

The conversion method is as in Ryan and Bulutoglu. The resulting graph
is bi-partite. The graph representation can be used for isomorphism
testing. ";

%feature("docstring")  nauty::oagraph2transformation "array_transformation_t oagraph2transformation(const std::vector< int >
&pp, const arraydata_t &arrayclass, int verbose=1)

From a relabelling of the graph return the corresponding array
transformation. ";

%feature("docstring")  nauty::unittest_nautynormalform "int
unittest_nautynormalform(const array_link &al, int verbose=1) ";


// File: InfInt_8h.xml


// File: lmc_8h.xml
%feature("docstring")  algorithm_t_list "std::string
algorithm_t_list() ";

%feature("docstring")  algnames "std::string algnames(algorithm_t m)

return name of the algorithm ";

%feature("docstring")  apply_hadamard "void apply_hadamard(const
arraydata_t *ad, array_t *array, colindex_t hcol)

Apply Hadamard transformation to orthogonal array. ";

%feature("docstring")  apply_hadamard "void apply_hadamard(array_link
&al, colindex_t hcol)

Apply Hadamard transformation to orthogonal array. ";

%feature("docstring")  getGlobalStaticIndexed "LMC_static_struct_t*
getGlobalStaticIndexed(int n) ";

%feature("docstring")  cleanGlobalStaticIndexed "void
cleanGlobalStaticIndexed() ";

%feature("docstring")  getGlobalStatic "LMC_static_struct_t*
getGlobalStatic()

return static structure from dynamic global pool, return with
releaseGlobalStatic ";

%feature("docstring")  releaseGlobalStatic "void
releaseGlobalStatic(LMC_static_struct_t *p) ";

%feature("docstring")  cleanGlobalStatic "void cleanGlobalStatic() ";

%feature("docstring")  getGlobalStaticOne "LMC_static_struct_t&
getGlobalStaticOne() ";

%feature("docstring")  rootPlus "array_link rootPlus(const
arraydata_t &ad)

helper function ";

%feature("docstring")  insertUnique "void insertUnique(std::vector<
Type > &cp, const Type &cpv)

helper function ";

%feature("docstring")  check_root_update "int
check_root_update(array_t *array, const arraydata_t &ad)

return 0 if target is equal to original, otherwise return 1 and copy
root initialization + 1 ";

%feature("docstring")  check_root_form "bool check_root_form(const
array_t *array, const arraydata_t &ad)

return true if target is in root form, otherwise return false ";

%feature("docstring")  check_root_update "int
check_root_update(carray_t *original, const arraydata_t &ad, array_t
*target)

return 0 if target is equal to original, otherwise return 1 and copy
root initialization + 1 ";

%feature("docstring")  jj45val "jj45_t jj45val(carray_t *array,
rowindex_t N, int jj, const colperm_t comb, int j5val=-1, int
dosort=1)

return value based on J4-J5 ordering ";

%feature("docstring")  random_transformation "void
random_transformation(array_t *array, const arraydata_t *adp)

Apply a random transformation to an array ";

%feature("docstring")  LMCreduction_train "lmc_t
LMCreduction_train(const array_link &al, const arraydata_t *ad,
LMCreduction_t *reduction, const OAextend &oaextend) ";

%feature("docstring")  LMCreduction_train "lmc_t
LMCreduction_train(const array_t *original, const arraydata_t *ad,
const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend) ";

%feature("docstring")  LMCreduce "lmc_t LMCreduce(array_t const
*original, array_t const *array, const arraydata_t *ad, const
dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend)

helper function ";

%feature("docstring")  LMCreduceFull "lmc_t LMCreduceFull(carray_t
*original, const array_t *array, const arraydata_t *ad, const
dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend, LMC_static_struct_t &tmpStatic)

Perform reduction or LMC check without root trick. ";

%feature("docstring")  LMCcheck "lmc_t LMCcheck(const array_t *array,
const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t
&reduction)

generic LMCcheck function ";

%feature("docstring")  LMCcheck "lmc_t LMCcheck(const array_link &al,
const arraydata_t &ad, const OAextend &oaextend, LMCreduction_t
&reduction)

generic LMCcheck function ";

%feature("docstring")  LMCcheckOriginal "lmc_t LMCcheckOriginal(const
array_link &al)

direct LMC check using the original LMC check ";

%feature("docstring")  calculateSymmetryGroups "LMCreduction_t
calculateSymmetryGroups(const array_link &al, const arraydata_t
&adata, const OAextend &oaextend, int verbose=1, int hack=0)

helper function ";

%feature("docstring")  LMCcheckSymmetryMethod "lmc_t
LMCcheckSymmetryMethod(const array_link &al, const arraydata_t &ad,
const OAextend &oaextend, LMCreduction_t &reduction, LMCreduction_t
&reductionsub, int dverbose) ";

%feature("docstring")  combadd2perm "void combadd2perm(const larray<
numtype > &comb, int newidx, int n, larray< numtype > &target, larray<
numtype > &wtmp)

Convert selection of elements to extended permutation. ";

%feature("docstring")  comb2perm "larray<numtype> comb2perm(const
larray< numtype > comb, int n)

Convert selection of elements to extended permutation. ";

%feature("docstring")  comb2perm "std::vector<numtype>
comb2perm(const std::vector< numtype > comb, int n)

Convert selection of elements to extended permutation. ";

%feature("docstring")  reduceArraysGWLP "void reduceArraysGWLP(const
arraylist_t *input_arrays, arraylist_t &reduced_arrays, int verbose,
int dopruning=1, int strength=2, int dolmc=1)

reduce arrays to canonical form using delete-1-factor ordering ";

%feature("docstring")  reductionDOP "array_transformation_t
reductionDOP(const array_link &al, int verbose=0) ";

%feature("docstring")  selectUniqueArrays "void
selectUniqueArrays(arraylist_t &xlist, arraylist_t &earrays, int
verbose=1) ";

%feature("docstring")  reduceLMCform "array_link reduceLMCform(const
array_link &al)

reduce an array to canonical form using LMC ordering ";

%feature("docstring")  reduceDOPform "array_link reduceDOPform(const
array_link &al, int verbose=0)

reduce an array to canonical form using delete-1-factor ordering ";

%feature("docstring")  LMCcheckLex "std::vector<int>
LMCcheckLex(arraylist_t const &list, arraydata_t const &ad, int
verbose=0)

Apply LMC check (original mode) to a list of arrays ";

%feature("docstring")  LMCcheckLex "lmc_t LMCcheckLex(array_link
const &al, arraydata_t const &ad)

Perform LMC check lexicographically. ";

%feature("docstring")  LMCcheckj4 "lmc_t LMCcheckj4(array_link const
&al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend
&oaextend, int jj=4) ";

%feature("docstring")  LMCcheckj5 "lmc_t LMCcheckj5(array_link const
&al, arraydata_t const &ad, LMCreduction_t &reduction, const OAextend
&oaextend, int hack=0)

Perform LMC check with J5 ordering. ";

%feature("docstring")  root_row_permutation_from_index "void
root_row_permutation_from_index(int permindex, const arraydata_t *ad,
levelperm_t *lperms) ";

%feature("docstring")  create_root_permutations_index "rowperm_t*
create_root_permutations_index(const arraydata_t *ad, int
&totalpermsr) ";

%feature("docstring")  create_root_permutations_index_helper "void
create_root_permutations_index_helper(rowperm_t *rperms, levelperm_t
*lperms, const arraydata_t *ad, int level, int &permcounter) ";

%feature("docstring")  print_rowsort "void print_rowsort(rowsort_t
*rowsort, int N) ";

%feature("docstring")  print_column_rowsort "void
print_column_rowsort(const array_t *arraycol, rowsort_t *rowsort, int
N) ";

%feature("docstring")  print_fracs "void print_fracs(int logl=NORMAL)
";

%feature("docstring")  clear_fracs "void clear_fracs() ";


// File: mathtools_8h.xml
%feature("docstring")  vectormax "Type vectormax(const std::vector<
Type > &v, Type defaultvalue)

Return maximum element of a std::vector. ";

%feature("docstring")  vectormin "Type vectormin(const std::vector<
Type > &v, Type defaultvalue)

Return minimum element of a std::vector. ";

%feature("docstring")  cumsum "std::vector<NumType> cumsum(const
std::vector< NumType > x)

calculate cumulative sum of a vector ";

%feature("docstring")  cumsum0 "std::vector<NumType> cumsum0(const
std::vector< NumType > x)

calculate cumulative sum of a vector with added zero ";

%feature("docstring")  cumsum0 "std::vector<Type>
cumsum0(std::vector< InputType > s) ";

%feature("docstring")  permutation "std::vector<NumType>
permutation(int n)

create permutation of specified length ";

%feature("docstring")  array2vector "std::vector<NumType>
array2vector(const NumTypeIn *x, int len)

convert array given by pointer to std::vector ";

%feature("docstring")  array2larray "larray<NumType>
array2larray(const NumTypeIn *x, int len) ";

%feature("docstring")  print_perm "static void
print_perm(std::ostream &out, const permutationType *s, const int len,
const int maxlen=256)

Print permutation.

Prints a permutation to output stream

Parameters:
-----------

out:  Output stream

s:  Pointer to start of array

len:  Length of array to be printed

maxlen:  (optional) Maximum length to print ";

%feature("docstring")  print_perm "static void
print_perm(std::ostream &out, const larray< permutationType > s, const
int maxlen=256, const bool ret=true) ";

%feature("docstring")  print_perm "static void
print_perm(std::ostream &out, const std::vector< permutationType > s,
const int maxlen=256, const bool ret=true) ";

%feature("docstring")  print_perm_int "static void
print_perm_int(const std::vector< permutationType > s, const int
maxlen=256, const bool ret=true) ";

%feature("docstring")  print_perm "static void print_perm(const
larray< permutationType > s, const int maxlen=256, const bool
ret=true) ";

%feature("docstring")  print_perm "static void print_perm(const
std::vector< permutationType > s, const int maxlen=256, const bool
ret=true) ";

%feature("docstring")  print_perm "static void print_perm(const
permutationType *s, const int len, const int maxlen=256) ";

%feature("docstring")  compare_matrix "int compare_matrix(const
numtype *A, const numtype *B, int r, int c)

Compare two arrays and return whether equal or not.

Parameters:
-----------

A:  Pointer to array

B:  Pointer to array

r:  Number of rows

c:  Number of columns ";

%feature("docstring")  fact "Type fact(const Type f)

Calculates factorial.

A small function that calculates the factorial of a number. Can be
inlined if the compiler decides it is faster. Returns one if the
number is smaller or equal to 1

Parameters:
-----------

f:  Number to calculate the factorial of ";

%feature("docstring")  factorial "static numtype factorial(const
argtype f)

Calculates factorial of type numtype.

This inline factorial function is the same as the standard factorial
calculation, except that the return type is generic

Parameters:
-----------

f:  number to calculate factorial of ";

%feature("docstring")  ncombs "Type ncombs(const Type n, const Type
k)

Calculates number of combinations.

The number of combinations is calculated using the an addapted formula

Parameters:
-----------

n:  Total number of entries to choose from

k:  Number of entries in a certain combination

n!/(k! * (n-k)!) = (n - k + 1) * ..... * n/k! ";

%feature("docstring")  ncombscacheNumber "int ncombscacheNumber() ";

%feature("docstring")  initncombscache "void initncombscache(int N)

initialize datastructure, this function is not thread safe ";

%feature("docstring")  ncombscache "long ncombscache(int n, int k)

return number of combinations from previously calculated results

The results should be calculated with initncombscache ";

%feature("docstring")  ncombsm "Type ncombsm(const Type &n, const
Type &k)

calculate using multiplicative formula, seehttp://en.wikipedia.org/wik
i/Binomial_coefficient#Computing_the_value_ of_binomial_coefficients
";

%feature("docstring")  next_perm "bool next_perm(std::vector<
permutationType > &s) ";

%feature("docstring")  fastrand "int fastrand() ";

%feature("docstring")  seedfastrand "void seedfastrand(int s) ";

%feature("docstring")  fastrandK "int fastrandK(int k) ";

%feature("docstring")  set_srand "void set_srand(unsigned int s)

set the random number seed using srand ";

%feature("docstring")  next_perm_twoperm "void
next_perm_twoperm(permutationType *s, const int len) ";

%feature("docstring")  next_perm "void next_perm(permutationType *s,
const int len) ";

%feature("docstring")  permutationLex "numtype*
permutationLex(numtype k, objecttype *s, numtype n)

See
alsohttp://en.wikipedia.org/wiki/Permutation#Numbering_permutations

Parameters: ";

%feature("docstring")  random_perm "void random_perm(objecttype *s,
numtype len)

Create random permutation using Fisher-Yates shuffle, or Knuth shuffle
";

%feature("docstring")  random_perm "void random_perm(std::vector<
objecttype > &s)

Create random permutation using Fisher-Yates shuffle, or Knuth shuffle
";

%feature("docstring")  new_comb_init "numtype* new_comb_init(int len)

Create a new combination and initialize. ";

%feature("docstring")  delete_comb "void delete_comb(numtype *comb)

Delete combination. ";

%feature("docstring")  init_comb "int init_comb(numtype *comb, int k,
int n)

Initialize a combination

Parameters:
-----------

comb:  Pointer to combination array

k:  Number of elements

n:  Numbers to choose from

Number of combinations possible ";

%feature("docstring")  next_combination "numtype
next_combination(numtype *comb, int k, int n)

Gives combination nr k.

Gives combination number k based on an algorithm from wikipedia.

Parameters:
-----------

comb:  Pointer to combination

k:  Number of the current combination

n:  Number of elements in combination ";

%feature("docstring")  next_combination_fold "numtype
next_combination_fold(numtype *comb, int k, int n) ";

%feature("docstring")  print_combinations "void
print_combinations(int n, int k) ";

%feature("docstring")  perm_is_ordered "int perm_is_ordered(numtype
*perm, int len)

Check whether a permutation is ordered or not.

Parameters:
-----------

perm:

len:

Return 1 of an ordered permutation, 0 otherwise ";

%feature("docstring")  new_perm "numtype* new_perm(int len)

Create a new permutation. ";

%feature("docstring")  clone_perm "numtype* clone_perm(numtype
*source, int len)

Create a new permutation. ";

%feature("docstring")  delete_perm "void delete_perm(numtype *perm)

Delete a permutation. ";

%feature("docstring")  invert_permutation "std::vector<numtype>
invert_permutation(const std::vector< numtype > perm)

Invert a permutation.

Parameters:
-----------

perm:  Permutation as integer type std::vector

New permutation that is the inverse of the argument ";

%feature("docstring")  invert_permutation "void
invert_permutation(const std::vector< numtype > perm, std::vector<
numtype > &iperm)

Invert a permutation.

Parameters:
-----------

perm:  Permutation as integer type std::vector

permout:  Output permutation ";

%feature("docstring")  invert_permutation "void
invert_permutation(numtype *perm, int len, numtype *iperm)

Invert a permutation.

Parameters:
-----------

perm:  Pointer to permutation

len:

Pointer to new permutation that is the inverse of the argument ";

%feature("docstring")  invert_permutation "numtype*
invert_permutation(numtype *perm, int len)

Invert a permutation.

Parameters:
-----------

perm:  Pointer to permutation

len:

Pointer to new permutation that is the inverse of the argument ";

%feature("docstring")  perform_level_perm "void
perform_level_perm(numtype *src, int n, const numtype *perm)

Perform level permutation on an array

Parameters:
-----------

src:  Pointer to array

n:  Length of array

perm:  Permutation to perform ";

%feature("docstring")  composition_perm "void composition_perm(const
numtype *A, const numtype *B, int n, numtype *C)

Calculate composition of 2 permutations.

Calculates C = B  A

Parameters:
-----------

A:

B:

n:

C:  ";

%feature("docstring")  composition_perm "void composition_perm(const
std::vector< numtype > &A, const std::vector< numtype > &B,
std::vector< numtype > &C)

Calculate composition of 2 permutations.

Calculates C = B  A

Parameters:
-----------

A:

B:

n:

C:  ";

%feature("docstring")  perform_perm "void perform_perm(const object
*const src, object *const target, const int n, const numtype *perm)

Perform a permutation on a set of objects.

Parameters:
-----------

src:

target:

n:

perm:  ";

%feature("docstring")  perform_perm "std::vector<object>
perform_perm(const std::vector< object > src, const std::vector<
numtype > perm)

Perform a permutation on a set of objects. ";

%feature("docstring")  perform_inv_perm "void perform_inv_perm(const
std::vector< object > src, std::vector< object > &target, const int n,
const std::vector< numtype > perm)

Perform inverse permutation. ";

%feature("docstring")  perform_inv_perm "void perform_inv_perm(const
object *const src, object *const target, const int n, const numtype
*perm)

Perform inverse permutation. ";

%feature("docstring")  perform_inv_perm "void perform_inv_perm(const
object *const src, object *const target, const int n, const
std::vector< numtype > &perm)

Perform inverse permutation. ";

%feature("docstring")  perform_level_perm "void
perform_level_perm(const numtype *const src, numtype *const target,
const int n, const numtype *perm)

Perform a permutation on a set of data elements.

Parameters:
-----------

src:

target:

n:

perm:  ";

%feature("docstring")  init_perm "void init_perm(std::vector< numtype
> &perm)

Initialiaze a permutation

Parameters:
-----------

perm:

len:  ";

%feature("docstring")  init_perm "void init_perm(numtype *perm, int
len)

Initialiaze a permutation

Parameters:
-----------

perm:

len:  ";

%feature("docstring")  init_signperm "void init_signperm(std::vector<
numtype > &signperm)

Initialiaze a permutation

Parameters:
-----------

perm:

len:  ";

%feature("docstring")  compare_perm "bool compare_perm(const numtype
*permA, const numtype *permB, int len) ";

%feature("docstring")  copy_perm "void copy_perm(const numtype
*source, numtype *target, int len)

copy a permuntation ";

%feature("docstring")  init_perm_n "outtype init_perm_n(numtype
*perm, int len)

initialize a permutation and return the number of permutations ";

%feature("docstring")  new_perm_init "numtype* new_perm_init(int len)
";

%feature("docstring")  issorted "bool issorted(_ForwardIterator
first, const _ForwardIterator last) ";

%feature("docstring")  new_valueindex "returntype*
new_valueindex(const basetype *bases, const numtype n) ";

%feature("docstring")  init_valueindex_forward "numtype*
init_valueindex_forward(numtype *valueindex, const numtype *bases,
const numtype n) ";

%feature("docstring")  init_valueindex "numtype*
init_valueindex(numtype *valueindex, const numtype *bases, const
numtype n) ";

%feature("docstring")  argsort "std::vector<int> argsort(const
std::vector< numtype > vv) ";

%feature("docstring")  permute "std::vector<Type> permute(const
std::vector< Type > x, const std::vector< IndexType > indices)

Permute a std::vector. ";

%feature("docstring")  permuteback "std::vector<Type>
permuteback(const std::vector< Type > x, const std::vector< IndexType
> indices)

Permute a std::vector with inverse permutation. ";

%feature("docstring")  symm_group_index_plain "int
symm_group_index_plain(const numtype *vec, const int n, itype *&idx,
itype *&gstart, itype *&gsize)

Calculate symmetry groups of a list of integers under permutations.

Parameters:
-----------

vec:

n:

idx:

gstart:

gsize:

Number of groups found ";

%feature("docstring")  ipow "long ipow(long x, long y)

Power of two integers. ";

%feature("docstring")  ipow "unsigned int ipow(unsigned int x,
unsigned int p)

Power of two unsigned integers. ";

%feature("docstring")  powmo "int powmo(int n)

-1 to the power n (integer) ";

%feature("docstring")  krawtchouk "Type krawtchouk(Type j, Type x,
Type n, Type s, int verbose=0)

calculate value of Krawtchouk polynomial ";

%feature("docstring")  krawtchouksCache "Type krawtchouksCache(Type
j, Type x, Type n)

calculate value of Krawtchouk polynomial ";

%feature("docstring")  krawtchouks "Type krawtchouks(Type j, Type x,
Type n)

calculate value of Krawtchouk polynomial ";

%feature("docstring")  conditionNumber "double conditionNumber(const
Eigen::Matrix< Type,-1,-1 > A)

return the condition number of a matrix ";


// File: md5_8h.xml
%feature("docstring")  md5 "std::string md5(void *data, int numbytes)

calculate md5 sum of a data block in memory ";

%feature("docstring")  md5 "std::string md5(const std::string
filename)

calculate md5 sum of a file on disk ";


// File: msstdint_8h.xml


// File: nonroot_8h.xml
%feature("docstring")  LMCreduce_non_root "lmc_t
LMCreduce_non_root(const array_t *original, const arraydata_t *ad,
dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend, const LMC_static_struct_t &tmpStatic)

default reduction function for non-root stage ";

%feature("docstring")  LMCreduce_non_root_j4 "lmc_t
LMCreduce_non_root_j4(const array_t *original, const arraydata_t *ad,
const dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend, LMC_static_struct_t &tmpStatic) ";

%feature("docstring")  LMCreduce_non_root_2level "lmc_t
LMCreduce_non_root_2level(const array_t *original, const arraydata_t
*ad, dyndata_t *dyndata, LMCreduction_t *reduction, const OAextend
&oaextend, const LMC_static_struct_t &tmpStatic) ";

%feature("docstring")  LMC_check_col_rowsymm "lmc_t
LMC_check_col_rowsymm(const array_t *arraycol, const arraydata_t *ad,
const symmdata &sd, int col, int dverbose=0) ";


// File: oaoptions_8h.xml
%feature("docstring")  compile_information "std::string
compile_information()

Print the compile-time options to string. String with information ";

%feature("docstring")  version "std::string version()

Print version. ";

%feature("docstring")  print_copyright "void print_copyright()

Print copyright statement. ";

%feature("docstring")  print_copyright_old "void
print_copyright_old()

Print copyright statement. ";

%feature("docstring")  print_copyright_light "void
print_copyright_light()

Print copyright statement. ";

%feature("docstring")  print_options "void print_options(std::ostream
&outx)

Print compile time options. ";

%feature("docstring")  print_options "void print_options() ";

%feature("docstring")  oadevelop "int oadevelop() ";


// File: pareto_8h.xml


// File: printfheader_8h.xml


// File: strength_8h.xml
%feature("docstring")  new_strength_freq_table "strength_freq_table
new_strength_freq_table(int ncolcombs, int *nvalues, int &nelements)

Constructor.

Create a table with strength frequencies

Parameters:
-----------

ncolcombs:  Number of column combinations to store

nvalues:  Number of tuples that can occur for each column combination

nelements:  Return variable with the size of the allocated array ";

%feature("docstring")  create_reverse_colcombs_fixed "rev_index*
create_reverse_colcombs_fixed(const int ncolcombs) ";

%feature("docstring")  create_reverse_colcombs "rev_index*
create_reverse_colcombs(colindex_t **colcombs, const int ncols, const
int strength) ";

%feature("docstring")  free_colcombs_fixed "void
free_colcombs_fixed(colindex_t **colcombs, int *lambda, int *nvalues)
";

%feature("docstring")  copy_freq_table "void
copy_freq_table(strength_freq_table source, strength_freq_table
target, int ftsize)

Copy the frequency count table. ";

%feature("docstring")  check_divisibility "bool
check_divisibility(const arraydata_t *)

check whether an array passes divisibility test ";

%feature("docstring")  print_frequencies "void print_frequencies(int
**frequencies, const int nelements, const int *lambda, const int N) ";

%feature("docstring")  set_colcombs_fixed "colindex_t**
set_colcombs_fixed(int *&xlambda, int *&nvalues, int &ncolcombs, const
array_t *s, const int strength, const int fixedcol, const int N)

Return all column combinations including a fixed column. At the same
time allocate space for the number of values these columns have

Parameters:
-----------

xlambda:

nvalues:

ncolcombs:

s:

strength:

fixedcol:

N:  ";

%feature("docstring")  add_element_freqtable "void
add_element_freqtable(extend_data_t *es, rowindex_t activerow,
carray_t *array, strength_freq_table freqtable)

Add row to frequency table using cache system. ";

%feature("docstring")  add_element_freqtable_col "void
add_element_freqtable_col(extend_data_t *es, rowindex_t activerow,
carray_t *arraycol, strength_freq_table freqtable)

fast version of add_element_freqtable ";

%feature("docstring")  init_frequencies "void
init_frequencies(extend_data_t *es, array_t *array)

Initialize the table of t-tuple frequencies. ";

%feature("docstring")  recount_frequencies "void
recount_frequencies(int **frequencies, extend_data_t *es, colindex_t
currentcol, rowindex_t rowstart, rowindex_t rowlast, carray_t *array)
";

%feature("docstring")  set_indices "vindex_t** set_indices(colindex_t
**colcombs, basetype *bases, const int k, colindex_t ncolcombs)

Helper function. ";

%feature("docstring")  strength_check "bool strength_check(const
arraydata_t &ad, const array_link &al, int verbose=1)

perform strength check on an array

Special case for extension of an array with proper strength ";

%feature("docstring")  strength_check "bool strength_check(const
array_link &al, int strength, int verbose=0)

perform strength check on an array ";

%feature("docstring")  valid_element "bool valid_element(const
extend_data_t *es, const extendpos *p, carray_t *array)

Determine whether an element passes the strength test.

Parameters:
-----------

es:

p:

array:  ";

%feature("docstring")  valid_element_2level "bool
valid_element_2level(const extend_data_t *es, const extendpos *p)

Determine whether an element passes the strength test, specialized for
2-level array ";


// File: timsort_8hpp.xml


// File: tools_8h.xml
%feature("docstring")  base_name "std::string base_name(std::string
const &path) ";

%feature("docstring")  printfd_handler "void printfd_handler(const
char *file, const char *func, int line, const char *message,...) ";

%feature("docstring")  log_print "int log_print(const int level,
const char *message,...) ";

%feature("docstring")  getloglevel "int getloglevel() ";

%feature("docstring")  setloglevel "void setloglevel(int n) ";

%feature("docstring")  checkloglevel "bool checkloglevel(int l) ";

%feature("docstring")  logstream "std::ostream& logstream(int level)
";

%feature("docstring")  system_uname "std::string system_uname() ";

%feature("docstring")  mycheck_handler "void mycheck_handler(const
char *file, const char *func, int line, int condition, const char
*message,...) ";

%feature("docstring")  myassert "void myassert(int condition, const
char *str=0) ";

%feature("docstring")  myassertdebug2 "void myassertdebug2(int
condition, const char *str) ";

%feature("docstring")  cprintf "int cprintf(int check, const char
*message,...) ";

%feature("docstring")  ff "void ff()

flush to stdout ";

%feature("docstring")  safedelete "void safedelete(A *p)

Delete a pointer and set to zero. ";

%feature("docstring")  safedeletearray "void safedeletearray(A *p)

Delete array and set pointer to zero

Parameters:
-----------

p:  ";

%feature("docstring")  next_comb "int next_comb(std::vector< Type >
&comb, int k, int n)

Gives next combination.

Gives next combination for k elements out of n based on an algorithm
from wikipedia. The generate is sorted.

Parameters:
-----------

comb:  Pointer to combination

k:  Number of the current combination

n:  Number of elements in combination ";

%feature("docstring")  next_comb "int next_comb(int *comb, int k, int
n)

Go to next combination in sequence ";

%feature("docstring")  next_comb_s "int next_comb_s(int *comb, int k,
int n)

Go to next combination in sequence ";

%feature("docstring")  row_rank "static int row_rank(array_t *array,
const int n_columns, const int n_rows, const int *index)

Calculates row rank.

Calculates the rank of a row, in order to determine the sorting
order.Todo Combine n_columns and n_rows into *nrs, to reduce overhead

Parameters:
-----------

array:

n_columns:

n_rows:

index:  ";

%feature("docstring")  row_rank_partial "static int
row_rank_partial(const carray_t *array, rowindex_t n_rows, const
vindex_t *index, colindex_t start_idx, colindex_t end_idx, const
colperm_t &colperm, rowindex_t row)

Returns the value of (part of) a row.

Parameters:
-----------

array:

start_idx:

end_idx:

n_rows:

index:  Value index for each of the columns of the array ";

%feature("docstring")  row_rank_partial "static array_t
row_rank_partial(carray_t *array, colindex_t start_idx, colindex_t
end_idx, rowindex_t n_rows, const vindex_t *index)

Returns the value of (part of) a row.

Parameters:
-----------

array:

start_idx:

end_idx:

n_rows:

index:  Value index for each of the columns of the array ";

%feature("docstring")  row_rank_partial "static array_t
row_rank_partial(carray_t *array, const colindex_t start_idx, const
colindex_t end_idx, const rowindex_t row, const rowindex_t n_rows,
const int *index)

Returns the value of (part of) a row.

Parameters:
-----------

array:

start_idx:

end_idx:

n_rows:

index:  Value index for each of the columns of the array ";

%feature("docstring")  swap_object "void swap_object(Object &a,
Object &b)

Template to swap two objects of arbitrary datatype Please use
std::swap instead.

Parameters:
-----------

a:

b:  ";

%feature("docstring")  malloc2d_nelements "int
malloc2d_nelements(const int nrows, const rtype *rowsizes)

Calculate the number of elements in a 2D table with rows with
different sizes ";

%feature("docstring")  malloc2d_irr "DataType** malloc2d_irr(const
int nrows, const rtype *rowsizes)

Allocate a 2-dimensional array with non-uniform rows.

Parameters:
-----------

nrows:  Number of rows in the table

rowsizes:  Size of each row ";

%feature("docstring")  malloc2d_irr "DataType** malloc2d_irr(const
int nrows, const rtype *rowsizes, int &nelements)

Allocate a 2-dimensional array with non-uniform rows, return size of
allocated space.

Parameters:
-----------

nrows:  Number of rows in the table

rowsizes:  Size of each row

nelements:  This parameter is initialized with the size of the array
allocated ";

%feature("docstring")  malloc2d "DataType** malloc2d(const numtype
nrows, const int rowsize)

Allocate a 2-dimensional array of specified size.

Parameters:
-----------

nrows:

rowsize:  ";

%feature("docstring")  free2d "void free2d(DataType **data, const int
nrows)

Release a 2-dimensional array.

Parameters:
-----------

data:

nrows:  ";

%feature("docstring")  free2d "void free2d(DataType **data)

Release a 2-dimensional array.

Parameters:
-----------

data:  ";

%feature("docstring")  free2d_irr "void free2d_irr(DataType **data)

Release a 2-dimensional non-uniform array.

Parameters:
-----------

data:  ";

%feature("docstring")  free2d_irr "void free2d_irr(DataType **data,
const int nrows)

Release a 2-dimensional non-uniform array.

Parameters:
-----------

data:

nrows:  ";

%feature("docstring")  print_array "void print_array(const char *str,
const array_t *array, const rowindex_t r, const colindex_t c) ";

%feature("docstring")  print_array "void print_array(const array_t
*array, const rowindex_t r, const colindex_t c) ";

%feature("docstring")  print_array "void print_array(const array_link
&A)

Print array to stdout. ";

%feature("docstring")  display_vector "void display_vector(const
std::vector< atype > &v)

print vector ";

%feature("docstring")  printf_vector "void printf_vector(const
std::vector< atype > &v, const char *format)

print vector ";

%feature("docstring")  show_array_dyn "void show_array_dyn(const
atype *array, const int x, const int y) ";

%feature("docstring")  countelements "void countelements(carray_t
*array, const int nelements, const int maxval, int *elements)

Counts the number of occurences of each value in an array. ";

%feature("docstring")  addelement "void addelement(const array_t
elem, int *elements)

Add element to element counter.

Parameters:
-----------

elem:

elements:  ";

%feature("docstring")  get_time_ms "double get_time_ms()

return time with milisecond precision ";

%feature("docstring")  get_time_ms "double get_time_ms(double t0)

return time difference with milisecond precision ";

%feature("docstring")  trim "void trim(std::string &str, const
std::string &trimChars=\"\")

trim a string by removing the specified characters from the left and
right ";

%feature("docstring")  currenttime "std::string currenttime()

return the current time as a string ";

%feature("docstring")  oafilestring "std::string oafilestring(const
arraydata_t *ad)

return string describing array ";

%feature("docstring")  oafilestring "std::string
oafilestring(rowindex_t rows, colindex_t cols, array_t *s)

return string describing array ";

%feature("docstring")  itos "std::string itos(numtype i)

Convert integer to C++ string.

Parameters:
-----------

i:  Integer

String representation of the integer ";

%feature("docstring")  printfstring "std::string printfstring(const
char *message,...)

printf-style function that returns std::string

printf-style function that returns std::string

Parameters:
-----------

message:  ";

%feature("docstring")  printtime "std::string printtime() ";

%feature("docstring")  insertionSort "void insertionSort(Object x[],
int length)

Tempalate for insertionSort.

Parameters:
-----------

x[]:

length:  ";

%feature("docstring")  bubbleSort "void bubbleSort(itemType a[],
indexType l, indexType r)

sort arrays using bubbleSort ";

%feature("docstring")  flipSort "void flipSort(itemType a[],
indexType l, indexType r)

sorting similar to bubblesort but fast for sorted arrays

The indices l and r are inclusive. ";

%feature("docstring")  bubbleSort2 "void bubbleSort2(Object obj[],
indexType array_size)

Template for bubble sort.

Parameters:
-----------

obj[]:

array_size:  ";

%feature("docstring")  quickSort "void quickSort(T a[], const int
&leftarg, const int &rightarg)

sort list using quickSort ";

%feature("docstring")  shellSort "void shellSort(itemType a[],
indexType l, indexType r)

sort list using shellSort ";

%feature("docstring")  replaceString "std::string
replaceString(std::string subject, const std::string &search, const
std::string &replace)

replace all occurces of a substring in a string ";

%feature("docstring")  printdoubleasbits "void
printdoubleasbits(double decker)

print a double value as bits ";


// File: version_8h.xml


// File: dir_68267d1309a1af8e8297ef4c3efbcdba.xml


// File: indexpage.xml

