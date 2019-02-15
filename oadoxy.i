
// File: index.xml

// File: struct__my__lldiv__t.xml


%feature("docstring") _my_lldiv_t "

C++ includes: InfInt.h
";

// File: classAnyOption.xml


%feature("docstring") AnyOption "

C++ includes: anyoption.h
";

%feature("docstring") AnyOption::AnyOption "
";

%feature("docstring") AnyOption::AnyOption "
";

%feature("docstring") AnyOption::AnyOption "
";

%feature("docstring") AnyOption::~AnyOption "
";

%feature("docstring") AnyOption::setCommandPrefixChar "
";

%feature("docstring") AnyOption::setCommandLongPrefix "
";

%feature("docstring") AnyOption::setFileCommentChar "
";

%feature("docstring") AnyOption::setFileDelimiterChar "
";

%feature("docstring") AnyOption::useCommandArgs "
";

%feature("docstring") AnyOption::useFileName "
";

%feature("docstring") AnyOption::noPOSIX "
";

%feature("docstring") AnyOption::setVerbose "
";

%feature("docstring") AnyOption::setOption "
";

%feature("docstring") AnyOption::setOption "
";

%feature("docstring") AnyOption::setOption "
";

%feature("docstring") AnyOption::setFlag "
";

%feature("docstring") AnyOption::setFlag "
";

%feature("docstring") AnyOption::setFlag "
";

%feature("docstring") AnyOption::setCommandOption "
";

%feature("docstring") AnyOption::setCommandOption "
";

%feature("docstring") AnyOption::setCommandOption "
";

%feature("docstring") AnyOption::setCommandFlag "
";

%feature("docstring") AnyOption::setCommandFlag "
";

%feature("docstring") AnyOption::setCommandFlag "
";

%feature("docstring") AnyOption::setFileOption "
";

%feature("docstring") AnyOption::setFileOption "
";

%feature("docstring") AnyOption::setFileOption "
";

%feature("docstring") AnyOption::setFileFlag "
";

%feature("docstring") AnyOption::setFileFlag "
";

%feature("docstring") AnyOption::setFileFlag "
";

%feature("docstring") AnyOption::processOptions "
";

%feature("docstring") AnyOption::processCommandArgs "
";

%feature("docstring") AnyOption::processCommandArgs "
";

%feature("docstring") AnyOption::processCommandArgs "
";

%feature("docstring") AnyOption::processCommandArgs "
";

%feature("docstring") AnyOption::processFile "
";

%feature("docstring") AnyOption::processFile "
";

%feature("docstring") AnyOption::getValue "
";

%feature("docstring") AnyOption::getValue "
";

%feature("docstring") AnyOption::getFlag "
";

%feature("docstring") AnyOption::getFlag "
";

%feature("docstring") AnyOption::getDoubleValue "
";

%feature("docstring") AnyOption::getDoubleValue "
";

%feature("docstring") AnyOption::getIntValue "
";

%feature("docstring") AnyOption::getIntValue "
";

%feature("docstring") AnyOption::getLongValue "
";

%feature("docstring") AnyOption::getStringValue "
";

%feature("docstring") AnyOption::getStringValue "
";

%feature("docstring") AnyOption::printUsage "
";

%feature("docstring") AnyOption::printAutoUsage "
";

%feature("docstring") AnyOption::addUsage "
";

%feature("docstring") AnyOption::addUsage "
";

%feature("docstring") AnyOption::printHelp "
";

%feature("docstring") AnyOption::autoUsagePrint "
";

%feature("docstring") AnyOption::getArgc "
";

%feature("docstring") AnyOption::getArgv "
";

%feature("docstring") AnyOption::hasOptions "
";

// File: structarray__link.xml


%feature("docstring") array_link "

C++ includes: arraytools.h
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

The array is intialized with zeros.  

Parameters
----------
* `nrows` :  
    Number of rows  
* `ncols` :  
    Number of columns  
* `index` :  
    Number to keep track of lists of designs  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

Initialize with data from a pointer.  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

Initialize with data from another array_link object.  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

Initialize with data from an Eigen matrix.  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

The array is initialized by permuting the columns of another array  

Parameters
----------
* `array` :  
    Source to copy from  
* `column_permutation` :  
    The permuntation to apply  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  
";

%feature("docstring") array_link::array_link "

A class representing an integer valued array  

The array is initialized by copying the values from a vector.  
";

%feature("docstring") array_link::~array_link "
";

%feature("docstring") array_link::clone "
";

%feature("docstring") array_link::showarray "

print array to stdout  
";

%feature("docstring") array_link::showarrayString "

print array to string  
";

%feature("docstring") array_link::showarraycompact "

print array to stdout in compact format (no whitespace between elemenents)  
";

%feature("docstring") array_link::showproperties "

print array properties to stdout  
";

%feature("docstring") array_link::is2level "

return true if the array is a 2-level array (e.g. only contains values 0 and 1)  
";

%feature("docstring") array_link::is_mixed_level "

return true is the array is a mixel-level array  
";

%feature("docstring") array_link::is_orthogonal_array "

return true is the array is array with values in 0, 1, ..., for each column  
";

%feature("docstring") array_link::is_conference "

return true if the array is a +1, 0, -1 valued array  
";

%feature("docstring") array_link::is_conference "

return true if the array is a +1, 0, -1 valued array, with specified number of
zeros in each column  
";

%feature("docstring") array_link::isSymmetric "

return true if the array is symmetric  
";

%feature("docstring") array_link::makeSymmetric "

make the array symmetric by copying the upper-right to the lower-left  
";

%feature("docstring") array_link::deleteColumn "

return array with selected column removed  
";

%feature("docstring") array_link::selectFirstRows "

return array with first number_of_arrays rows  
";

%feature("docstring") array_link::selectFirstColumns "

return array with first number_of_arrays columns selected  
";

%feature("docstring") array_link::selectLastColumns "

return array with last number_of_arrays columns selected  
";

%feature("docstring") array_link::selectColumns "

select columns from an array  
";

%feature("docstring") array_link::selectColumns "

select single column from an array  
";

%feature("docstring") array_link::setColumn "

set a column of the array to the given vector  
";

%feature("docstring") array_link::setColumn "

set a column of the array to the given vector  
";

%feature("docstring") array_link::transposed "

return transposed array  
";

%feature("docstring") array_link::Defficiency "

calculate D-efficiency  
";

%feature("docstring") array_link::DsEfficiency "

calculate main effect robustness (or Ds-optimality)  
";

%feature("docstring") array_link::Defficiencies "

calculate D-efficiency, calculate main effect robustness (or Ds-optimality) and
D1-efficiency for an orthogonal array  
";

%feature("docstring") array_link::VIFefficiency "
";

%feature("docstring") array_link::Aefficiency "

calculate A-efficiency  
";

%feature("docstring") array_link::Eefficiency "

calculate E-efficiency  
";

%feature("docstring") array_link::Fvalues "

Calculate F-values of a 2-level matrix.  

This assumes the strength is at least 3. Otherwise use the jstruct_t object  
";

%feature("docstring") array_link::FvaluesConference "

Calculate F-values of a conference design  

       \\param number_of_columns Number of columns to use
       \\return The Fk vector with k the number of columns specified
  
";

%feature("docstring") array_link::Jcharacteristics "

Calculate the Jk-characteristics of the matrix (the values are signed)  

Parameters
----------
* `jj` :  
    Number of columns to use  

Returns
-------
Vector with calculated Jk values  
";

%feature("docstring") array_link::PECsequence "

Calculate the projective estimation capacity sequence.  
";

%feature("docstring") array_link::PICsequence "

Calculate the projective information capacity sequence.  
";

%feature("docstring") array_link::rank "

calculate rank of array  
";

%feature("docstring") array_link::GWLP "

Calculate generalized wordlength pattern  

See also: GWLP  
";

%feature("docstring") array_link::strength "

calculate strength of an array  
";

%feature("docstring") array_link::foldover "

return true if the array is a foldover array  
";

%feature("docstring") array_link::min "
";

%feature("docstring") array_link::max "
";

%feature("docstring") array_link::CL2discrepancy "

Calculate centered L2 discrepancy  

The method is from \"A connection between uniformity and aberration in regular
fractions of two-level factorials\", Fang and Mukerjee, 2000  
";

%feature("docstring") array_link::randomperm "

apply a random permutation of rows, columns and levels of an orthogonal array  
";

%feature("docstring") array_link::randomcolperm "

apply a random permutation of columns of an orthogonal array  
";

%feature("docstring") array_link::randomrowperm "

apply a random permutation of rows of an orthogonal array  
";

%feature("docstring") array_link::getModelMatrix "

Caculate model matrix of an orthogonal array  

Parameters
----------
* `order` :  
    For 0 return only the intercept; for 1 return intercept and main effects;
    for 2 return intercept, main effects and interaction effects.  
* `intercept` :  
    If 1, then include the intercept in the output.  
* `verbose` :  
    Verbosity level  

Returns
-------
Calculated model matrix  

This function uses array2eigenModelMatrixMixed for the calculation.  
";

%feature("docstring") array_link::deepcopy "
";

%feature("docstring") array_link::shallowcopy "
";

%feature("docstring") array_link::equalsize "

return true of two array have the same dimensions  
";

%feature("docstring") array_link::atfast "

get element from array, no error checking, inline version  
";

%feature("docstring") array_link::atfast "

get element from array, no error checking, inline version  
";

%feature("docstring") array_link::_at "

get element at specified position, no bounds checking  
";

%feature("docstring") array_link::_at "

get element at specified position, no bounds checking  
";

%feature("docstring") array_link::at "

get element at specified position  
";

%feature("docstring") array_link::at "

get element at specified position  
";

%feature("docstring") array_link::at "

get element at specified position  
";

%feature("docstring") array_link::setconstant "

set all elements in the array to a value  
";

%feature("docstring") array_link::setvalue "

set value of an array  
";

%feature("docstring") array_link::setvalue "

set value of an array  
";

%feature("docstring") array_link::_setvalue "

set value of an array, no bounds checking!  
";

%feature("docstring") array_link::negateRow "

multiply a row by -1  
";

%feature("docstring") array_link::show "

print information about array  
";

%feature("docstring") array_link::showstr "

return string describing the array  
";

%feature("docstring") array_link::md5 "

return md5 sum of array representation (as represented with 32bit int datatype
in memory)  
";

%feature("docstring") array_link::columnEqual "

return true if two columns are equal  
";

%feature("docstring") array_link::firstColumnDifference "

return index of first different column  
";

%feature("docstring") array_link::firstDiff "

Calculate row and column index of first difference between two arrays  

The difference is according to the column-major ordering.  
";

%feature("docstring") array_link::create_root "

create root in arraylink  
";

%feature("docstring") array_link::nonzero_fraction "

return fraction of nonzero elements in array  
";

%feature("docstring") array_link::clear "

fill array with zeros  
";

%feature("docstring") array_link::getarraydata "
";

%feature("docstring") array_link::setarraydata "

internal function  
";

%feature("docstring") array_link::setarraydata "

special method for SWIG interface  
";

%feature("docstring") array_link::setarraydata "

internal function  
";

%feature("docstring") array_link::setcolumn "

set column to values  
";

%feature("docstring") array_link::init "
";

%feature("docstring") array_link::row_symmetry_group "

return the row_symmetry group of an array  
";

%feature("docstring") array_link::reduceLMC "

return the LMC form of the array  
";

%feature("docstring") array_link::reduceDOP "

return the delete-one-factor-projection form of the array  
";

%feature("docstring") array_link::getEigenMatrix "

return the array as an Eigen matrix  
";

%feature("docstring") array_link::columnGreater "

return true of specified column is smaller than column in another array  
";

%feature("docstring") array_link::debug "
";

// File: classarray__transformation__t.xml


%feature("docstring") array_transformation_t "

Contains a transformation of an array.  

Contains an array transformation. The transformation consists of column, row and
level permutations. The level and column permutations are not commutative (since
the level permutations are tied to a particular column). We apply the column
permutations first.  

C++ includes: arraytools.h
";

%feature("docstring") array_transformation_t::array_transformation_t "
";

%feature("docstring") array_transformation_t::array_transformation_t "
";

%feature("docstring") array_transformation_t::array_transformation_t "
";

%feature("docstring") array_transformation_t::array_transformation_t "

copy constructor  
";

%feature("docstring") array_transformation_t::~array_transformation_t "
";

%feature("docstring") array_transformation_t::show "

show the array transformation  
";

%feature("docstring") array_transformation_t::show "
";

%feature("docstring") array_transformation_t::isIdentity "

return true if the transformation is equal to the identity  
";

%feature("docstring") array_transformation_t::inverse "

return the inverse transformation  
";

%feature("docstring") array_transformation_t::reset "

return the transformation to the identity transformation  
";

%feature("docstring") array_transformation_t::randomize "

initialize to a random transformation  
";

%feature("docstring") array_transformation_t::randomizecolperm "

initialize with a random column permutation  
";

%feature("docstring") array_transformation_t::randomizerowperm "

initialize with a random row permutation  
";

%feature("docstring") array_transformation_t::apply "

apply transformation to an array_link object  
";

%feature("docstring") array_transformation_t::apply "

apply transformation to an array (inplace)  
";

%feature("docstring") array_transformation_t::apply "

apply transformation to an array  
";

%feature("docstring") array_transformation_t::print_transformed "

apply transformation and show resulting array  
";

%feature("docstring") array_transformation_t::rowperm "

return the row permutation of the transformation  
";

%feature("docstring") array_transformation_t::colperm "

return the column permutation of the transformation  
";

%feature("docstring") array_transformation_t::lvlperm "

return the level permutations of the transformation  
";

%feature("docstring") array_transformation_t::setrowperm "

set the row permutation of the transformation  
";

%feature("docstring") array_transformation_t::setcolperm "

set the column permutation of the transformation  
";

%feature("docstring") array_transformation_t::setlevelperm "

set the level permutation of the transformation  
";

// File: structarraydata__t.xml


%feature("docstring") arraydata_t "

Specifies a class of arrays.  

The specification includes the number of rows, number of columns, factor levels
and strength.  

C++ includes: arraytools.h
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  

Parameters
----------
* `s` :  
    Factor levels  
* `N` :  
    Number of rows  
* `strength` :  
    Strength for class  
* `ncols` :  
    Number of columns for the class  
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  

Parameters
----------
* `s` :  
    Factor levels  
* `N` :  
    Number of rows  
* `strength` :  
    Strength for class  
* `ncols` :  
    Number of columns for the class  
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  
";

%feature("docstring") arraydata_t::arraydata_t "

Specifies a class of orthogonal arrays  

The specification includes the number of rows, number of columns, factor levels
and strength.  

An orthogonal array of strength t, N runs, k factors (columns) and factor levels
s[i] is an N times k array with symbols 0, 1, ..., s[i]-1 in column i such that
for every t columns every t-tuple of elements occurs equally often.  
";

%feature("docstring") arraydata_t::~arraydata_t "
";

%feature("docstring") arraydata_t::ismixed "

return true if the class represents mixed-level arrays  
";

%feature("docstring") arraydata_t::is2level "

return true if the class represents a 2-level array  
";

%feature("docstring") arraydata_t::randomarray "

return random array from the class. this operation is only valid for strength 0
or 1  
";

%feature("docstring") arraydata_t::writeConfigFile "

Write file with specification of orthognal array class.  

Parameters
----------
* `filename` :  
    Filename to write to  
";

%feature("docstring") arraydata_t::idstr "
";

%feature("docstring") arraydata_t::idstrseriesfull "
";

%feature("docstring") arraydata_t::fullidstr "
";

%feature("docstring") arraydata_t::latexstr "

return latex string describing the class  
";

%feature("docstring") arraydata_t::reduceColumns "
";

%feature("docstring") arraydata_t::showstr "
";

%feature("docstring") arraydata_t::show "
";

%feature("docstring") arraydata_t::complete_arraydata "

Calculate derived data such as the index and column groups from a design.  
";

%feature("docstring") arraydata_t::lmc_overflow_check "

check whether the LMC calculation will overflow  
";

%feature("docstring") arraydata_t::complete_arraydata_fixlast "
";

%feature("docstring") arraydata_t::complete_arraydata_splitn "
";

%feature("docstring") arraydata_t::set_colgroups "
";

%feature("docstring") arraydata_t::set_colgroups "

set column group equal to that of a symmetry group  
";

%feature("docstring") arraydata_t::show_colgroups "

show column groups in the array class  
";

%feature("docstring") arraydata_t::calculate_oa_index "

calculate the index of the orthogonal arrays in this class  
";

%feature("docstring") arraydata_t::create_root "

return the root array for the class  
";

%feature("docstring") arraydata_t::getfactorlevel "

return the factor level for the specified column return -1 if the column index
is invalid  
";

%feature("docstring") arraydata_t::getS "

return factor levels  
";

%feature("docstring") arraydata_t::factor_levels "

return factor levels  
";

%feature("docstring") arraydata_t::reset_strength "

Reset strength of arraydata.  

Parameters
----------
* `strength` :  
    The strength to reset the structure to  
";

%feature("docstring") arraydata_t::get_col_group "

Return index of the column group for a column.  
";

%feature("docstring") arraydata_t::is_factor_levels_sorted "

Return True if the factor levels are sorted from large to small.  
";

// File: structarrayfile_1_1arrayfile__t.xml


%feature("docstring") arrayfile::arrayfile_t "

Structure for reading or writing a file with arrays.  

The format of the file is determined by the `arrayfilemode_t` The format
described in detail in the documentation of the OApackage
https://oapackage.readthedocs.io/en/latest/.  

C++ includes: arraytools.h
";

%feature("docstring") arrayfile::arrayfile_t::arrayfile_t "

Structure for reading or writing a file with arrays  
";

%feature("docstring") arrayfile::arrayfile_t::arrayfile_t "

Structure for reading or writing a file with arrays  

Parameters
----------
* `filename` :  
    File to open for reading  
* `verbose` :  
    Verbosity level  
";

%feature("docstring") arrayfile::arrayfile_t::arrayfile_t "

Structure for reading or writing a file with arrays  

Open new array file for writing  

Parameters
----------
* `filename` :  
    File to open  
* `nrows` :  
    Number of rows  
* `ncols` :  
    Number of columns  
* `narrays` :  
    Specify a number of arrays, or -1 to add dynamically  
* `mode` :  
    File mode  
* `number_of_bits` :  
    Number of bits to use for storage. For 2-level arrays only 1 bit is needed  
";

%feature("docstring") arrayfile::arrayfile_t::~arrayfile_t "

destructor function, closes all filehandles  
";

%feature("docstring") arrayfile::arrayfile_t::createfile "

Open a new file for writing and (if opened) close the current file.  
";

%feature("docstring") arrayfile::arrayfile_t::closefile "

close the array file  
";

%feature("docstring") arrayfile::arrayfile_t::isopen "

return true if file is open  
";

%feature("docstring") arrayfile::arrayfile_t::seek "

seek to specified array position  
";

%feature("docstring") arrayfile::arrayfile_t::read_array "

read array and return index  
";

%feature("docstring") arrayfile::arrayfile_t::read_array "

read array and return index  
";

%feature("docstring") arrayfile::arrayfile_t::readnext "

read next array from the file  
";

%feature("docstring") arrayfile::arrayfile_t::readarrays "

read set of array from the file  
";

%feature("docstring") arrayfile::arrayfile_t::flush "

flush any open file pointer  
";

%feature("docstring") arrayfile::arrayfile_t::isbinary "

return true if the file has binary format  
";

%feature("docstring") arrayfile::arrayfile_t::append_arrays "

append list of arrays to the file  
";

%feature("docstring") arrayfile::arrayfile_t::append_array "

append a single array to the file  
";

%feature("docstring") arrayfile::arrayfile_t::swigcheck "

return True if code is wrapper by SWIG  
";

%feature("docstring") arrayfile::arrayfile_t::showstr "

return string describing the object  
";

%feature("docstring") arrayfile::arrayfile_t::pos "

return current position in file  
";

%feature("docstring") arrayfile::arrayfile_t::hasrandomaccess "

return true of the file format has random access mode  
";

%feature("docstring") arrayfile::arrayfile_t::updatenumbers "
";

%feature("docstring") arrayfile::arrayfile_t::finisharrayfile "
";

%feature("docstring") arrayfile::arrayfile_t::setVerbose "

set verbosity level  
";

%feature("docstring") arrayfile::arrayfile_t::getnbits "
";

%feature("docstring") arrayfile::arrayfile_t::parseModeString "

parse string to determine the file mode  
";

%feature("docstring") arrayfile::arrayfile_t::arrayNbits "

return number of bits necessary to store an array  
";

%feature("docstring") arrayfile::arrayfile_t::arrayNbits "

return number of bits necessary to store an array  
";

// File: structarraywriter__t.xml


%feature("docstring") arraywriter_t "

structure to write arrays to disk, thread safe  

C++ includes: arraytools.h
";

%feature("docstring") arraywriter_t::arraywriter_t "
";

%feature("docstring") arraywriter_t::~arraywriter_t "
";

%feature("docstring") arraywriter_t::flush "

flush all output files  
";

%feature("docstring") arraywriter_t::writeArray "

write a single array to disk  
";

%feature("docstring") arraywriter_t::writeArray "

write a list of arrays to disk  
";

%feature("docstring") arraywriter_t::initArrayFiles "

initialize the result files  
";

%feature("docstring") arraywriter_t::nArraysWritten "

return the total number arrays written to disk  
";

%feature("docstring") arraywriter_t::closeafiles "
";

// File: classCandidateGeneratorBase.xml


%feature("docstring") CandidateGeneratorBase "

Class to generate candidate extensions with caching  

We assume that the designs to be extended are run ordered, so that the caching
has maximal effect.  

The key idea used is that any valid extension of a design A with k columns is a
permutation of a valid extension of the design B obtained by taking the first l
< k columns of A. The permutations that are allowed are called the symmetry
inflations. All the j2 checks performed for the extension of B do not have to be
repeated for the permutations of this extension.  

C++ includes: conference.h
";

%feature("docstring") CandidateGeneratorBase::CandidateGeneratorBase "
";

%feature("docstring") CandidateGeneratorBase::showCandidates "

Show the candidate extensions for each column  
";

%feature("docstring") CandidateGeneratorBase::candidates "

return all candidates for the kth column  
";

// File: classCandidateGeneratorConference.xml


%feature("docstring") CandidateGeneratorConference "

Class to generate conference candidate extensions.  

C++ includes: conference.h
";

%feature("docstring") CandidateGeneratorConference::CandidateGeneratorConference "
";

%feature("docstring") CandidateGeneratorConference::generateCandidates "

Generate a list of candidate extensions for the specified design.  
";

%feature("docstring") CandidateGeneratorConference::generateCandidatesZero "

generate all candidate extensions with a zero at the specified position  
";

// File: classCandidateGeneratorDouble.xml


%feature("docstring") CandidateGeneratorDouble "

Class to generate double conference candidate extensions with caching.  

C++ includes: conference.h
";

%feature("docstring") CandidateGeneratorDouble::CandidateGeneratorDouble "
";

%feature("docstring") CandidateGeneratorDouble::generateCandidates "

Generate a list of candidate extensions for the specified design  

This method uses symmetry inflation, assumes j1=0 and j2=0. Optimal performance
is achieved when the arrays to be extended have identical first columns.  
";

// File: classCombinations.xml


%feature("docstring") Combinations "

C++ includes: mathtools.h
";

%feature("docstring") Combinations::~Combinations "
";

%feature("docstring") Combinations::number_combinations_max_n "

return max number of N that can be calculated with number_combinations  
";

%feature("docstring") Combinations::initialize_number_combinations "

initialize datastructure for number_combinations, this function is not thread
safe  
";

%feature("docstring") Combinations::number_combinations "

Return number of combinations from previously calculated results  

The results should be initialized with initialize_number_combinations  
";

// File: classgfx_1_1Compare.xml


%feature("docstring") gfx::Compare "

C++ includes: timsort.hpp
";

%feature("docstring") gfx::Compare::Compare "
";

%feature("docstring") gfx::Compare::Compare "
";

%feature("docstring") gfx::Compare::lt "
";

%feature("docstring") gfx::Compare::le "
";

%feature("docstring") gfx::Compare::gt "
";

%feature("docstring") gfx::Compare::ge "
";

%feature("docstring") gfx::Compare::less_function "
";

// File: classconference__t.xml


%feature("docstring") conference_t "

Structure representing the type of conference designs.  

C++ includes: conference.h
";

%feature("docstring") conference_t::conference_t "

Structure representing the type of conference designs  
";

%feature("docstring") conference_t::conference_t "

Structure representing the type of conference designs  

Parameters
----------
* `N` :  
    Number of rows  
* `k` :  
    Number of columns  
* `j1zero` :  
    If True then require the J1-characteristics to be zero  
";

%feature("docstring") conference_t::conference_t "
";

%feature("docstring") conference_t::idstr "
";

%feature("docstring") conference_t::create_root "

create the unique representative of the 2 column conference design in LMC0 form  
";

%feature("docstring") conference_t::create_root_three_columns "

create the unique representative of the 3 column conference design in LMC0 form  
";

%feature("docstring") conference_t::createDoubleConferenceRootArrays "

create the root arrays with 1 column for the double conference matrices  
";

%feature("docstring") conference_t::createRootArrays "

return the list of root arrays for the class of conference designs  
";

%feature("docstring") conference_t::__repr__ "

return string representation of the object  
";

// File: classconference__transformation__t.xml


%feature("docstring") conference_transformation_t "

Contains a transformation of a conference matrix.  

Contains an array transformation. The transformation consists of column
permutations, row permutations and sign switches for both the rows and columns.  

The sign switches and the permutations are not commutative. We apply the
permutations first and then the sign flips.  

C++ includes: arraytools.h
";

%feature("docstring") conference_transformation_t::conference_transformation_t "
";

%feature("docstring") conference_transformation_t::conference_transformation_t "

default constructor  
";

%feature("docstring") conference_transformation_t::conference_transformation_t "
";

%feature("docstring") conference_transformation_t::conference_transformation_t "
";

%feature("docstring") conference_transformation_t::show "

show the array transformation  
";

%feature("docstring") conference_transformation_t::isIdentity "

return true if the transformation is equal to the identity  
";

%feature("docstring") conference_transformation_t::inverse "

return the inverse transformation  
";

%feature("docstring") conference_transformation_t::reset "

return the transformation to the identity transformation  
";

%feature("docstring") conference_transformation_t::randomize "

initialize to a random transformation  
";

%feature("docstring") conference_transformation_t::randomizecolperm "

initialize with a random column permutation  
";

%feature("docstring") conference_transformation_t::randomizerowperm "

initialize with a random row permutation  
";

%feature("docstring") conference_transformation_t::randomizecolflips "

initialize with random col switches  
";

%feature("docstring") conference_transformation_t::randomizerowflips "

initialize with random row switches  
";

%feature("docstring") conference_transformation_t::apply "

apply transformation to an array_link object  
";

%feature("docstring") conference_transformation_t::setrowperm "
";

%feature("docstring") conference_transformation_t::setcolperm "
";

// File: structcounter__t.xml


%feature("docstring") counter_t "

structure to count and show number of arrays generated, the structure is thread
safe  

C++ includes: evenodd.h
";

%feature("docstring") counter_t::counter_t "
";

%feature("docstring") counter_t::addNfound "
";

%feature("docstring") counter_t::nArrays "
";

%feature("docstring") counter_t::addNumberFound "
";

%feature("docstring") counter_t::addNumberFound "
";

%feature("docstring") counter_t::clearNumberFound "
";

%feature("docstring") counter_t::showcountscompact "

show information about the number of arrays found  
";

%feature("docstring") counter_t::showcounts "

show information about the number of arrays found  
";

%feature("docstring") counter_t::showcounts "

show information about the number of arrays found  
";

// File: classDconferenceFilter.xml


%feature("docstring") DconferenceFilter "

class to filter single or double conference designs  

C++ includes: conference.h
";

%feature("docstring") DconferenceFilter::DconferenceFilter "
";

%feature("docstring") DconferenceFilter::show "

print object to stdout  
";

%feature("docstring") DconferenceFilter::filterList "

filter a list of columns using the filter method  
";

%feature("docstring") DconferenceFilter::filterListJ2last "
";

%feature("docstring") DconferenceFilter::filterListZero "

filter a list of cperms using the filterZero method  
";

%feature("docstring") DconferenceFilter::filter "

return True if the extension satisfies all checks  
";

%feature("docstring") DconferenceFilter::filterJpartial "

Filter on partial column (only last col)  

Parameters
----------
* `column` :  
    Extension column  
* `maxrow` :  
    the number of rows that are valid  
";

%feature("docstring") DconferenceFilter::filterJ "

return True if the extension satisfies all J-characteristic checks  
";

%feature("docstring") DconferenceFilter::filterJlast "

return True if the extension satisfies all J-characteristic checks for the last
columns  
";

%feature("docstring") DconferenceFilter::filterReason "

return True if the extension satisfies all checks. prints the reason for
returning True or False to stdout  
";

%feature("docstring") DconferenceFilter::filterJ3 "

return True if the candidate satisfies the J3 check  
";

%feature("docstring") DconferenceFilter::filterJ3s "

return True if the candidate satisfies the J3 check for specified pairs  
";

%feature("docstring") DconferenceFilter::filterJ3inline "

return True if the candidate satisfies the J3 check  
";

%feature("docstring") DconferenceFilter::filterSymmetry "

return True of the candidate satisfies the symmetry check  
";

%feature("docstring") DconferenceFilter::filterJ2 "

return True of the candidate extension satisfies the J2 check  
";

%feature("docstring") DconferenceFilter::filterJ2last "

return True of the candidate extension satisfies the J2 check for the last
column of the array checked against  
";

%feature("docstring") DconferenceFilter::filterZero "

return True of the candidate extension satisfies the zero check  

This means that the first entries of the extension do not contain a zero.  
";

// File: structdepth__extend__sub__t.xml


%feature("docstring") depth_extend_sub_t "

Helper structure for dynamic extension  

In this structure we keep track of pointers to valid column extensions  

C++ includes: evenodd.h
";

%feature("docstring") depth_extend_sub_t::depth_extend_sub_t "
";

%feature("docstring") depth_extend_sub_t::resize "
";

%feature("docstring") depth_extend_sub_t::n "
";

%feature("docstring") depth_extend_sub_t::updateExtensionPointers "
";

%feature("docstring") depth_extend_sub_t::initialize "

initialize the new list of extension columns  
";

%feature("docstring") depth_extend_sub_t::selectArraysZ "

select the arrays with are LMC and hence need to be written to disk  
";

%feature("docstring") depth_extend_sub_t::selectArraysXX "
";

%feature("docstring") depth_extend_sub_t::info "
";

// File: structdepth__extend__t.xml


%feature("docstring") depth_extend_t "

Helper structure for dynamic extension.  

This structure allows for writing the generated arrays to disk. It also contains
functions to print progress of the extension.  

Multiple copies of this class are made, but they all share the same counter_t
and arraywriter_t object. Also t0 and tp are shared  

C++ includes: evenodd.h
";

%feature("docstring") depth_extend_t::depth_extend_t "
";

%feature("docstring") depth_extend_t::depth_extend_t "
";

%feature("docstring") depth_extend_t::~depth_extend_t "
";

%feature("docstring") depth_extend_t::show "
";

%feature("docstring") depth_extend_t::setNarraysMax "
";

%feature("docstring") depth_extend_t::maxArrayCheck "
";

%feature("docstring") depth_extend_t::showsearchpath "
";

%feature("docstring") depth_extend_t::showprogress "

show information about the progress of the loop  
";

%feature("docstring") depth_extend_t::info "
";

%feature("docstring") depth_extend_t::setposition "

set the position in the dextend structure  
";

%feature("docstring") depth_extend_t::setpositionGEC "

set the position in the dextend structure  
";

// File: structdepth__extensions__storage__t.xml


%feature("docstring") depth_extensions_storage_t "

Helper structure for the even-odd depth extension.  

C++ includes: evenodd.h
";

%feature("docstring") depth_extensions_storage_t::resize "
";

%feature("docstring") depth_extensions_storage_t::set "
";

// File: structdepth__path__t.xml


%feature("docstring") depth_path_t "

structure containing current position in search tree  

C++ includes: evenodd.h
";

%feature("docstring") depth_path_t::depth_path_t "
";

%feature("docstring") depth_path_t::updatePositionGEC "
";

%feature("docstring") depth_path_t::updatePosition "
";

%feature("docstring") depth_path_t::show "
";

%feature("docstring") depth_path_t::init "
";

// File: structdextend__t.xml


%feature("docstring") dextend_t "

Structure for dynamic extension of arrays based on D-efficiencies.  

C++ includes: extend.h
";

%feature("docstring") dextend_t::dextend_t "
";

%feature("docstring") dextend_t::resize "
";

%feature("docstring") dextend_t::DefficiencyFilter "

perform filtering using D-efficiency  
";

%feature("docstring") dextend_t::filterArrays "

filter the arrays based on values in filter  
";

// File: structDoptimReturn.xml


%feature("docstring") DoptimReturn "

Structure containing results of the Doptimize function  

C++ includes: Deff.h
";

// File: structdyndata__t.xml


%feature("docstring") dyndata_t "

Contains dynamic data of an array.  

The dynamic data are used in the inner loops of the LMC algorithm. In particular
they keep track of the current row ordering and column permutation. By not
applying these transformations to the array we can save calculation time.  

We try to prevent copying the object, so it is re-used at different levels in
the algorithm.  

*   N: static
    -   col: changes at each column level  
*   rowsort: changes at each column level, used mainly in non-root stage  
*   colperm: changes at all levels  

    See also: arraydata_t  

C++ includes: lmc.h
";

%feature("docstring") dyndata_t::dyndata_t "
";

%feature("docstring") dyndata_t::dyndata_t "
";

%feature("docstring") dyndata_t::dyndata_t "
";

%feature("docstring") dyndata_t::~dyndata_t "
";

%feature("docstring") dyndata_t::show "
";

%feature("docstring") dyndata_t::reset "
";

%feature("docstring") dyndata_t::setColperm "
";

%feature("docstring") dyndata_t::setColperm "
";

%feature("docstring") dyndata_t::setColperm "
";

%feature("docstring") dyndata_t::getRowperm "

get lightweight row permutation  
";

%feature("docstring") dyndata_t::getRowperm "

get row permutation  
";

%feature("docstring") dyndata_t::getRowperm "

return lightweight row permutation  
";

%feature("docstring") dyndata_t::getColperm "

return column permutation  
";

%feature("docstring") dyndata_t::getColperm "

set column permutation  
";

%feature("docstring") dyndata_t::allocate_rowsortl "

allocate lightweight rowsort structure  
";

%feature("docstring") dyndata_t::deleterowsortl "
";

%feature("docstring") dyndata_t::initrowsortl "

initialize rowsortl from rowsort  
";

%feature("docstring") dyndata_t::rowsortl2rowsort "

copy rowsortl variable to rowsrt  
";

%feature("docstring") dyndata_t::copydata "
";

// File: structextend__data__t.xml


%feature("docstring") extend_data_t "

Contains static data for the extend loop.  

C++ includes: strength.h
";

%feature("docstring") extend_data_t::extend_data_t "
";

%feature("docstring") extend_data_t::~extend_data_t "
";

%feature("docstring") extend_data_t::init_frequencies "

Initialize the table of t-tuple frequencies.  
";

// File: classindexsort.xml


%feature("docstring") indexsort "

Class to sort data without moving the data in memory.  

The data is sorted by using a list of indices. A stable sort is being used.  

C++ includes: mathtools.h
";

%feature("docstring") indexsort::indexsort "
";

%feature("docstring") indexsort::indexsort "

Constructor for deque class.  
";

%feature("docstring") indexsort::indexsort "

Constructor for vector class.  
";

%feature("docstring") indexsort::init "

initialize sorting structure with specified values  
";

%feature("docstring") indexsort::init "

initialize sorting structure with specified values  
";

%feature("docstring") indexsort::sort "

sort values and store the indices  
";

%feature("docstring") indexsort::sort "

sort values and store the indices  
";

%feature("docstring") indexsort::sort "

sort values and store the indices  
";

%feature("docstring") indexsort::sortdescending "

sort values and store the indices  
";

%feature("docstring") indexsort::show "
";

%feature("docstring") indexsort::sorted "

return array sorted using the order from the indexsort structure  
";

%feature("docstring") indexsort::sorted "

return array sorted using the order from the indexsort structure  
";

%feature("docstring") indexsort::issorted "

Returns true of the data is sorted ascending.  
";

%feature("docstring") indexsort::issorteddescending "

Returns true of the data is sorted descending.  
";

// File: classInfInt.xml


%feature("docstring") InfInt "

C++ includes: InfInt.h
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::InfInt "
";

%feature("docstring") InfInt::intSqrt "
";

%feature("docstring") InfInt::digitAt "
";

%feature("docstring") InfInt::numberOfDigits "
";

%feature("docstring") InfInt::size "
";

%feature("docstring") InfInt::toString "
";

%feature("docstring") InfInt::toInt "
";

%feature("docstring") InfInt::toLong "
";

%feature("docstring") InfInt::toUnsignedInt "
";

%feature("docstring") InfInt::toUnsignedLong "
";

// File: classJcounter.xml


%feature("docstring") Jcounter "

object to hold counts of maximum J_k-values  

C++ includes: evenodd.h
";

%feature("docstring") Jcounter::Jcounter "
";

%feature("docstring") Jcounter::Jcounter "
";

%feature("docstring") Jcounter::validData "
";

%feature("docstring") Jcounter::hasColumn "

return true if specified column is in the data  
";

%feature("docstring") Jcounter::isOpen "
";

%feature("docstring") Jcounter::showPerformance "
";

%feature("docstring") Jcounter::narrays "
";

%feature("docstring") Jcounter::show "

show statistics of the object  
";

%feature("docstring") Jcounter::maxCols "
";

%feature("docstring") Jcounter::getCount "
";

%feature("docstring") Jcounter::getTotalsJvalue "
";

%feature("docstring") Jcounter::getTotals "
";

%feature("docstring") Jcounter::showcompact "

show statistics of the object  
";

%feature("docstring") Jcounter::addArrays "

add list of arrays to object  
";

%feature("docstring") Jcounter::addArray "

add single array to statistics object  
";

// File: structjindex__t.xml


%feature("docstring") jindex_t "

helper class for indexing statistics of designs  

The index consists of the number of columns and the value for the
J-characteristic  

C++ includes: evenodd.h
";

%feature("docstring") jindex_t::jindex_t "
";

%feature("docstring") jindex_t::toString "
";

// File: classjstruct__t.xml


%feature("docstring") jstruct_t "

struct to hold data of an array, e.g. J-characteristic, rank  

See papers: Minimum G2-aberration properties of two-level foldover designs,
Butler, 2004 Design Selection and Classification for Hadamard Matrices Using
Generalized Minimum Aberration Criteria, Deng and Tang  

C++ includes: arraytools.h
";

%feature("docstring") jstruct_t::jstruct_t "

Create an object to calculate J-characteristics.  
";

%feature("docstring") jstruct_t::jstruct_t "

Create an object to calculate J-characteristics.  
";

%feature("docstring") jstruct_t::jstruct_t "
";

%feature("docstring") jstruct_t::jstruct_t "
";

%feature("docstring") jstruct_t::~jstruct_t "
";

%feature("docstring") jstruct_t::maxJ "

calculate maximum J value  
";

%feature("docstring") jstruct_t::Fval "

calculate possible values in F vector  
";

%feature("docstring") jstruct_t::calculateF "

calculate histogram of J values for a 2-level array  
";

%feature("docstring") jstruct_t::calculateAberration "

Calculate aberration value  

This is equal to the sum of the squares of all Jk values, divided by the number
of rows squared.  
";

%feature("docstring") jstruct_t::show "

Show contents of structure.  
";

%feature("docstring") jstruct_t::showdata "
";

%feature("docstring") jstruct_t::showstr "
";

%feature("docstring") jstruct_t::allzero "

return 1 if all J values are zero, otherwise return 0  
";

// File: classjstructbase__t.xml


%feature("docstring") jstructbase_t "

struct to hold data of an array, e.g. J-characteristic. Abstract base class  

C++ includes: arraytools.h
";

%feature("docstring") jstructbase_t::maxJ "

calculate maximum J value  
";

%feature("docstring") jstructbase_t::Jvalues "

calculate possible values in F vector  
";

%feature("docstring") jstructbase_t::calculateF "

Calculate histogram of J values  

       The histogram bins are given by the values of @ref Jvalues

       \\returns Histogram of J values
  
";

%feature("docstring") jstructbase_t::calc "

Calculate the J-values for a given array.  
";

%feature("docstring") jstructbase_t::show "

Show contents of structure.  
";

%feature("docstring") jstructbase_t::showdata "
";

%feature("docstring") jstructbase_t::showstr "
";

%feature("docstring") jstructbase_t::allzero "

return 1 if all vals are zero  
";

// File: classjstructconference__t.xml


%feature("docstring") jstructconference_t "

Calculate J-characteristics of conference designs  

C++ includes: arraytools.h
";

%feature("docstring") jstructconference_t::jstructconference_t "

Create structure to calculate J-characteristics of conference designs  

Parameters
----------
* `N` :  
    Number of rows  
* `jj` :  
    Number of columns to use for the Jk-characteristics  
";

%feature("docstring") jstructconference_t::jstructconference_t "

Calculate J-characteristics of a conference design  

Parameters
----------
* `array` :  
    Array to calculate the J-characteristics for  
* `jj` :  
    Number of columns to use for the Jk-characteristics  
";

// File: classlarray.xml


%feature("docstring") larray "

lightweight array class  

C++ includes: mathtools.h
";

%feature("docstring") larray::larray "

Create unallocated array.  
";

%feature("docstring") larray::larray "
";

%feature("docstring") larray::larray "
";

%feature("docstring") larray::larray "
";

%feature("docstring") larray::begin "
";

%feature("docstring") larray::~larray "
";

%feature("docstring") larray::resize "
";

%feature("docstring") larray::size "
";

%feature("docstring") larray::at "
";

%feature("docstring") larray::addelement "

add constant value to the elements of the array  
";

// File: structLMCreduction__helper__t.xml


%feature("docstring") LMCreduction_helper_t "

Contains structures used by the LMC reduction or LMC check.  

Part of the allocations is for structures that are constant and are re-used each
time an LMC calculation is performed. Some other structures are temporary
buffers that are written to all the time.  

C++ includes: lmc.h
";

%feature("docstring") LMCreduction_helper_t::LMCreduction_helper_t "
";

%feature("docstring") LMCreduction_helper_t::~LMCreduction_helper_t "
";

%feature("docstring") LMCreduction_helper_t::show "
";

%feature("docstring") LMCreduction_helper_t::init "
";

%feature("docstring") LMCreduction_helper_t::freeall "
";

%feature("docstring") LMCreduction_helper_t::update "

update structure with new design specification  
";

%feature("docstring") LMCreduction_helper_t::needUpdate "
";

%feature("docstring") LMCreduction_helper_t::init_root_stage "
";

%feature("docstring") LMCreduction_helper_t::init_nonroot_stage "
";

%feature("docstring") LMCreduction_helper_t::init_rootrowperms "

Static initialization of root row permutations.  
";

%feature("docstring") LMCreduction_helper_t::init_rootrowperms_full "

Static initialization of root row permutations (full group)  
";

// File: structLMCreduction__t.xml


%feature("docstring") LMCreduction_t "

Class to describe an LMC reduction.  

The most important variable is the transformation itself, contained in
transformation. The state contains information about how the reduction was
performed.  

C++ includes: lmc.h
";

%feature("docstring") LMCreduction_t::LMCreduction_t "
";

%feature("docstring") LMCreduction_t::LMCreduction_t "

copy constructor  
";

%feature("docstring") LMCreduction_t::~LMCreduction_t "
";

%feature("docstring") LMCreduction_t::getArray "

Assignment operator.  
";

%feature("docstring") LMCreduction_t::setArray "
";

%feature("docstring") LMCreduction_t::setArray "
";

%feature("docstring") LMCreduction_t::updateSDpointer "

update the pointer to the symmetry data based on the specified array  
";

%feature("docstring") LMCreduction_t::releaseStatic "

release internal LMCreduction_helper_t object  
";

%feature("docstring") LMCreduction_t::initStatic "

acquire a reference to a LMCreduction_helper_t object  
";

%feature("docstring") LMCreduction_t::getReferenceReductionHelper "

return a reference to a object with LMC reduction data  
";

%feature("docstring") LMCreduction_t::reset "

reset the reduction: clears the symmetries and sets the transformation to zero  
";

%feature("docstring") LMCreduction_t::show "
";

%feature("docstring") LMCreduction_t::__repr__ "
";

%feature("docstring") LMCreduction_t::updateFromLoop "

called whenever we find a reduction  
";

%feature("docstring") LMCreduction_t::updateTransformation "
";

%feature("docstring") LMCreduction_t::updateLastCol "
";

// File: structmvalue__t.xml


%feature("docstring") mvalue_t "

Multi-value type.  

This object represents a multi-valued object. The objects are ordered using
lexicographic ordering.  

C++ includes: mathtools.h
";

%feature("docstring") mvalue_t::mvalue_t "

Create multi-valued object  

The object consists of a vector of elements.  
";

%feature("docstring") mvalue_t::mvalue_t "

Create multi-valued object  

The object consists of a vector of elements.  

Parameters
----------
* `element` :  
    Single element to add to the vector  
* `dd` :  
    Ordering to use  
";

%feature("docstring") mvalue_t::mvalue_t "

Create multi-valued object  

The object consists of a vector of elements.  

Parameters
----------
* `elements` :  
    Vector to use for initalization of the object  
* `dd` :  
    Ordering to use  
";

%feature("docstring") mvalue_t::mvalue_t "

Create multi-valued object  

The object consists of a vector of elements.  

Parameters
----------
* `elements` :  
    Vector to use for initalization of the object  
* `dd` :  
    Ordering to use  
";

%feature("docstring") mvalue_t::~mvalue_t "
";

%feature("docstring") mvalue_t::raw_values "

Return vector with the raw values in this object.  
";

%feature("docstring") mvalue_t::size "
";

%feature("docstring") mvalue_t::show_integer "

Show the object on stdout by casting to integer type objects.  
";

%feature("docstring") mvalue_t::string_representation "

return a string representation of the object  
";

// File: classOAextend.xml


%feature("docstring") OAextend "

Options for the extend code.  

class containing parameters of the extension and LMC algorithm  

C++ includes: extend.h
";

%feature("docstring") OAextend::OAextend "

Options for the extension algorithm  
";

%feature("docstring") OAextend::OAextend "

Options for the extension algorithm  
";

%feature("docstring") OAextend::OAextend "

Options for the extension algorithm  

The algorithm is automatically determined from the specified arrayclass.  
";

%feature("docstring") OAextend::setAlgorithm "

Set the algorithm to use for LMC checks.  
";

%feature("docstring") OAextend::setAlgorithmAuto "

Set the algorithm automatically.  
";

%feature("docstring") OAextend::getAlgorithm "

Return algorithm used.  
";

%feature("docstring") OAextend::getAlgorithmName "

Return algorithm used (as string)  
";

%feature("docstring") OAextend::updateArraydata "

update the options structuer with the specified class of designs  
";

%feature("docstring") OAextend::info "

print configuration to stdout  
";

%feature("docstring") OAextend::__repr__ "
";

%feature("docstring") OAextend::getPreferredAlgorithm "

return preferred extension algorithm  

       \\param arrayclass Class of designs to extend
       \\param verbose Verbosity level
  
";

// File: classobject__pool.xml


%feature("docstring") object_pool "

Class to make a pool of objects that can be re-used.  

C++ includes: mathtools.h
";

%feature("docstring") object_pool::object_pool "

Create a pool of objects that can be re-used.  
";

%feature("docstring") object_pool::reset "

assume the pool is filled with pointers: remove them all and call the destructor  
";

%feature("docstring") object_pool::begin "
";

%feature("docstring") object_pool::begin "
";

%feature("docstring") object_pool::end "
";

%feature("docstring") object_pool::end "
";

%feature("docstring") object_pool::at "
";

%feature("docstring") object_pool::at "
";

%feature("docstring") object_pool::size "
";

%feature("docstring") object_pool::New "
";

%feature("docstring") object_pool::Delete "
";

// File: classPareto.xml


%feature("docstring") Pareto "

Class to the calculate Pareto optimal elements.  

The class is templated by the type of values to be compared and an index type.
The index type is used to index the elements.  

For elements added to the Pareto structure larger is better.  

C++ includes: pareto.h
";

%feature("docstring") Pareto::Pareto "

constructor  
";

%feature("docstring") Pareto::~Pareto "
";

%feature("docstring") Pareto::number "

return the total number of Pareto optimal values  
";

%feature("docstring") Pareto::numberindices "

return the toal number Pareto optimal objects  
";

%feature("docstring") Pareto::__repr__ "
";

%feature("docstring") Pareto::show "

show the current set of Pareto optimal elements  
";

%feature("docstring") Pareto::allindicesdeque "

return all indices of the Pareto optimal elements as a std::deque  
";

%feature("docstring") Pareto::allindices "

return all indices of the Pareto optimal elements  
";

%feature("docstring") Pareto::allvalues "

return all Paretop optimal elements  
";

%feature("docstring") Pareto::addvalue "

add a new element  
";

%feature("docstring") Pareto::showvalue "
";

// File: structpareto__element.xml


%feature("docstring") pareto_element "

helper class for the Pareto class to hold elements  

C++ includes: pareto.h
";

%feature("docstring") pareto_element::dominates "

return true of the argument element dominates this value  
";

%feature("docstring") pareto_element::isdominated "

return true of the argument element is dominated by this value  
";

%feature("docstring") pareto_element::equal "

return true of the argument element is equal to this element  
";

// File: classrankStructure.xml


%feature("docstring") rankStructure "

Structure to efficiently calculate the rank of the second order interaction
matrix of many arrays  

The efficiency is obtained if the arrays share a common subarray. The theory is
described in \"Efficient rank calculation for matrices with a common
submatrix\", Eendebak, 2016  

C++ includes: arrayproperties.h
";

%feature("docstring") rankStructure::rankStructure "

constructor  
";

%feature("docstring") rankStructure::rankStructure "

constructor  
";

%feature("docstring") rankStructure::info "

print information about the rank structure  
";

%feature("docstring") rankStructure::updateStructure "

update the structure cache with a new array  
";

%feature("docstring") rankStructure::rankdirect "

calculate the rank of an array directly, uses special threshold  
";

%feature("docstring") rankStructure::rankxfdirect "

calculate the rank of the second order interaction matrix of an array directly  
";

%feature("docstring") rankStructure::rankxf "

calculate the rank of the second order interaction matrix of an array using the
cache system  
";

// File: classrowsorter__t.xml


%feature("docstring") rowsorter_t "

Structure to sort rows of arrays.  

C++ includes: lmc.h
";

%feature("docstring") rowsorter_t::rowsorter_t "
";

%feature("docstring") rowsorter_t::~rowsorter_t "
";

// File: structgfx_1_1TimSort_1_1run.xml

// File: classsort__indices.xml


%feature("docstring") sort_indices "

Helper class.  

C++ includes: mathtools.h
";

%feature("docstring") sort_indices::sort_indices "
";

// File: classsort__indices__container.xml


%feature("docstring") sort_indices_container "

Helper class.  

C++ includes: mathtools.h
";

%feature("docstring") sort_indices_container::sort_indices_container "
";

// File: classsort__indices__deque.xml


%feature("docstring") sort_indices_deque "

Helper class.  

C++ includes: mathtools.h
";

%feature("docstring") sort_indices_deque::sort_indices_deque "
";

// File: classsort__indices__vector.xml


%feature("docstring") sort_indices_vector "

Helper class.  

C++ includes: mathtools.h
";

%feature("docstring") sort_indices_vector::sort_indices_vector "
";

// File: structsymmdata.xml


%feature("docstring") symmdata "

structure containing data related to symmetries of arrays  

C++ includes: arraytools.h
";

%feature("docstring") symmdata::symmdata "
";

%feature("docstring") symmdata::show "
";

%feature("docstring") symmdata::checkIdx "

list with indices set to check for symmetry reductions  
";

// File: classsymmetry__group.xml


%feature("docstring") symmetry_group "

Class to describe the symmetry group of a list of elements.  

The class assumes the list is sorted. The symmetry group is then a direct
product of full permutation groups.  

We do not implement this using templates because we want to export to Python.  

C++ includes: mathtools.h
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "
";

%feature("docstring") symmetry_group::symmetry_group "

default constructor  
";

%feature("docstring") symmetry_group::permsize "

Return size of the group of all permutations respecting the symmetry  

The return type can overflow quickly. For larger group sizes use permsize_large  
";

%feature("docstring") symmetry_group::permsize_large "

return size of the group of all permutations respecting the symmetry  
";

%feature("docstring") symmetry_group::checkIndices "

list with indices set to check for symmetry reductions  
";

%feature("docstring") symmetry_group::__repr__ "

representation function (for python interface)  
";

%feature("docstring") symmetry_group::show "

show the symmetry group  
";

// File: classsymmetry__group__walker.xml


%feature("docstring") symmetry_group_walker "

Class to walk over all elements of a symmetry group  

The elements are generated by walking over all product permutations.  

C++ includes: mathtools.h
";

%feature("docstring") symmetry_group_walker::symmetry_group_walker "
";

%feature("docstring") symmetry_group_walker::show "

show all elements in the symmetry group  
";

%feature("docstring") symmetry_group_walker::next "

go to next element of the symmetry group  
";

%feature("docstring") symmetry_group_walker::fullperm "

return the full permutation corresponding to the current permutation in the
walker  
";

// File: classgfx_1_1TimSort.xml


%feature("docstring") gfx::TimSort "

C++ includes: timsort.hpp
";

// File: namespacearrayfile.xml

// File: namespacegfx.xml

%feature("docstring") gfx::timsort "

Same as std::stable_sort(first, last).  
";

%feature("docstring") gfx::timsort "

Same as std::stable_sort(first, last, c).  
";

// File: namespacenauty.xml

%feature("docstring") nauty::reduceNauty "

Reduce a colored graph to Nauty minimal form  

The transformation returned is from the normal form to the specified graph.  

Parameters
----------
* `graph` :  
    Graph in incidence matrix form  
* `colors` :  
    Colors of the graph nodes  
* `verbose` :  
    Verbosity level  

Returns
-------
Relabelling of the graph vertices  
";

// File: namespacestd.xml

// File: anyoption_8h.xml

// File: arrayproperties_8h.xml

%feature("docstring") DAEefficiencyWithSVD "

Calculate D-efficiency and VIF-efficiency and E-efficiency values using SVD.  
";

%feature("docstring") array2rank_Deff_Beff "

Calculate the rank of the second order interaction matrix of an orthogonal array  

The model is the intercept, main effects and interaction effects The rank,
D-efficiency, VIF-efficiency and E-efficiency are appended to the second
argument  

The return vector is filled with the rank, Defficiency, VIF efficiency and
Eefficiency  
";

%feature("docstring") Defficiency "

Calculate D-efficiency for a 2-level array using symmetric eigenvalue
decomposition.  
";

%feature("docstring") Defficiencies "

Calculate efficiencies for an array  

Parameters
----------
* `array` :  
    Array to use in calculation  
* `arrayclass` :  
    Specification of the array class  
* `verbose` :  
    Verbosity level  
* `addDs0` :  
    If True, then add the Ds0-efficiency to the output  

Returns
-------
Vector with the calculate D-efficiency, the main effect robustness (or Ds-
optimality) and D1-efficiency for an orthogonal array  
";

%feature("docstring") VIFefficiency "

Calculate VIF-efficiency of matrix.  
";

%feature("docstring") Aefficiency "

Calculate A-efficiency of matrix.  
";

%feature("docstring") Eefficiency "

Calculate E-efficiency of matrix (1 over the VIF-efficiency)  
";

%feature("docstring") Aefficiencies "

calculate various A-efficiencies  
";

%feature("docstring") projDeff "

Calculate D-efficiencies for all projection designs  

Parameters
----------
* `array` :  
    Design to calculate D-efficiencies for  
* `number_of_factors` :  
    Number of factors into which to project  
* `verbose` :  
    Verbosity level  

Returns
-------
Vector with calculated D-efficiencies  
";

%feature("docstring") PECsequence "

Calculate the projection estimation capacity sequence for a design  

Parameters
----------
* `array` :  
    Input array  
* `verbose` :  
    Verbosity level  

Returns
-------
Vector with the caculated PEC sequence  

The PECk of a design is the fraction of estimable second-order models in k
factors. The vector (PEC1, PEC2, ..., ) is called the projection estimation
capacity sequence. See \"Ranking Non-regular Designs\", J.L. Loeppky, 2004.  
";

%feature("docstring") PICsequence "

Calculate the projection information capacity sequence for a design.  

Parameters
----------
* `array` :  
    Input array  
* `verbose` :  
    Verbosity level  

Returns
-------
Vector with the caculated PIC sequence  

The PICk of a design is the average D-efficiency of estimable second-order
models in k factors. The vector (PIC1, PIC2, ..., ) is called the PIC sequence.  
";

%feature("docstring") distance_distribution "

Return the distance distribution of a design  

The distance distribution is described in \"Generalized minimum aberration for
asymmetrical fractional factorial designs\", Wu and Xu, 2001  
";

%feature("docstring") Jcharacteristics "

Calculate Jk-characteristics of a matrix  

The calcualted Jk-values are signed.  

Parameters
----------
* `array` :  
    Array to calculate Jk-characteristics for  
* `number_of_columns` :  
    Number of columns  
* `verbose` :  
    Verbosity level  

Returns
-------
Vector with calculated Jk-characteristics  
";

%feature("docstring") GWLP "

Calculate GWLP (generalized wordlength pattern)  

The method used for calculation is from Xu and Wu (2001), \"Generalized minimum
aberration for asymmetrical
fractional factorial desings\". For non-symmetric arrays see \"Algorithmic
Construction of Efficient Fractional Factorial Designs With Large Run
Sizes\", Xu, Technometrics, 2009.  

Parameters
----------
* `array` :  
    Array to calculate the GWLP value for  
* `verbose` :  
    Verbosity level  
* `truncate` :  
    If True then round values near zero to solve double precision errors  

Returns
-------
Vector with calculated generalized wordlength pattern  

A more detailed description of the generalized wordlength pattern can also be
found in the documentation at https://oapackage.readthedocs.io/.  
";

%feature("docstring") GWLPmixed "

Calculate GWLP (generalized wordlength pattern) for mixed-level arrays.  

The method used for calculation is from \"Algorithmic Construction of Efficient
Fractional Factorial Designs With Large Run
Sizes\", Xu, Technometrics, 2009.  

Parameters
----------
* `array` :  
    Array to calculate the GWLP value for  
* `verbose` :  
    Verbosity level  
* `truncate` :  
    If True then round values near zero to solve double precision errors  

Returns
-------
Vector with calculated generalized wordlength pattern  
";

%feature("docstring") projectionGWLPs "

calculate delete-one-factor GWLP (generalized wordlength pattern) projections  
";

%feature("docstring") sortGWLP "

sort a list of GWLP values and return the sorted list  
";

%feature("docstring") CL2discrepancy "

Calculate centered L2-discrepancy of a design  

The method is from \"A connection between uniformity and aberration in regular
fractions of two-level factorials\", Fang and Mukerjee, 2000  
";

%feature("docstring") array2secondorder "

Calculate second order interaction model for 2-level array  

Parameters
----------
* `array` :  
    Array to calculate second order interaction model from  

Returns
-------
Array interaction effects  
";

%feature("docstring") array2xf "

calculate second order interaction model for 2-level array  

Parameters
----------
* `array` :  
    Array to calculate second order interaction model from  

Returns
-------
Array with intercept, main effects and interaction effects  
";

%feature("docstring") conference_design2modelmatrix "

Calculate model matrix for a conference design  

Parameters
----------
* `conference_design` :  
    Conference design  
* `mode` :  
    Can be 'm' for main effects, 'i' for interaction effects or 'q' for
    quadratic effects  
* `verbose` :  
    Verbosity level  

Returns
-------
Calculated model matrix  
";

%feature("docstring") array2modelmatrix "

Convert orthogonal array or conference design to model matrix  

The model matrix consists of the intercept, main effects and (optionally) the
interaction effects and quadratic effects. The order in the interaction effects
is (c1, c2)=(0,0), (1,0), (2,0), (2,1), ... with c2<c1 for columns c1, c2. The
size of the model matrix calculated by this function is given by
array2modelmatrix_sizes.  

Parameters
----------
* `array` :  
    Orthogonal array or conference design  
* `mode` :  
    Type of model matrix to calculate. Can be 'm' for main effects, 'i' for
    interaction effects or 'q' for quadratic effects  
* `verbose` :  
    Verbosity level  

Returns
-------
Calculated model matrix  

For conference designs the method conference_design2modelmatrix is used. For
orthogonal array the calculated is performed with array2eigenModelMatrixMixed.  
";

%feature("docstring") array2modelmatrix_sizes "

Return the sizes of the model matrices calculated  

Parameters
----------
* `array` :  
    Orthogonal array or conference designs  

Returns
-------
List with the sizes of the model matrix for: only intercept; intercept, main;
intercept, main, and iteraction terms, intercept, main and full second order  
";

%feature("docstring") array2xfeigen "

calculate second order interaction model for 2-level array  

Parameters
----------
* `array` :  
    Array to calculate second order interaction model from  

Returns
-------
Array with intercept, main effects and interaction effects  
";

%feature("docstring") arrayrankFullPivQR "

return rank of an array based on Eigen::FullPivHouseholderQR  
";

%feature("docstring") arrayrankColPivQR "

return rank of an array based on Eigen::ColPivHouseholderQR  
";

%feature("docstring") arrayrankFullPivLU "

return rank of an array based on Eigen::FullPivLU  
";

%feature("docstring") arrayrankSVD "

return rank of an array based on Eigen::JacobiSVD  
";

%feature("docstring") arrayrank "

calculate the rank of an array  
";

%feature("docstring") arrayrankInfo "

Return rank of an array. Information about the different methods for rank
calculation is printed to stdout.  
";

%feature("docstring") arrayrankInfo "

Return rank of an array. Information about the different methods for rank
calculation is printed to stdout.  
";

%feature("docstring") arraylink2eigen "

convert array_link to Eigen matrix  
";

%feature("docstring") conditionNumber "

Return the condition number of a matrix.  
";

%feature("docstring") calculateParetoEvenOdd "

Calculate the Pareto optimal arrays from a list of array files  

Pareto optimality is calculated according to (rank; A3,A4; F4)  
";

%feature("docstring") parsePareto "
";

%feature("docstring") A3A4 "

calculate A3 and A4 value for array  

Parameters
----------
* `al` :  
    Array for which to calculate A3 and A4  

Returns
-------
Object with A3 and A4  
";

%feature("docstring") F4 "

calculate F4 value for 2-level array  
";

%feature("docstring") calculateArrayParetoRankFA "

Calculate properties of an array and create a Pareto element  

The values calculated are:  

1) Rank (higher is better) 2) A3, A4 (lower is better) 3) F4 (lower is better,
sum of elements is constant)  

Valid for 2-level arrays of strength at least 3  
";

%feature("docstring") addJmax "

add Jmax criterium to Pareto set  
";

%feature("docstring") calculateArrayParetoJ5 "

Calculate Pareto element with J5 criterium.  
";

%feature("docstring") parseArrayPareto "

Add array to list of Pareto optimal arrays  

The values to be optimized are:  

1) Rank (higher is better) 2) A3, A4 (lower is better) 3) F4 (lower is better,
sum of elements is constant)  
";

%feature("docstring") Cvalue2Dvalue "

convert C value to D-efficiency value  
";

%feature("docstring") Dvalue2Cvalue "

convert D-efficiency value to C value  
";

// File: arraytools_8h.xml

%feature("docstring") arrayfile::throw_runtime_exception "
";

%feature("docstring") arrayfile::eigenInfo "

Print information about an Eigen matrix  

Parameters
----------
* `m` :  
    Matrix about which to print information  
* `str` :  
    String to prepend in output  
* `verbose` :  
    Verbosity level  
";

%feature("docstring") arrayfile::print_eigen_matrix "

Print Eigen matrix to stdout  
";

%feature("docstring") arrayfile::eigen2numpyHelper "
";

%feature("docstring") arrayfile::sizeof_array_t "

return size in bytes of array_t type  
";

%feature("docstring") arrayfile::possible_F_values "

possible values for J-values of 2-level design  
";

%feature("docstring") arrayfile::file_exists "

return true if the specified file exists  
";

%feature("docstring") arrayfile::file_exists "

return true if the specified file exists  
";

%feature("docstring") arrayfile::oa_file_exists "

return true if the specified oa file exists  
";

%feature("docstring") arrayfile::oa_file_exists "

return true if the specified oa file exists  
";

%feature("docstring") arrayfile::readConfigFile "

Read array configuration from file.  
";

%feature("docstring") arrayfile::printfstring "

Function similar to printf returning C++ style string.  

Parameters
----------
* `message` :  

Returns
-------  
";

%feature("docstring") arrayfile::copy_array "

Make a copy of an array.  
";

%feature("docstring") arrayfile::destroy_array "

Delete an array.  

Parameters
----------
* `array` :  

Returns
-------  
";

%feature("docstring") arrayfile::create_array "

Create an array.  

Parameters
----------
* `nrows` :  
    Number of rows  
* `ncols` :  
    Number of columns  

Returns
-------  
";

%feature("docstring") arrayfile::create_array "

Create an array from an arraydata_t structure.  
";

%feature("docstring") arrayfile::clone_array "

Clone an array.  
";

%feature("docstring") arrayfile::exampleArray "

Return example array  

Parameters
----------
* `idx` :  
    Index of example array to return  
* `verbose` :  
    If True, then print information about the array to stdout  
";

%feature("docstring") arrayfile::Jcharacteristics_conference "

Calculate Jk-characteristics for a conference design  

Parameters
----------
* `array` :  
    Conference design  
* `number_of_columns` :  
    Specifies the number of columns to use  
* `verbose` :  
    Verbosity level  

Returns
-------
A vector of calculated inner products between all combinations of k columns.  
";

%feature("docstring") arrayfile::hstack "

concatenate 2 arrays in vertical direction  

concatenate 2 arrays in horizontal direction  
";

%feature("docstring") arrayfile::hstack "

concatenate array and conference_column  
";

%feature("docstring") arrayfile::hstacklastcol "

concatenate the last column of array B to array A  
";

%feature("docstring") arrayfile::vstack "

concatenate two columns  
";

%feature("docstring") arrayfile::perform_column_permutation "

perform column permutation for an array  
";

%feature("docstring") arrayfile::perform_row_permutation "

perform row permutation for an array  
";

%feature("docstring") arrayfile::arraylink2arraydata "

create arraydata_t structure from array  

Parameters
----------
* `array` :  
    Array to use as input specifiction for array class  
* `extracols` :  
    Number of extra columns to add to the number of columns of the array  
* `strength` :  
    Strength to set in the array class. If -1, then use the strength of the
    array  
";

%feature("docstring") arrayfile::addConstant "

add a constant value to all arrays in a list  
";

%feature("docstring") arrayfile::getJcounts "

Return number of arrays with j_{2n+1}=0 for number_of_arrays<m  
";

%feature("docstring") arrayfile::create_root "

set first columns of an array to root form  
";

%feature("docstring") arrayfile::create_root "

Creates the root of an orthogonal array. The root is appended to the list of
arrays.  
";

%feature("docstring") arrayfile::array_diff "

Compare 2 arrays and return position of first difference.  
";

%feature("docstring") arrayfile::fastJupdate "

helper function to calculate J-values  
";

%feature("docstring") arrayfile::jvalue "

Calculate J-value for a 2-level array  
";

%feature("docstring") arrayfile::jvaluefast "

Calculate J-value for a column combination of a 2-level array  

We assume the array has values 0 and 1. No boundary checks are performed.  
";

%feature("docstring") arrayfile::analyseArrays "

Analyse a list of arrays.  
";

%feature("docstring") arrayfile::showArrayList "

print a list of arrays to stdout  
";

%feature("docstring") arrayfile::nArrays "

return number of arrays in an array file  
";

%feature("docstring") arrayfile::arrayfileinfo "

return information about file with arrays  

Parameters
----------
* `filename` :  
    Filename of array file  
* `number_of_arrays` :  
    Variable is set with number of arrays  
* `number_of_rows` :  
    Variable is set with number of rows  
* `number_of_columns` :  
    Variable is set with number of columns  
";

%feature("docstring") arrayfile::readarrayfile "

Read all arrays in a file  

Parameters
----------
* `fname` :  
    Filename to read from  
* `verbose` :  
    Verbosity level  
* `setcols` :  
    Pointer to return number of columns from array file  

Returns
-------
List of arrays  
";

%feature("docstring") arrayfile::readarrayfile "

Read all arrays in a file and append then to an array list  

Parameters
----------
* `filename` :  
    Filename to read from  
* `arraylist` :  
    Pointer to list of arrays  
* `verbose` :  
    Verbosity level  
* `setcols` :  
    Reference that is set with the number of columns from the file  
* `setrows` :  
    Reference that is set with the number of rows from the file  
* `setbits` :  
    Reference that is set with the number of bits from the file  

Returns
-------  
";

%feature("docstring") arrayfile::writearrayfile "

Write a list of arrays to file on disk  

Parameters
----------
* `filename` :  
    Filename to use  
* `arraylist` :  
    List of arrays to write  
* `mode` :  
    Mode for the file with designs  
* `nrows` :  
    If the list of arrays is empty, use this number of rows for the design file  
* `ncols` :  
    If the list of arrays is empty, use this number of rows for the design file  

Returns
-------
Value zero if succesfull  
";

%feature("docstring") arrayfile::writearrayfile "

Write a single array to file.  
";

%feature("docstring") arrayfile::append_arrayfile "

Append a single array to an array file. creates a new file if no file exists.  
";

%feature("docstring") arrayfile::selectArrays "

Make a selection of arrays from binary array file, append to list.  
";

%feature("docstring") arrayfile::selectArrays "

Select a single array from a file.  
";

%feature("docstring") arrayfile::selectArrays "

Make a selection of arrays.  
";

%feature("docstring") arrayfile::selectArrays "

Make a selection of arrays.  
";

%feature("docstring") arrayfile::selectArrays "

Make a selection of arrays, append to list.  
";

%feature("docstring") arrayfile::selectArrays "

Make a selection of arrays, append to list.  
";

%feature("docstring") arrayfile::keepElements "

From a container keep all elements with specified indices.  
";

%feature("docstring") arrayfile::removeElements "

From a container remove all elements with specified indices.  
";

%feature("docstring") arrayfile::selectArraysMask "

Make a selection of arrays from a list, append to list.  
";

%feature("docstring") arrayfile::appendArrays "

Append selection of arrays to existing list.  
";

%feature("docstring") arrayfile::appendArrays "

Append set of arrays to existing list.  
";

%feature("docstring") arrayfile::write_array_format "

Write a formatted array  
";

%feature("docstring") arrayfile::write_array_format "

Write an array to a file pointer.  
";

%feature("docstring") arrayfile::write_array_latex "

write an array in latex style  
";

%feature("docstring") arrayfile::convert_array_file "

Convert a file with arrays to a different format  
";

%feature("docstring") arrayfile::readbinheader "

Read header for binary data file. Return true if valid header file  

The header consists of 4 integers: 2 magic numbers, then the number of rows and
columns  
";

%feature("docstring") arrayfile::writebinheader "

Write header for binary data file.  
";

%feature("docstring") arrayfile::vector2doublebinfile "

Write a vector of numeric elements to binary file as double values.  
";

%feature("docstring") arrayfile::vectorvector2binfile "

Write a vector of vector elements to binary file.  
";

%feature("docstring") arrayfile::array2eigenX1 "

Convert 2-level array to main effects in Eigen format  

Parameters
----------
* `array` :  
    Array to convert  
* `intercept` :  
    If True, then include the intercept  

Returns
-------
The main effects model  
";

%feature("docstring") arrayfile::array2eigenX2 "

Convert 2-level array to second order interaction matrix in Eigen format  

The intercept and main effects are not included.  

Parameters
----------
* `array` :  
    Array to convert  

Returns
-------
The second order interaction model  
";

%feature("docstring") arrayfile::array2eigenModelMatrix "

Convert 2-level array to second order interaction model matrix (intercept, main
effects, interaction effects)  

Parameters
----------
* `array` :  
    Design of which to calculate the model matrix  

Returns
-------
Eigen matrix with the model matrix  
";

%feature("docstring") arrayfile::array2eigenModelMatrixMixed "

Create first and second order model matrix for mixed-level orthogonal array  

Parameters
----------
* `array` :  
    Input array  
* `verbose` :  
    Verbosity level  

Returns
-------
Pair with main effects and two-factor interaction model  

For 2-level arrays a direct calculation is used. For mixel-level arrays Helmert
contrasts are used.  
";

%feature("docstring") arrayfile::numberModelParams "

Calculate number of parameters in the model matrix  

A list of integers is returned, with the number of columns in:  

*   The intercept (always 1)  
*   The main effects  
*   The interaction effects (second order interaction terms without quadratics)  
*   The quadratic effects  

Parameters
----------
* `array` :  
    Orthogonal array or conference design  
* `order` :  
    Not used any more  

Returns
-------
List of sizes  
";

%feature("docstring") arrayfile::arrayInFile "

return index of specified array in a file. returns -1 if array is not found  

Parameters
----------
* `array` :  
    Array to find  
* `array_file` :  
    Location if file with arrays  
* `verbose` :  
    Verbosity level  

Returns
-------
Position of array in list  
";

%feature("docstring") arrayfile::arrayInList "

return index of specified array in a list. returns -1 if array is not found  

Parameters
----------
* `array` :  
    Array to find  
* `arrays` :  
    List of arrays  
* `verbose` :  
    Verbosity level  

Returns
-------
Position of array in list  
";

// File: conference_8h.xml

%feature("docstring") print_column "

print a candidate extension  
";

%feature("docstring") showCandidates "

Show a list of candidate extensions  

Parameters
----------
* `column_candidates` :  
    List of candidates to show  
";

%feature("docstring") conference2DSD "

Convert conference design to definitive screening design  

The DSD is created by appending the negated design to the conference design and
then appending a row of zeros.  

Parameters
----------
* `conference_design` :  
    Array with the conference design  
* `add_zeros` :  
    If True, then append a row of zeros  

Returns
-------
The DSD generated from the conference design  
";

%feature("docstring") reduceConference "

Reduce conference matrix to normal form using Nauty  

See also: reduceConferenceTransformation  
";

%feature("docstring") reduceConferenceTransformation "

Reduce conference matrix to normal form using Nauty  

The design is converted to a graph representation. The graph is then reduced
using Nauty to normal form and the resulting graph translated back to a
conference design.  

Parameters
----------
* `conference_design` :  
    Design to be reduced to normal form  
* `verbose` :  
    Verbosity level  

Returns
-------
A transformation that converts the input design to normal form  
";

%feature("docstring") extend_conference "

Extend a list of conference designs with a single column.  

The list of conference designs is extended by adding to each design the
candidate extentions generated by CandidateGenerator.  

Parameters
----------
* `lst` :  
    List of conference designs  
* `conference_type` :  
    Type specification for the conference designs  
* `verbose` :  
    Verbosity level  
* `select_isomorphism_classes` :  
    If True then select only a single design for each isomorphism class
    specified by the conference type.  

Returns
-------
List of generated conference designs  

The extension algorithm tried to generate designs in LMC0 normal form and prune
any designs that are not in LMC0 form.  
";

%feature("docstring") extend_conference_plain "

Extend a list of conference designs with a single column, plain version without
caching  

Research function.  
";

%feature("docstring") extend_conference_restricted "

Extend a list of conference designs with a single column  

Research function.  
";

%feature("docstring") extend_double_conference "

Extend a list of double conference matrices with an additional column  

The list of designs is extended by adding each design with the candidate
extentions generated by CandidateGeneratorDouble.  

Parameters
----------
* `lst` :  
    List of double conference designs  
* `conference_type` :  
    Type specification for the double conference designs  
* `verbose` :  
    Verbosity level  

Returns
-------
List of generated double conference designs  
";

%feature("docstring") selectConferenceIsomorpismClasses "

select representatives for the isomorphism classes of a list of conference
arrays  
";

%feature("docstring") selectConferenceIsomorpismIndices "

select representatives for the isomorphism classes of a list of conference
arrays, return indices of classes  
";

%feature("docstring") selectLMC0doubleconference "

Select double conference designs in LMC0 form  

Parameters
----------
* `list` :  
    List of double conference designs  
* `verbose` :  
    Verbosity level  
* `ctype` :  
    Specifiation of the class of designs  

Returns
-------
List with only the designs in the input list that are in LMC0 normal form.  
";

%feature("docstring") selectLMC0 "

Select conference designs in LMC0 form  

Parameters
----------
* `list` :  
    List of conference designs  
* `verbose` :  
    Verbosity level  
* `ctype` :  
    Specifiation of the class of designs  

Returns
-------
List with only the designs in the input list that are in LMC0 normal form.  
";

%feature("docstring") generateConferenceExtensions "

Generate candidate extensions for a conference design  

Parameters
----------
* `array` :  
    Design to be extended  
* `conference_type` :  
    Class of conference designs  
* `zero_index` :  
    index of zero in candidate column  
* `verbose` :  
    Verbosity level  
* `filtersymm` :  
    If True, filter based on symmetry  
* `filterj2` :  
    If True, filter based on J2 values  

Returns
-------
List of generated extensions  
";

%feature("docstring") generateConferenceRestrictedExtensions "

Generate candidate extensions for restricted isomorphism classes  
";

%feature("docstring") generateDoubleConferenceExtensions "

generate extensions for double conference matrices in LMC0 form  
";

%feature("docstring") generateSingleConferenceExtensions "

generate extensions for conference matrices in LMC0 form  
";

%feature("docstring") maxz "

return max position of zero in array, returns -1 if no zero is found  

The parameter k specifies the column to search in. For k=-1 all columns are
searched.  
";

%feature("docstring") compareLMC0 "

Return true if the first array is smaller in LMC-0 ordering than the second
array  
";

%feature("docstring") sortLMC0 "

sort list of conference designs according to LMC0 ordering  
";

%feature("docstring") LMC0checkDC "
";

%feature("docstring") LMC0check "
";

%feature("docstring") isConferenceFoldover "

return true if the design is a foldover array  
";

%feature("docstring") double_conference_foldover_permutation "

For a double conference design return a row permutation to a single conference
design  

If the design is not a foldover design then the first element of the returned
permutation is -1.  

Parameters
----------
* `double_conference` :  
    A double conference design  

Returns
-------
Permutation  
";

%feature("docstring") minz "

return minimal position of zero in specified column of a design  
";

// File: Deff_8h.xml

%feature("docstring") scoreD "

Calculate score from a set of efficiencies  

The score is the weighted sum of the efficiencies.  

Parameters
----------
* `efficiencies` :  
    Vector with calculated efficiencies  
* `weights` :  
    Weights for the efficiencies  

Returns
-------
Weighted sum of the efficiencies  
";

%feature("docstring") Doptimize "

Generates optimal designs for the specified class of designs  

The method uses a coordinate-exchange algorithm to optimze a target function
defined by the optimziation paramaters. The optimization is performed multiple
times to prevent finding a design in a local minmum of the target function.  

The method is described in more detail in \"Two-Level Designs to Estimate All
Main Effects and Two-Factor Interactions\", Eendebak et al., 2015,
Technometrics, https://doi.org/10.1080/00401706.2016.1142903.  

Parameters
----------
* `arrayclass` :  
    Class of designs to optimize  
* `nrestarts` :  
    Number of restarts to perform  
* `alpha` :  
    Optimization parameters. The target function is alpha_1 D + alpha_2 D_s +
    alpha D_1  
* `verbose` :  
    Verbosity level  
* `method` :  
    Method for optimization algorithm  
* `niter` :  
    Maximum number of iterations for each restart  
* `maxtime` :  
    Maximum calculation time. If this time is exceeded, the function is aborted  
* `nabort` :  
    Maximum number of iterations when no improvement is found  

Returns
-------
A structure with the generated optimal designs  
";

%feature("docstring") DoptimizeMixed "

Function to generate optimal designs with mixed optimization approach  

This function is beta code. See Doptimize for detauls of the parameters.  
";

%feature("docstring") optimDeff "

Optimize a design according to the optimization function specified.  

Arguments:  

Parameters
----------
* `array` :  
    Array to be optimized  
* `arrayclass` :  
    Structure describing the design class  
* `alpha` :  
    3x1 array with optimization parameters  
* `verbose` :  
    Verbosity level  
* `optimmethod` :  
    Optimization method to use  
* `niter` :  
    Number of iterations  
* `nabort` :  
    Number of iterations after which to abort when no improvements are found  

Returns
-------
Optimized designs  
";

// File: evenodd_8h.xml

%feature("docstring") processDepth "

Extend arrays using a depth-first or breadth-first approach  

Parameters
----------
* `goodarrays` :  
    List of arrays to extend  
* `depthalg` :  
    Extend using depth-first or breadth-first  
* `dextend` :  
    Option structure for the extension  
* `dextendsublight` :  
    Data structure for the extensions  
* `extensioncol` :  
    Column to extend  
* `verbose` :  
    Verbosity level  
";

%feature("docstring") depth_extend_hybrid "

depth-first extension of arrays. depending on the symmetry group of the array to
be extended a direct method is used or a method with caching of candidate
columns  
";

%feature("docstring") depth_extend_direct "

variation of depth_extend for arrays with large symmetry groups  
";

%feature("docstring") depth_extend_array "

depth extend a single array  
";

%feature("docstring") calculateArrayParetoJ5Cache "
";

%feature("docstring") addArraysToPareto "

add arrays to set of Pareto results  
";

%feature("docstring") addArraysToPareto "

add arrays to set of Pareto results  
";

%feature("docstring") readStatisticsFile "

read statistics object from disk  
";

%feature("docstring") writeStatisticsFile "

write statistics object to disk  
";

%feature("docstring") calculateJstatistics "

calculate J-value statistics  
";

// File: extend_8h.xml

%feature("docstring") extend_arraylist "

Extend a list of orthogonal arrays  

Parameters
----------
* `array_list` :  
    The list of arrays to be extended  
* `array_class` :  
    Class of arrays to generate  
* `oaextend_options` :  
    Parameters for the extension algorithm  

Returns
-------
List of all generated arrays  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  
";

%feature("docstring") extend_arraylist "

Extend a list of arrays with default options  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  
";

%feature("docstring") extend_arraylist "

Extend a list of orthogonal arrays  

Parameters
----------
* `array_list` :  
    The list of arrays to be extended  
* `array_class` :  
    Class of arrays to generate  
* `oaextend_options` :  
    Parameters for the extension algorithm  

Returns
-------
List of all generated arrays  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  

Parameters
----------
* `extensioncol` :  
    Index of column to be added to the designs  
* `extensions` :  
    List to append generated designs to  

Returns
-------
Number of candidate arrays generated  
";

%feature("docstring") extend_array "

Extend a single orthogonal array  

Parameters
----------
* `array` :  
    The array to be extended  
* `array_class` :  
    Class of arrays to generate  
* `oaextend` :  
    Parameters for the extension algorithm  
";

%feature("docstring") extend_array "

Extend a single orthogonal array with the default LMC algorithm  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  
";

%feature("docstring") extend_array "

Extend an orthogonal array with a single column  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  

Parameters
----------
* `array` :  
    Array to extend  
* `arrayclass` :  
    Array data for the full array  
* `extension_column` :  
    Column to extend  
* `extensions` :  
    List to which generated valid extensions are added  
* `oaextend` :  
    Structure with options  

Returns
-------
Number of candidate extensions generated  
";

%feature("docstring") runExtendRoot "

Run the LMC extension algorithm starting with the root array  

See also: extend_array(const array_link &, arraydata_t &, OAextend const &)  
";

// File: graphtools_8h.xml

%feature("docstring") nauty::transformGraph "

Apply a vertex permutation to a graph.  
";

%feature("docstring") nauty::reduceOAnauty "

Reduce an orthogonal array to Nauty minimal form. the array transformation is
returned.  
";

%feature("docstring") nauty::reduceOAnauty "

Reduce an orthogonal array to Nauty minimal form. the array transformation is
returned.  
";

%feature("docstring") nauty::array2graph "

Convert orthogonal array to graph representation  

The conversion method is as in Ryan and Bulutoglu. The resulting graph is bi-
partite. The graph representation can be used for isomorphism testing.  
";

%feature("docstring") nauty::array2graph "

Convert orthogonal array to graph representation  

The conversion method is as in Ryan and Bulutoglu. The resulting graph is bi-
partite. The graph representation can be used for isomorphism testing.  
";

%feature("docstring") nauty::oagraph2transformation "

From a relabelling of the graph return the corresponding array transformation.  
";

// File: InfInt_8h.xml

%feature("docstring") my_div "
";

%feature("docstring") my_ldiv "
";

%feature("docstring") my_lldiv "
";

// File: lmc_8h.xml

%feature("docstring") algorithm_t_list "
";

%feature("docstring") algnames "

return name of the algorithm  
";

%feature("docstring") apply_hadamard "

Apply Hadamard transformation to orthogonal array.  
";

%feature("docstring") acquire_LMCreduction_object "

return static structure from dynamic global pool, return with
releaseGlobalStatic  
";

%feature("docstring") release_LMCreduction_object "
";

%feature("docstring") clear_LMCreduction_pool "

release all objects in the pool  
";

%feature("docstring") insert_if_not_at_end_of_vector "

Append element to vector if the element the element is not at the end of vector.  
";

%feature("docstring") is_root_form "

Return True if the array is in root form  

Parameters
----------
* `array` :  
    Array to check  
* `strength` :  
    Strength to use  

Returns
-------
True if the array is in root form for the specified strength  
";

%feature("docstring") LMCreduction_train "

helper function for LMC reduction  
";

%feature("docstring") LMCcheck "

Perform LMC check or reduction on an array.  
";

%feature("docstring") LMCcheck "

Perform LMC check or reduction on an array.  
";

%feature("docstring") LMCcheck "

Perform LMC check on an orthogonal array  

Parameters
----------
* `array` :  
    Array to be checked for LMC minimal form  

Returns
-------
Result of the LMC check  
";

%feature("docstring") LMCcheckOriginal "

Perform LMC check on a 2-level orthogonal array  

The algorithm used is the original algorithm from \"Complete enumeration of
pure-level and mixed-level orthogonal arrays\", Schoen et al, 2009  

Parameters
----------
* `array` :  
    Array to be checked for LMC minimal form  

Returns
-------
Result of the LMC check  
";

%feature("docstring") reduceArraysGWLP "

reduce arrays to canonical form using delete-1-factor ordering  
";

%feature("docstring") reductionDOP "

Caculate the transformation reducing an array to delete-on-factor normal  

The normal form is described in \"A canonical form for non-regular arrays based
on generalized wordlength pattern values of delete-one-factor projections\",
Eendebak, 2014  

Parameters
----------
* `array` :  
    Orthogonal array  
* `verbose` :  
    Verbosity level  

Returns
-------
The transformation that reduces the array to normal form  
";

%feature("docstring") reduceDOPform "

Reduce an array to canonical form using delete-1-factor ordering  

The normal form is described in \"A canonical form for non-regular arrays based
on generalized wordlength pattern values of delete-one-factor projections\",
Eendebak, 2014  

Parameters
----------
* `array` :  
    Orthogonal array  
* `verbose` :  
    Verbosity level  

Returns
-------
The array transformed to normal form  
";

%feature("docstring") selectUniqueArrays "

select the unique arrays in a list, the original list is sorted in place. the
unique arrays are append to the output list  
";

%feature("docstring") projectionDOFvalues "

Calculate projection values for delete-of-factor algorithm  
";

%feature("docstring") reduceLMCform "

reduce an array to canonical form using LMC ordering  
";

%feature("docstring") LMCcheckLex "

Apply LMC check (original mode) to a list of arrays  
";

%feature("docstring") LMCcheckLex "

Perform minimal form check with LMC orderin.  
";

%feature("docstring") LMCcheckj4 "

Perform minimal form check with J4 ordering.  
";

%feature("docstring") LMCcheckj5 "

Perform minimal form check for J5 ordering.  
";

%feature("docstring") print_rowsort "

Print the contents of a rowsort structure.  

Parameters
----------
* `rowsort` :  
    Pointer to rowsort structure  
* `N` :  
    Number of elements  
";

%feature("docstring") print_column_rowsort "
";

// File: mathtools_8h.xml

%feature("docstring") vectormax "

Return maximum element of a std::vector.  
";

%feature("docstring") vectormin "

Return minimum element of a std::vector.  
";

%feature("docstring") cumsum "

calculate cumulative sum of a vector  
";

%feature("docstring") cumsum0 "

calculate cumulative sum of a vector with added zero  
";

%feature("docstring") cumsum0 "

calculate cumulative sum of a vector with added zero  
";

%feature("docstring") permutation "

create permutation of specified length  
";

%feature("docstring") array2vector "

convert array given by pointer to std::vector  
";

%feature("docstring") array2larray "

convert array given by pointer to larray  
";

%feature("docstring") print_perm "

Print permutation.  

Prints a permutation to output stream  

Parameters
----------
* `out` :  
    Output stream  
* `s` :  
    Pointer to start of array  
* `len` :  
    Length of array to be printed  
* `maxlen` :  
    (optional) Maximum length to print  
";

%feature("docstring") print_perm "
";

%feature("docstring") print_perm "
";

%feature("docstring") print_perm "

print permutation with string in front  
";

%feature("docstring") print_perm "
";

%feature("docstring") print_perm "
";

%feature("docstring") print_perm "
";

%feature("docstring") print_perm_int "
";

%feature("docstring") compare_matrix "

Compare two arrays and return whether equal or not.  

Parameters
----------
* `A` :  
    Pointer to array  
* `B` :  
    Pointer to array  
* `number_of_rows` :  
    Number of rows  
* `number_of_columns` :  
    Number of columns  

Returns
-------  
";

%feature("docstring") factorial "

Calculates factorial A small function that calculates the factorial of a number.
Returns one if the argument is smaller or equal to 1.  

Parameters
----------
* `number` :  
    Number to calculate the factorial of  

Returns
-------
Factorial of specified number  
";

%feature("docstring") factorial_return_argument "

A small function that calculates the factorial of a number. This inline
factorial function is the same as the standard factorial calculation, except
that the return type is generic Returns one if the argument is smaller or equal
to 1  

Parameters
----------
* `number` :  
    number to calculate factorial of  
";

%feature("docstring") ncombs "

Calculates number of combinations.  

The number of combinations is calculated using the an addapted formula  

Parameters
----------
* `n` :  
    Total number of entries to choose from  
* `k` :  
    Number of entries in a certain combination  
";

%feature("docstring") ncombsm "

calculate using multiplicative formula, see
http://en.wikipedia.org/wiki/Binomial_coefficient#Computing_the_value_of_binomial_coefficients  
";

%feature("docstring") next_perm "
";

%feature("docstring") next_perm "
";

%feature("docstring") fastrand "
";

%feature("docstring") seedfastrand "
";

%feature("docstring") fastrandK "
";

%feature("docstring") set_srand "

seed the C rand method with the srand function  
";

%feature("docstring") next_perm_twoperm "
";

%feature("docstring") random_perm "

Create random permutation using Fisher-Yates shuffle, or Knuth shuffle  

The permutation is peformed inplace.  

Parameters
----------
* `array` :  
    Array of objects  
* `length` :  
    Length of array  
";

%feature("docstring") random_perm "

Create random permutation using Fisher-Yates shuffle, or Knuth shuffle  
";

%feature("docstring") new_comb_init "

Create a new combination and initialize.  
";

%feature("docstring") delete_comb "

Delete combination.  
";

%feature("docstring") init_comb "

Initialize a combination  

Parameters
----------
* `comb` :  
    Pointer to combination array  
* `k` :  
    Number of elements  
* `n` :  
    Numbers to choose from  

Returns
-------
Number of combinations possible  
";

%feature("docstring") next_combination "

Gives combination nr k.  

Gives combination number k based on an algorithm from wikipedia.  

Parameters
----------
* `comb` :  
    Pointer to combination  
* `k` :  
    Number of the current combination  
* `n` :  
    Number of elements in combination  
";

%feature("docstring") next_combination_fold "
";

%feature("docstring") perm_is_ordered "

Check whether a permutation is ordered or not.  

Parameters
----------
* `perm` :  
* `len` :  

Returns
-------
Return 1 of an ordered permutation, 0 otherwise  
";

%feature("docstring") new_perm "

Create a new permutation.  
";

%feature("docstring") clone_perm "

Create a new permutation.  
";

%feature("docstring") delete_perm "

Delete a permutation.  
";

%feature("docstring") invert_permutation "

Invert a permutation.  

Parameters
----------
* `perm` :  
    Permutation as integer type std::vector  

Returns
-------
New permutation that is the inverse of the argument  
";

%feature("docstring") invert_permutation "

Invert a permutation.  

Parameters
----------
* `perm` :  
    Permutation as integer type std::vector  
* `iperm` :  
    Output permutation  
";

%feature("docstring") invert_permutation "

Invert a permutation.  

Parameters
----------
* `perm` :  
    Pointer to permutation  
* `len` :  
* `iperm` :  
    Pointer to new permutation that is the inverse of the argument  
";

%feature("docstring") invert_permutation "

Invert a permutation.  

Parameters
----------
* `perm` :  
    Pointer to permutation  
* `len` :  

Returns
-------
Pointer to new permutation that is the inverse of the argument  
";

%feature("docstring") perform_level_perm "

Perform level permutation on an array  

Parameters
----------
* `src` :  
    Pointer to array  
* `n` :  
    Length of array  
* `perm` :  
    Permutation to perform  
";

%feature("docstring") perform_level_perm "

Perform a permutation on a set of data elements.  

Parameters
----------
* `src` :  
* `target` :  
* `n` :  
* `perm` :  
";

%feature("docstring") composition_perm "

Calculate composition of 2 permutations.  

Calculates C = B &#9675; A  

Parameters
----------
* `A` :  
* `B` :  
* `n` :  
    Length of permutations  
* `C` :  
";

%feature("docstring") composition_perm "

Calculate composition of 2 permutations.  

Calculates C = B &#9675; A  

Parameters
----------
* `A` :  
* `B` :  
* `C` :  
";

%feature("docstring") perform_perm "

Perform a permutation on a set of objects.  

Parameters
----------
* `src` :  
* `target` :  
* `n` :  
* `perm` :  
";

%feature("docstring") perform_perm "

Perform a permutation on a set of objects.  
";

%feature("docstring") perform_inv_perm "

Perform inverse permutation.  
";

%feature("docstring") perform_inv_perm "

Perform inverse permutation.  
";

%feature("docstring") perform_inv_perm "

Perform inverse permutation.  
";

%feature("docstring") init_perm "

Initialize a permutation  

Parameters
----------
* `perm` :  
";

%feature("docstring") init_perm "

Initialiaze a permutation  

Parameters
----------
* `perm` :  
* `len` :  
";

%feature("docstring") init_signperm "

Initialiaze a sign permutation with all +1s  

Parameters
----------
* `signperm` :  
    Permutation  
";

%feature("docstring") compare_perm "

return true if two permutations are equal  
";

%feature("docstring") copy_perm "

copy a permuntation  
";

%feature("docstring") init_perm_n "

initialize a permutation and return the number of permutations  
";

%feature("docstring") new_perm_init "

create a new permutation and initialize  
";

%feature("docstring") issorted "
";

%feature("docstring") new_valueindex "
";

%feature("docstring") init_valueindex "
";

%feature("docstring") argsort "
";

%feature("docstring") permute "

Permute a std::vector.  
";

%feature("docstring") permuteback "

Permute a std::vector with inverse permutation.  
";

%feature("docstring") symm_group_index_plain "

Calculate symmetry groups of a list of integers under permutations.  

Parameters
----------
* `vec` :  
* `n` :  
* `idx` :  
* `gstart` :  
* `gsize` :  

Returns
-------
Number of groups found  
";

%feature("docstring") ipow "

Power of two integers.  
";

%feature("docstring") ipow "

Power of two unsigned integers.  
";

%feature("docstring") power_minus_one "

-1 to the power n (integer)  
";

%feature("docstring") power_minus_one "

-1 to the power n (integer)  
";

%feature("docstring") krawtchouk "

calculate value of Krawtchouk polynomial  
";

%feature("docstring") krawtchouksCache "

calculate value of Krawtchouk polynomial  
";

%feature("docstring") krawtchouks "

calculate value of Krawtchouk polynomial  

See also: https://en.wikipedia.org/wiki/Kravchuk_polynomials  
";

%feature("docstring") conditionNumber "

return the condition number of a matrix  
";

%feature("docstring") average "

calculate the average value of a vector of numbers  
";

%feature("docstring") fraction_nonzero "

calculate the fraction of non-zero elemenets of a vector  
";

// File: md5_8h.xml

%feature("docstring") md5 "

calculate md5 sum of a data block in memory  
";

%feature("docstring") md5 "

calculate md5 sum of a file on disk  
";

// File: msstdint_8h.xml

// File: nonroot_8h.xml

%feature("docstring") LMCreduce_non_root "

default reduction function for non-root stage  
";

%feature("docstring") LMCreduce_non_root_j4 "

default reduction function for non-root stage with J4 algorithm  
";

%feature("docstring") LMCreduce_non_root_2level "

specialized reduction function for 2-level arrays  
";

// File: oaoptions_8h.xml

%feature("docstring") compile_information "

Print the compile-time options to string.  

Returns
-------
String with compile time information  
";

%feature("docstring") version "

Print version.  
";

%feature("docstring") print_copyright "

Print copyright statement.  
";

%feature("docstring") print_copyright_light "

Print copyright statement.  
";

%feature("docstring") print_options "

Print compile time options to output stream.  
";

%feature("docstring") print_options "

Print compile time options to stdout.  
";

%feature("docstring") oadevelop "
";

// File: pareto_8h.xml

// File: printfheader_8h.xml

// File: strength_8h.xml

%feature("docstring") create_reverse_colcombs_fixed "
";

%feature("docstring") check_divisibility "

Checks on the divisibility of the number of runs by the product of the levels in
the factors for all t-tuples.  
";

%feature("docstring") check_divisibility "

check whether an array passes divisibility test  
";

%feature("docstring") add_element_freqtable "

Add row to frequency table using cache system.  
";

%feature("docstring") add_element_freqtable_col "

fast version of add_element_freqtable  
";

%feature("docstring") recount_frequencies "
";

%feature("docstring") strength_check "

Perform strength check on an array  

Special case for extension of an array with proper strength  
";

%feature("docstring") strength_check "

perform strength check on an array  
";

%feature("docstring") valid_element "

Determine whether an element passes the strength test.  

Parameters
----------
* `es` :  
* `position` :  
* `array` :  
    Pointer to array  

Returns
-------  
";

%feature("docstring") valid_element_2level "

Determine whether an element passes the strength test, specialized for 2-level
array  
";

// File: timsort_8hpp.xml

// File: tools_8h.xml

%feature("docstring") printfd_handler "

function to print debugging messages  
";

%feature("docstring") log_print "
";

%feature("docstring") getloglevel "

return current level of logging  
";

%feature("docstring") setloglevel "

reset the level of logging  
";

%feature("docstring") checkloglevel "

return True if the current logging level is smaller or equal than the specified
level  
";

%feature("docstring") logstream "

log a stream to stdout if level specied is higher than the current logging level  
";

%feature("docstring") system_uname "

Return string describing the system.  
";

%feature("docstring") path_separator "

return path separator symbol for the current platform  
";

%feature("docstring") mycheck_handler "

handler for error messages. throws an std::runtime_error exception  
";

%feature("docstring") myassert "

Check whether the condition is true and throw an expception otherwise.  
";

%feature("docstring") throw_runtime_exception "

Throw a runtime_error exception with specified error message  

This exception is caught in the SWIG interface.  
";

%feature("docstring") cprintf "

conditional printf  
";

%feature("docstring") flush_stdout "

flush to stdout  
";

%feature("docstring") safedelete "

Delete object given by a pointer and set to zero.  
";

%feature("docstring") safedeletearray "

Delete array and set pointer to zero  

Parameters
----------
* `pointer` :  
    Pointer to allocated array  
";

%feature("docstring") next_comb "

Gives next combination.  

Gives next combination for k elements out of n based on an algorithm from
wikipedia. The generate is sorted.  

Parameters
----------
* `comb` :  
    Pointer to combination  
* `k` :  
    Number of the current combination  
* `n` :  
    Number of elements in combination  
";

%feature("docstring") next_comb "

Go to next combination in sequence  
";

%feature("docstring") next_comb_s "

Go to next combination in sequence  
";

%feature("docstring") swap_object "

Template to swap two objects of arbitrary datatype Please use std::swap instead.  

Parameters
----------
* `a` :  
* `b` :  
";

%feature("docstring") malloc2d_nelements "

Calculate the number of elements in a 2D table with rows with different sizes  
";

%feature("docstring") malloc2d_irr "

Allocate a 2-dimensional array with non-uniform rows.  

Parameters
----------
* `nrows` :  
    Number of rows in the table  
* `rowsizes` :  
    Size of each row  

Returns
-------  
";

%feature("docstring") malloc2d_irr "

Allocate a 2-dimensional array with non-uniform rows, return size of allocated
space.  

Parameters
----------
* `nrows` :  
    Number of rows in the table  
* `rowsizes` :  
    Size of each row  
* `nelements` :  
    This parameter is initialized with the size of the array allocated  

Returns
-------  
";

%feature("docstring") malloc2d "

Allocate a 2-dimensional array of specified size.  

Parameters
----------
* `nrows` :  
    Number of rows  
* `rowsize` :  
    Size of each row  

Returns
-------  
";

%feature("docstring") free2d "

Release a 2-dimensional array.  

Parameters
----------
* `data` :  
* `nrows` :  
";

%feature("docstring") free2d "

Release a 2-dimensional array.  

Parameters
----------
* `data` :  
";

%feature("docstring") free2d_irr "

Release a 2-dimensional non-uniform array.  

Parameters
----------
* `data` :  
    Pointer to allocated array  
* `nrows` :  
    Not used at the moment  
";

%feature("docstring") print_array "
";

%feature("docstring") print_array "
";

%feature("docstring") display_vector "

print vector using generic std::cout print functionality  
";

%feature("docstring") printf_vector "

print vector using printf function  

Parameters
----------
* `vector` :  
    Vector to be displayed  
* `format` :  
    Format to use in printf  
* `separator` :  
    Separator symbol to use  
";

%feature("docstring") get_time_ms "

return time with milisecond precision  
";

%feature("docstring") get_time_ms "

return time difference with milisecond precision  
";

%feature("docstring") trim "

trim a string by removing the specified characters from the left and right  
";

%feature("docstring") currenttime "

return the current time as a string  
";

%feature("docstring") oafilestring "

return string describing array  
";

%feature("docstring") itos "

Convert integer to C++ string.  

Parameters
----------
* `integer_value` :  
    Integer  

Returns
-------
String representation of the integer  
";

%feature("docstring") printfstring "

printf-style function that returns std::string  
";

%feature("docstring") insertionSort "

Template for insertionSort.  

Parameters
----------
* `array` :  
    Data to be sorted  
* `length` :  
    Length of array  
";

%feature("docstring") bubbleSort "

sort arrays using bubbleSort  
";

%feature("docstring") flipSort "

Sorting similar to bubblesort but fast for sorted arrays  

The indices left and right are inclusive.  
";

%feature("docstring") bubbleSort2 "

Template for bubble sort.  

Parameters
----------
* `array` :  
    Array to be sorted  
* `array_size` :  
    Size of the array  
";

%feature("docstring") quickSort "

sort list using quickSort  
";

%feature("docstring") shellSort "
";

%feature("docstring") replaceString "

replace all occurces of a substring in a string  
";

%feature("docstring") printdoubleasbits "

print a double value as bits  
";

%feature("docstring") splitDir "

calculate directory name for job splitted into parts  
";

%feature("docstring") splitFile "

calculate file name of job splitted into parts  
";

%feature("docstring") splitTag "

calculate tag for job splitted into parts  
";

// File: unittests_8h.xml

%feature("docstring") unittest_reduceConferenceTransformation "
";

%feature("docstring") unittest_nautynormalform "
";

%feature("docstring") checkTransformationComposition "
";

%feature("docstring") test_array_manipulation "
";

%feature("docstring") checkConferenceComposition "
";

%feature("docstring") test_conference_candidate_generators "
";

%feature("docstring") checkTransformationInverse "
";

%feature("docstring") checkConferenceInverse "
";

%feature("docstring") testLMC0checkDC "
";

// File: version_8h.xml

// File: dir_68267d1309a1af8e8297ef4c3efbcdba.xml

// File: indexpage.xml

