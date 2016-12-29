
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(OApackage)

# ----------------------------------------------------------------------------
#   External packages
# ----------------------------------------------------------------------------

if(WIN32)
	#message(STATUS "Set build type...")
	set(CMAKE_CONFIGURATION_TYPES Release CACHE TYPE INTERNAL FORCE) # no Debug
	#set(CMAKE_BUILD_TYPE Release CACHE STRING)
endif()

if (NOT WIN32)
  find_package(Boost QUIET COMPONENTS system filesystem  )
  message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
  if (Boost_FOUND)
    message(STATUS "Boost was found")
  endif()
endif()

find_package(ZLIB)	# Needed for Python modules
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
endif()
message(STATUS "ZLIB_FOUND: ${ZLIB_FOUND}")
if(ZLIB_FOUND)
	#message(STATUS "ZLIB_INCLUDE_DIRS: ${ZLIB_INCLUDE_DIRS}")
	#message(STATUS "ZLIB_LIBRARIES: ${ZLIB_LIBRARIES}")
	include_directories(${ZLIB_INCLUDE_DIRS})
endif()

find_package(SWIG)	# Needed for Python modules
find_package(NumPy QUIET)	# Needed for Python modules
if (NumPy_FOUND)
  message(STATUS "NumPy was found")
endif()

#find_package(MPI)	# Needed for legacy code and OpenMP functionality

# ----------------------------------------------------------------------------
#   Code
# ----------------------------------------------------------------------------

include_directories(.)
include_directories(src)
include_directories(nauty)
link_directories(${CMAKE_SOURCE_DIR})

if (NOT EXISTS "${CMAKE_SOURCE_DIR}/dev/oadevelop.h") 
  set(OADEV 0 CACHE BOOL "Enable development code")
else()
  set(OADEV 1 CACHE BOOL "Enable development code")
endif()

if (NOT EXISTS "${CMAKE_SOURCE_DIR}/src/extend.h") 
  set(OAEXTEND 0 CACHE BOOL "Enable OA generation code")
else()
  set(OAEXTEND 1 CACHE BOOL "Enable OA generation code")
endif()


message(STATUS "OA development code: ${OADEV}")
message(STATUS "OA extend code: ${OAEXTEND}")

set(headers src/pareto.h src/analysis.h src/anyoption.h src/arrayproperties.h src/arraytools.h  src/graphtools.h src/mathtools.h src/oaoptions.h src/tools.h src/Deff.h src/evenodd.h)
set(srcs src/pareto.cpp src/anyoption.cpp src/md5.cpp src/strength.cpp src/arrayproperties.cpp src/graphtools.cpp src/arraytools.cpp src/mathtools.cpp src/tools.cpp src/oaoptions.cpp src/bitarray/bit_array.cpp src/Deff.cpp src/evenodd.cpp)

set(headersextend ${headers} src/extend.h src/strength.h src/lmc.h src/nonroot.h )
set(srcsextend ${srcs} src/extend.cpp src/lmc.cpp src/nonroot.cpp )

set(headersextend ${headersextend} src/InfInt.h)
set(srcsextend ${srcsextend} src/InfInt.cpp)

set(headersextend ${headersextend} src/conference.h)
set(srcsextend ${srcsextend} src/conference.cpp)

if(OADEV)
  set(headersextend ${headersextend} dev/oadevelop.h)
  set(srcsextend ${srcsextend} dev/oadevelop.cpp)
  include_directories(dev)
endif()

set(OUTPUTDIR ${CMAKE_SOURCE_DIR})

if(OADEV)
  # if development code, try to use openmp
  find_package(OpenMP )

if (0)
set(OPENMP_FOUND False)
message(STATUS "Setting OpenMP to False!")
endif()

if(OPENMP_FOUND AND OADEV)
   message(STATUS "OpenMP found: OpenMP_C_FLAGS ${OpenMP_C_FLAGS}, OpenMP_EXE_LINKER_FLAGS ${OpenMP_EXE_LINKER_FLAGS}")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DDOOPENMP")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DDOOPENMP")
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
   message(STATUS "OpenMP not found! building code without OpenMP support")
endif()
endif(OADEV)

# for ccache
if(NOT WIN32 AND 0)
  if("$ENV{USER}" STREQUAL "eendebakpt")
    message(STATUS "Using ccache for compilation")
    set(CMAKE_CXX_COMPILER "ccache")
    set(CMAKE_CXX_COMPILER_ARG1 "g++")
  endif()
endif()

# used for HPC at Antwerpen
#option(INTEL_CC "Use the Intel compiler." OFF) 

site_name(HOSTNAME)
string(FIND ${HOSTNAME} hopper VSC)
if(VSC GREATER 0) 
  set(VSC TRUE)
else()
  set(VSC FALSE)
endif()
message(STATUS "Hostname ${HOSTNAME}, VSC ${VSC}")

if (VSC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") # -DMPICH_IGNORE_CXX_SEEK -axSSE4.2,SSE4.1
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -axSSE4.2,SSE4.1")
endif()

if (INTEL_CC) 
  find_program(CMAKE_C_COMPILER NAMES icc) 
  find_program(CMAKE_CXX_COMPILER NAMES icpc) 
  find_program(CMAKE_AR NAMES xiar) 
  find_program(CMAKE_LINKER NAMES xild) 

  if (CMAKE_C_COMPILER MATCHES CMAKE_C_COMPILER-NOTFOUND OR 
      CMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND OR 
      CMAKE_AR MATCHES CMAKE_AR-NOTFOUND OR 
      CMAKE_LINKER MATCHES CMAKE_LINKER-NOTFOUND) 
    message(FATAL_ERROR "Cannot find Intel compiler.  You may need to run `.  /opt/intel/bin/compilervars.sh intel64'") 
  endif () 

  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPICH_IGNORE_CXX_SEEK -axSSE4.2,SSE4.1 -O3")
  message(STATUS "At HPC cluster")
endif (INTEL_CC) 

set(LINKOPTS "")
#message(STATUS "ENV{USERNAME} $ENV{USERNAME}, USER $ENV{USER}")
if( "$ENV{USER}" STREQUAL "vsc20149")
  set(LINKOPTS "-s")	# strip executables for file size
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  #message(STATUS "CMAKE_COMPILER_IS_GNUCXX ${CMAKE_COMPILER_IS_GNUCXX}")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wreturn-type -Wno-unused-but-set-variable -Wno-unused-result -fPIC -g -fpermissive")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wreturn-type  -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unused-result -Wno-unknown-pragmas -fPIC  -g") # -fpermissive  # 
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
  set(COPTS "${COPTS}  -O3 ") # -fpermissive
  set(COPTS "${COPTS}  -DEIGEN_MPL2_ONLY")
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings") # remove warnings

  set(COPTS "${COPTS} -Wlong-long -pedantic -Wno-variadic-macros")	# make code comply to R
endif ()

# for debugging
if("$ENV{USER}" STREQUAL "eendebakpt" )
  #set(COPTS "-g -O2 -march=native" ) # -std=c++0x ") # -O3 does not always give the fastest results for executables!
  #set(COPTS "-g -O0" ) # -std=c++0x ")
endif()

if("$ENV{OADEBUG}" STREQUAL "1")
  #set(COPTS "-g -O2 -march=native" ) # -std=c++0x ") # -O3 does not always give the fastest results for executables!
  set(COPTS "-g -O0 -DOADEBUG" ) # -std=c++0x -Weffc++") 
endif()

if(OADEV)
  if(WIN32)
  else()
    set(LINKOPTS "${LINKOPTS}")
  endif()
else()
  if(WIN32)
  else()
    set(LINKOPTS "${LINKOPTS} -s")
  endif()
endif()

if(WIN32)
  set(COPTS "/MT /wd4018 /wd4996")  
endif()
if (MSVC)
  message(STATUS "MSCV compiler")
  set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /EHsc /W0 /wd4996 /D_SCL_SECURE_NO_WARNINGS")
  ADD_DEFINITIONS(-DNODEFAULTLIB:library)
  set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  #set ( CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} /NODEFAULTLIB:MSVCRT /EHsc /W0 /MT /verbose" )
 # SET( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_STATIC_LINKER_FLAGS} "/NODEFAULTLIBS:LIBCMT" )
endif()

set(COPTS "${COPTS} -DFULLPACKAGE")

if(ZLIB_FOUND AND 1)
  set(COPTS "${COPTS} -DUSEZLIB")
  include_directories(ZLIB_INCLUDE_DIRS)
else()
  set(COPTS "${COPTS} -DNOZLIB")
endif()

if (VERBOSE)
	message(STATUS "copts ${COPTS}")
endif()


if(Boost_FOUND)
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
  set(COPTS "${COPTS} -DHAVE_BOOST")
endif()  

if(OADEV)
  set(COPTS "${COPTS} -DOADEV")
endif()

set(POG 0)	# Profile Guided Optimization

#set(COPTS "${COPTS} -ffast-math")

if(POG STREQUAL "1")
  message(STATUS "POG: profile")
  message(STATUS "POG: run from command line: ")
  set(COPTS "${COPTS} -fprofile-generate --coverage")
  set(LINKOPTS "${LINKOPTS} -fprofile-generate --coverage")
  
 # set(COPTS "${COPTS} -prof-gen -fprofile-generate --coverage ")
 # set(LINKOPTS "${LINKOPTS} -prof-gen -fprofile-generate --coverage")  
endif()
if(POG STREQUAL "2")
  message(STATUS "POG: use")
  set(COPTS "${COPTS} -fprofile-use")
  set(LINKOPTS "${LINKOPTS} -fprofile-use")
endif()


macro(add_compile_flags _target)
  # Provide user feedback
#  message(STATUS "  Adding compile flags to: ${_target}")
  set_target_properties(${_target} PROPERTIES COMPILE_DEFINITIONS "OAEXTEND_SINGLECORE;NEWINTERFACE " )
  set_target_properties(${_target} PROPERTIES COMPILE_FLAGS "${COPTS}" )
  set_target_properties(${_target} PROPERTIES LINK_FLAGS "${LINKOPTS}")

endmacro()


############### Programs ################################
  
include_directories(.)

set(nautysrc src/nauty/nauty.c src/nauty/nautinv.c src/nauty/nautil.c src/nauty/naurng.c src/nauty/naugraph.c src/nauty/schreier.c src/nauty/naugroup.c) #  nauty/gtools.c
include_directories(src/nauty)

set(graphsrc )

set(progs oacat oajoin oapareto oasplit oaanalyse oainfo oaconvert oafilter oareduceDOP pareto_example oaconference oaconfcheck)
set(progsextend )

if(OAEXTEND)
  list(APPEND progs oacheck oastreaming oaunittest)
endif()
if(OADEV)
  list(APPEND progs oatest oaranktest oa_depth_extend oaextendA oaclustergather)
endif()
   
if (WIN32)
  add_library(oalib STATIC ${srcsextend} ${headersextend} ${nautysrc} ${graphsrc})
  add_compile_flags(oalib)
  add_library(oalibshared SHARED ${srcsextend} ${headersextend} ${nautysrc} ${graphsrc})
  add_compile_flags(oalibshared)
  if(ZLIB_FOUND)
	TARGET_LINK_LIBRARIES(oalib ${ZLIB_LIBRARIES})
	TARGET_LINK_LIBRARIES(oalibshared ${ZLIB_LIBRARIES})
  endif()

else()
  add_library(oalib STATIC ${srcsextend} ${headersextend} ${nautysrc} ${graphsrc})
  add_compile_flags(oalib)
  add_library(oalibshared SHARED ${srcsextend} ${headersextend} ${nautysrc} ${graphsrc})
  add_compile_flags(oalib)
    TARGET_LINK_LIBRARIES(oalib m)
    TARGET_LINK_LIBRARIES(oalibshared m)
  if(ZLIB_FOUND)
	TARGET_LINK_LIBRARIES(oalib ${ZLIB_LIBRARIES})
	TARGET_LINK_LIBRARIES(oalibshared ${ZLIB_LIBRARIES})
  endif()

endif()

if (OAEXTEND)
  add_executable(oaextendsingle "utils/oaextend.cpp" ${headersextend})
  add_compile_flags(oaextendsingle)

  # get_directory_property( DirDefs COMPILE_DEFINITIONS ) 
  #message( "COMPILE_DEFINITIONS = ${DirDefs}" )
  
  TARGET_LINK_LIBRARIES(oaextendsingle oalib)
  set(extendprogs oaextendsingle )


  if (MPI_FOUND AND NOT WIN32 AND 1)
    message(STATUS "Found MPI package: adding oaextendmpi")
    if (VERBOSE)
	message(STATUS "MPI executable: ${MPIEXEC}, libs ${MPI_LIBRARIES}")
    endif()
    #message(STATUS "   libs ${MPI_LIBRARIES}; flags ${MPI_COMPILE_FLAGS}; libs ${MPI_LIBRARIES}")
    include_directories( ${MPI_INCLUDE_PATH} )
    include_directories( src/mpitools )

  add_executable(oaextendmpi EXCLUDE_FROM_ALL "utils/oaextend.cpp" ${headersextend} ${srcsextend} "src/mpitools/mpitools.cpp")
  target_link_libraries(oaextendmpi ${MPI_LIBRARIES})
  if(MPI_COMPILE_FLAGS)
    set_target_properties(oaextendmpi PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()
  if(MPI_LINK_FLAGS)
      if (WIN32)
		set_target_properties(oaextendmpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS} /MT")
		endif()
    set_target_properties(oaextendmpi PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
  #set_target_properties(oaextendmpi PROPERTIES COMPILE_FLAGS -I${MPIPATH})
  #TARGET_LINK_LIBRARIES(oaextendmpi /usr/lib/mpich2/lib/libmpichcxx)
  set_target_properties(oaextendmpi PROPERTIES COMPILE_DEFINITIONS "OAEXTEND_MULTICORE;NEWINTERFACE" )

#  list(APPEND extendprogs oaextendmpi) 
   if (WIN32)
    if(ZLIB_FOUND)
      TARGET_LINK_LIBRARIES(oaextendmpi ${ZLIB_LIBRARIES})
    endif()
   else()
      TARGET_LINK_LIBRARIES(oaextendmpi m ${ZLIB_LIBRARIES})
   endif()
  
endif(MPI_FOUND AND NOT WIN32 AND 1)

message(STATUS "Extend progs: ${extendprogs}")
endif()

message(STATUS "Progs: ${progs}")

###############
#message(STATUS "xxx: utils/oagraph.cpp ${headers} ${nautysrc}")
if (0)
ADD_EXECUTABLE(oagraph "utils/oagraph.cpp" ${headers} )  # ${nautysrc} 
  add_compile_flags(oagraph)
  TARGET_LINK_LIBRARIES(oagraph oalib)

  if (WIN32) 
      if(ZLIB_FOUND)
      TARGET_LINK_LIBRARIES(oagraph ${ZLIB_LIBRARIES})
    endif()
   #	message(STATUS "Win32 build: ${item}")
  else()
	TARGET_LINK_LIBRARIES(oagraph m)
  endif()
  if(ZLIB_FOUND)
	TARGET_LINK_LIBRARIES(oagraph ${ZLIB_LIBRARIES})
  endif()
endif()

###############

foreach(item ${progs})
  ADD_EXECUTABLE(${item} "utils/${item}.cpp" ${headers}  ) # ${srcs}
  add_compile_flags(${item})
  TARGET_LINK_LIBRARIES(${item} oalib)

  if (WIN32) 
    if(ZLIB_FOUND)
      TARGET_LINK_LIBRARIES(${item} ${ZLIB_LIBRARIES})
    endif()
  else()
	TARGET_LINK_LIBRARIES(${item} m)
  endif()
  if(ZLIB_FOUND)
	TARGET_LINK_LIBRARIES(${item} ${ZLIB_LIBRARIES})
  endif()
endforeach()

foreach(item ${extendprogs})
  if (WIN32) 
  else()
    set_target_properties(${item} PROPERTIES COMPILE_FLAGS "${COPTS}" )
	TARGET_LINK_LIBRARIES(${item} m)
  endif()
  if(ZLIB_FOUND)
	TARGET_LINK_LIBRARIES(${item} ${ZLIB_LIBRARIES})
  endif()

  set_target_properties(${item} PROPERTIES LINK_FLAGS "${LINKOPTS}" )
endforeach()

if(Boost_FOUND)
  foreach(item ${progs} ${extendprogs})
    target_link_libraries(${item} ${Boost_LIBRARIES})
  endforeach()
endif()


##### R test
site_name(HOSTNAME)
#message(STATUS "Hostname ${HOSTNAME}")

if("$ENV{USER}" STREQUAL "eendebakpt" AND "${HOSTNAME}" STREQUAL "marmot")
  set(headersR src/Deff.h src/arrayproperties.h src/arraytools.h src/strength.h src/mathtools.h src/oaoptions.h src/tools.h)
  set(srcsR  src/Deff.cpp src/arrayproperties.cpp src/arraytools.cpp src/strength.cpp src/mathtools.cpp src/tools.cpp src/oaoptions.cpp )
  
  add_executable(Rtest utils/Rtest.cpp ${srcsR} ${headersR})
  set_target_properties(Rtest PROPERTIES COMPILE_FLAGS "-DRPACKAGE  -Isrc -I/usr/share/R/include") # ${COPTS}
  #set_target_properties(Rtest PROPERTIES LINK_FLAGS "-lR")
  TARGET_LINK_LIBRARIES(Rtest R)
endif()

####### Octave / Matlab support ##########################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

find_package(Matlab)
message(STATUS "Matlab: ${MATLAB_FOUND}, Matlab_LIBRARIES ${Matlab_LIBRARIES}")
set(MATLAB_FOUND ${Matlab_FOUND})
set(MATLAB_LIBRARIES ${Matlab_LIBRARIES})

if(MATLAB_FOUND)
  message(STATUS ${srcs})
  set(Ditems Doptim)
  ADD_LIBRARY(Doptim EXCLUDE_FROM_ALL SHARED mex/Doptim.cpp)
  set_target_properties(${item} PROPERTIES COMPILE_FLAGS "-DMATLAB_MEX_FILE ${COPTS} ")

  ADD_EXECUTABLE(Doptimbin EXCLUDE_FROM_ALL mex/Doptim.cpp) #  ${srcs} ${srcsextend})
	      set_target_properties(${item} PROPERTIES COMPILE_FLAGS "-DMATLAB_MEX_FILE ${COPTS} -DMAINMEX")
  list(APPEND Ditems Doptimbin)
	
    message(STATUS "libs ${MATLAB_LIBRARIES}")

    foreach(item  ${Ditems})
      message(STATUS "adding ${item}")
	    if(WIN32)
	      TARGET_LINK_LIBRARIES(${item} m ${ZLIB_LIBRARIES} )
	    endif()
	    TARGET_LINK_LIBRARIES(${item} ${MATLAB_LIBRARIES})
	    TARGET_LINK_LIBRARIES(${item} oalib)

	    GET_TARGET_PROPERTY(itemtype ${item} TYPE)

	    if(WIN32)
		  set(MEXEXTX ".mexw${BITNESS}")
		  if (${itemtype} MATCHES "EXECUTABLE")
			  SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX "")
		  else()
		  if (CYGWIN)
			  SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX ${MEXEXTX})
		  else()
			  SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX ${MEXEXTX} LINK_FLAGS /export:mexFunction)
		  endif()

		  endif()

	    else(WIN32)
		     # set_target_properties(${item} PROPERTIES LINK_FLAGS "-s")	# strip
	    endif(WIN32)
	    
	    if (${itemtype} MATCHES "EXECUTABLE")
	    SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX "")
	  else()

	    #message(STATUS "${item} TYPE ${itemtype}: mex suffix")
	    if(WIN32)
	      SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX ".mexw${BITNESS}")
	    else()
	      SET_TARGET_PROPERTIES(${item} PROPERTIES SUFFIX ".mexa${BITNESS}")
	    endif()
	  endif()
	  set_target_properties(${item} PROPERTIES PREFIX "")
	endforeach()
endif()

# > matlab -nojvm -r "Doptim; exit"


####### Python and R support ##########################


#add_custom_target( Rswig )
#add_custom_target( R2 )
#
#add_custom_command(
#  TARGET Rswig POST_BUILD
#  COMMAND swig -c++ -DOADEV -Idev/ -Isrc/ -o oalib_wrapR.cpp -r oalib.i\; export PKG_CXXFLAGS="-Isrc/ -I./ -DSWIGR"\; echo "pre"
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#  DEPENDS 
#)
#
#add_custom_command(
#  TARGET R2 POST_BUILD
#  COMMAND rm -f oalibR.so\; R CMD SHLIB  -o oalibR.so oalib_wrapR.cpp src/pareto.cpp src/anyoption.cpp src/md5.cpp src/strength.cpp src/arrayproperties.cpp src/arraytools.cpp src/mathtools.cpp src/tools.cpp src/oaoptions.cpp src/bitarray/bit_array.cpp src/extend.cpp src/oadevelop.cpp src/lmc.cpp src/nonroot.cpp
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#  DEPENDS 
#)
#
#R CMD SHLIB  -o example.so oalib_wrapR.cpp liboalib.a-d
#Rscript example.r



#add_dependencies(pythoncheck python)
if(DOPYTHON)	
	add_custom_target( pythondoc )
	add_custom_command(
	  TARGET pythondoc POST_BUILD
	  COMMAND cd xml; xsltproc ../swig.xsl index.xml > swig_doc.i; cd ..
	  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

	add_custom_target( pinstall )
	add_custom_command(
	  TARGET pinstall
	  POST_BUILD
	  COMMAND python setup.py install --user
	#  COMMAND sudo python setup.py install
	  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	)

endif(DOPYTHON)


add_custom_target( doc )
add_custom_command(
  TARGET doc
  POST_BUILD
  COMMAND doxygen Doxyfile
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)


################ INSTALLATION OPTIONS ###########################

foreach(prog ${progs})
	install (TARGETS ${prog} DESTINATION bin)
endforeach()
foreach(prog ${extendprogs})
	install (TARGETS ${prog} DESTINATION bin)
endforeach()
install (TARGETS oalib DESTINATION lib)
install (FILES README.md FORMAT.txt LICENSE DESTINATION doc)
install (FILES oaconfig.txt DESTINATION bin)


SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Program to calculate orthogonal arrays")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")


SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "3")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

#SET(CPACK_STRIP_FILES "bin/MyExecutable")


#set(CPACK_SOURCE_PACKAGE_FILE_NAME "xxx-${VERSION}" CACHE INTERNAL "tarball basename")
set(CPACK_SOURCE_GENERATOR TGZ)

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-source")

set(OADEVPACK 0)
set(CPACK_SOURCE_IGNORE_FILES "tests.txt$;/vsc/;/old/;/Eigen3.2/;/Eigen-devel/;/Eigen3.1/;/buildvs/;a.out$;/TestingCommands;/testing/;tags$;_oalib.so$;Makefile$;setup.py~$;ChangeLog~$;CMakeFiles$;TODO$;oaconfig-speedtest.txt$;CMakeCache.txt;cmake_install.cmake;.pyc$;/oacode.kdev4;/OA.kde;.cpp~$;/CMakeLists-develop.txt;/build/;/matlabtools/;/xml/;/pythonscripts/;/miscx/;/performancetest/;/oldx/;/buildwin/;/html/;/testsuite/;/.svn/;${CPACK_SOURCE_IGNORE_FILES}")
set(CPACK_SOURCE_IGNORE_FILES "example.so;oatest2.cpp;/build2/;PrecompiledHeader.cmake;/Eigen3.2.1/" ${CPACK_SOURCE_IGNORE_FILES})
set(CPACK_SOURCE_IGNORE_FILES "_oalib.cpython-34m.so;mex/Doptim.mex;/.git/;x.txt;/appveyor/;/tmp/;/dist/;/build2/;/hopper/;/OAlib.egg-info/;/buildx/;/OApackage.egg-info/;x2.txt;appveyor.yml;oalib_wrap_dev.cxx;oalib_wrap.cpp;oalib_wrap.cxx;sandbox;" ${CPACK_SOURCE_IGNORE_FILES})
set(CPACK_SOURCE_IGNORE_FILES "\\\\.cpp~$" "\\\\.txt~$" ${CPACK_SOURCE_IGNORE_FILES})
set(CPACK_SOURCE_IGNORE_FILES "~$" ${CPACK_SOURCE_IGNORE_FILES})

if (VERBOSE)
  message(STATUS "CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES}")
endif()

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-extend-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-source")
message(STATUS "Building package with extend options")
if(OADEVPACK)
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-develop-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}-source")
	message(STATUS "Building development package!")
else()
  set(CPACK_SOURCE_IGNORE_FILES "/testdata/;/pythondevelop/;${CPACK_SOURCE_IGNORE_FILES}")
  set(CPACK_SOURCE_IGNORE_FILES "oaextendA.cpp$;oatest.cpp$;oa_depth_extend.cpp$;/dev/oadevelop.h$;/dev/oadevelop.cpp;${CPACK_SOURCE_IGNORE_FILES}")
  #message(STATUS "CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES}")
endif(OADEVPACK)
#message(STATUS "CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES}")

# cmake package type: CygwinBinary, DEB, TGZ, ZIP, etc.
if(WIN32)
  SET(CPACK_GENERATOR "ZIP")
else()
  SET(CPACK_GENERATOR "DEB")
  #SET(CPACK_GENERATOR "TGZ")
  SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Pieter Eendebak") #required
endif()
INCLUDE(CPack)

add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)


### for VCS
#module load SWIG/2.0.4-ictce-4.0.1-Python-2.7.3
#module load zlib/1.2.5-ictce-4.0.1 
#module load CMake/2.8.0-ictce-3.2.1.015.u1   # negeer de ERROR die je krijgt
#export CC=icc
#export CXX=icpc

# cmake -DCMAKE_INSTALL_PREFIX:PATH=~vsc20149 ..

# or
# module load GCC/4.7.0
# cmake -DCMAKE_INSTALL_PREFIX:PATH=~vsc20149 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

# mmlsquota

