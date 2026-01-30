cmake_minimum_required(VERSION 3.20)
project(OApackage LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Options
option(BUILD_PYTHON_MODULE "Build Python SWIG module" ON)
option(BUILD_EXECUTABLES "Build command-line executables" ON)
option(OADEV "Build development code" OFF)

# ----------------------------------------------------------------------------
# Find packages
# ----------------------------------------------------------------------------
find_package(ZLIB)
if(ZLIB_FOUND)
    message(STATUS "ZLIB found: ${ZLIB_LIBRARIES}")
else()
    message(STATUS "ZLIB not found - building without compression support")
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# ----------------------------------------------------------------------------
# Source files
# ----------------------------------------------------------------------------
set(OA_HEADERS
    src/pareto.h
    src/unittests.h
    src/anyoption.h
    src/arrayproperties.h
    src/arraytools.h
    src/graphtools.h
    src/mathtools.h
    src/oaoptions.h
    src/tools.h
    src/md5.h
    src/Deff.h
    src/evenodd.h
    src/version.h
    src/conference.h
    src/extend.h
    src/strength.h
    src/lmc.h
    src/nonroot.h
)

set(OA_SOURCES
    src/pareto.cpp
    src/unittests.cpp
    src/md5.cpp
    src/strength.cpp
    src/arrayproperties.cpp
    src/graphtools.cpp
    src/arraytools.cpp
    src/mathtools.cpp
    src/tools.cpp
    src/oaoptions.cpp
    src/bitarray/bit_array.cpp
    src/Deff.cpp
    src/evenodd.cpp
    src/conference.cpp
    src/extend.cpp
    src/lmc.cpp
    src/nonroot.cpp
)

set(NAUTY_SOURCES
    src/nauty/nauty.c
    src/nauty/nautinv.c
    src/nauty/nautil.c
    src/nauty/naurng.c
    src/nauty/naugraph.c
    src/nauty/schreier.c
    src/nauty/naugroup.c
)

# ----------------------------------------------------------------------------
# Compiler flags
# ----------------------------------------------------------------------------
set(OA_COMPILE_DEFINITIONS FULLPACKAGE)

if(ZLIB_FOUND)
    list(APPEND OA_COMPILE_DEFINITIONS USEZLIB)
else()
    list(APPEND OA_COMPILE_DEFINITIONS NOZLIB)
endif()

if(WIN32)
    list(APPEND OA_COMPILE_DEFINITIONS WIN32 _WIN32)
endif()

if(OADEV)
    list(APPEND OA_COMPILE_DEFINITIONS OADEV)
endif()

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -O3
        -Wno-unknown-pragmas
        -Wno-sign-compare
        -Wno-unused-variable
        -Wno-unused-result
    )
    if(NOT WIN32)
        add_compile_options(-fPIC)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        add_compile_options(-Wno-c++11-compat-deprecated-writable-strings)
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_compile_options(-Wno-return-type -Wno-date-time)
    endif()
elseif(MSVC)
    add_compile_options(/EHsc /W0 /bigobj /wd4996 /wd4244 /wd26451 /wd26495 /wd6255 /D_SCL_SECURE_NO_WARNINGS)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()



# ----------------------------------------------------------------------------
# Python SWIG module
# ----------------------------------------------------------------------------
if(BUILD_PYTHON_MODULE)
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
    find_package(SWIG 4.0 REQUIRED COMPONENTS python)

    message(STATUS "Python executable: ${Python_EXECUTABLE}")
    message(STATUS "Python include: ${Python_INCLUDE_DIRS}")
    message(STATUS "NumPy include: ${Python_NumPy_INCLUDE_DIRS}")
    message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
    message(STATUS "SWIG version: ${SWIG_VERSION}")

    include(${SWIG_USE_FILE})

    # Set SWIG properties
    set_property(SOURCE oalib.i PROPERTY CPLUSPLUS ON)
    set_property(SOURCE oalib.i PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)

    # SWIG flags
    set(SWIG_MODULE_FLAGS
        -c++
        -doxygen
        -w503,401,362,302,389,446,509,305
        -DSWIGCODE
        -DFULLPACKAGE
    )
    if(ZLIB_FOUND)
        list(APPEND SWIG_MODULE_FLAGS -DUSEZLIB)
    else()
        list(APPEND SWIG_MODULE_FLAGS -DNOZLIB)
    endif()
    if(WIN32)
        list(APPEND SWIG_MODULE_FLAGS -DWIN32 -D_WIN32)
    endif()
    list(APPEND SWIG_MODULE_FLAGS -DNOOMP)

    set_property(SOURCE oalib.i PROPERTY SWIG_FLAGS ${SWIG_MODULE_FLAGS})

    # Create SWIG module
    swig_add_library(oalib
        TYPE MODULE
        LANGUAGE python
        SOURCES oalib.i ${OA_SOURCES} ${NAUTY_SOURCES}
    )

    # Include directories
    target_include_directories(oalib PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nauty
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Eigen
        ${Python_INCLUDE_DIRS}
        ${Python_NumPy_INCLUDE_DIRS}
    )

    # Compile definitions
    target_compile_definitions(oalib PRIVATE
        ${OA_COMPILE_DEFINITIONS}
        SWIGCODE
        NOOMP
    )

    # Link libraries
    target_link_libraries(oalib PRIVATE Python::Module Python::NumPy)
    if(ZLIB_FOUND)
        target_link_libraries(oalib PRIVATE ZLIB::ZLIB)
    endif()
    if(NOT WIN32)
        target_link_libraries(oalib PRIVATE m)
    endif()

    # Set output name
    set_target_properties(oalib PROPERTIES
        OUTPUT_NAME "_oalib"
        PREFIX ""
    )
    if(WIN32)
        set_target_properties(oalib PROPERTIES SUFFIX ".pyd")
    endif()

    # Install targets for scikit-build-core
    install(TARGETS oalib LIBRARY DESTINATION .)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/oalib.py DESTINATION .)
endif()

# ----------------------------------------------------------------------------
# Optional: Command-line executables
# ----------------------------------------------------------------------------
if(BUILD_EXECUTABLES)
    # Additional sources for executables
    set(OA_EXEC_SOURCES
        ${OA_SOURCES}
        src/anyoption.cpp
        src/InfInt.cpp
        src/depth_extend.cpp
    )

    # Static library for executables
    add_library(oalib_static STATIC ${OA_EXEC_SOURCES} ${NAUTY_SOURCES})
    target_include_directories(oalib_static PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/nauty
    )
    target_compile_definitions(oalib_static PUBLIC ${OA_COMPILE_DEFINITIONS})
    if(ZLIB_FOUND)
        target_link_libraries(oalib_static PUBLIC ZLIB::ZLIB)
    endif()
    if(NOT WIN32)
        target_link_libraries(oalib_static PUBLIC m)
    endif()

    # List of executables
    set(PROGS
        oacat oajoin oapareto oasplit oaanalyse oainfo oaunittest
        oaconvert oafilter oareduceDOP pareto_example oaconference
        oaconfcheck oacheck oastreaming oaranktest oatest oaclustergather
        oa_depth_extend oa_select_maxj
    )

    foreach(prog ${PROGS})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utils/${prog}.cpp")
            add_executable(${prog} utils/${prog}.cpp)
            target_link_libraries(${prog} PRIVATE oalib_static)

            target_precompile_headers(${prog} PRIVATE
                src/oaconfig.h
                src/arraytools.h
                <Eigen/Core>
            )
            message(STATUS "Adding executable: ${prog}")

        endif()
    endforeach()

    # oaextendsingle
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utils/oaextend.cpp")
        add_executable(oaextendsingle utils/oaextend.cpp)
        target_link_libraries(oaextendsingle PRIVATE oalib_static)
        target_compile_definitions(oaextendsingle PRIVATE OAEXTEND_SINGLECORE NEWINTERFACE)
        message(STATUS "Adding executable: oaextendsingle")
    endif()
endif()
