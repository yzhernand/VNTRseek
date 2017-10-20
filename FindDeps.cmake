# Check for ZLIB
FIND_PACKAGE(ZLIB)
IF (ZLIB_FOUND)
    MESSAGE("ZLIB lib: ${ZLIB_LIBRARY}")
    MESSAGE("ZLIB include: ${ZLIB_INCLUDE_DIR}")
ELSE (ZLIB_FOUND)
    MESSAGE(FATAL_ERROR "ZLIB not found. Please install zlib to continue.")
ENDIF(ZLIB_FOUND)

# Check Perl minimum version
FUNCTION(PERL_REQ_VER VER)
    FIND_PACKAGE(Perl)
    IF (PERL_FOUND)
        IF (PERL_VERSION_STRING VERSION_GREATER ${VER} OR PERL_VERSION_STRING VERSION_EQUAL ${VER})
            MESSAGE(STATUS "Perl >= ${VER} (${PERL_VERSION_STRING})")
        ELSE ()
            MESSAGE(FATAL_ERROR "Error: Perl version must be >= ${VER}. You have ${PERL_VERSION_STRING}.")
        ENDIF ()
    ELSE (PERL_FOUND)
        MESSAGE(FATAL_ERROR "Error: ${PACKAGE_NAME} requires Perl version ${VER} in order to run.")
    ENDIF(PERL_FOUND)
ENDFUNCTION(PERL_REQ_VER)

# Check GCC minimum version
FUNCTION(GCC_REQ_VER VER)
    MESSAGE(STATUS "Checking GCC version...")
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION)
    STRING(REGEX REPLACE "(\r?\n)+$" "" GCC_VERSION "${GCC_VERSION}")

    IF(GCC_VERSION VERSION_GREATER ${VER} OR GCC_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "GCC version >= ${VER} (${GCC_VERSION})")
    ELSE()
        MESSAGE(FATAL_ERROR "Error: GCC version must be >= ${VER}. You have ${GCC_VERSION}")
    ENDIF()
ENDFUNCTION(GCC_REQ_VER)

# Check GLIBC minimum version
FUNCTION(GLIBC_REQ_VER VER)
    MESSAGE(STATUS "Checking GLIBC version...")
    EXECUTE_PROCESS(COMMAND /lib/libc.so.6
        OUTPUT_VARIABLE GLIBC_VERSION)
    STRING(REGEX MATCH "[0-9]\\.[0-9][0-9]?" GLIBC_VERSION "${GLIBC_VERSION}")
    MESSAGE(STATUS "GLIBC version: ${GLIBC_VERSION}")

    IF(GLIBC_VERSION VERSION_LESS ${VER} OR GLIBC_VERSION VERSION_EQUAL ${VER})
        MESSAGE(STATUS "Downloading legacy build of TRF...")
        SET(TRFDLName "trf${TRFVer}.legacy${ARCH}" PARENT_SCOPE)
    ENDIF()
ENDFUNCTION(GLIBC_REQ_VER)

# Check mysql binary is found, and is of correct version
FUNCTION(MYSQL_REQ_VER VER)
    MESSAGE(STATUS "Checking MySQL client version...")
    FIND_PROGRAM(MYSQL_LOC "mysql")
    FIND_PROGRAM(MYSQL_CONF_LOC "mysql_config")

    IF(MYSQL_CONF_LOC STREQUAL "MYSQL_CONF_LOC-NOTFOUND")
        MESSAGE(FATAL_ERROR "Error: Unable to detect your version of mysql. Please ensure that mysql_config is in your PATH.")
    ELSEIF(MYSQL_LOC STREQUAL "MYSQL_LOC-NOTFOUND")
        MESSAGE(FATAL_ERROR "The MySQL client is a requirement for using the TRF pipeline. You can install it for free at https://www.mysql.com/ or through your OS's distribution channel.")
    ELSE()
        EXECUTE_PROCESS(COMMAND ${MYSQL_CONF_LOC} --version
            OUTPUT_VARIABLE MYSQL_VERSION)
        STRING(REGEX REPLACE "(\r?\n)+$" "" MYSQL_VERSION "${MYSQL_VERSION}")

        IF(MYSQL_VERSION VERSION_GREATER ${VER} OR MYSQL_VERSION VERSION_EQUAL ${VER})
            MESSAGE(STATUS "MySQL client version >= ${VER} (${MYSQL_VERSION})")
        ELSE()
            MESSAGE(FATAL_ERROR "Error: MySQL version must be >= ${VER}. You have ${MYSQL_VERSION}")
        ENDIF()
    ENDIF()
ENDFUNCTION(MYSQL_REQ_VER)

# Check for seqtk


# Check for samtools


# Check for bedtools

