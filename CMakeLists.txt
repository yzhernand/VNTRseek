SET(PACKAGE VNTRseek)
SET(PACKAGE_NAME VNTRseek)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(VNTRseek)

# Set current VNTRseek version. Not user configurable
SET(VNTRVer "1.09")
SET(InstallSuffix "vntrseek${VNTRVer}")

# Set current TRF version. Not user configurable
SET(TRFVer "409")

## Set install directory
#IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#    SET(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix (eg, /opt, the default). All non-web components of the pipeline will be installed in a directory under this, including binaries and support scripts or files." FORCE)
#ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Set the build type -- Does nothing at the moment
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Flags for the different builds. So far, we only have a definition for flags
# in release mode
SET(CMAKE_C_FLAGS_RELEASE "-O2")
# Dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

## Detect arch and OS
INCLUDE("DetermineArch.cmake")
DETERMINE_ARCH(ARCH)

IF (ARCH STREQUAL "UNSUPPORTED")
    MESSAGE(FATAL_ERROR "Your current OS and architecture (${CMAKE_SYSTEM_NAME}, ${CMAKE_SYSTEM_PROCESSOR}) are unsupported. Installation cannot continue.")
ELSEIF (ARCH STREQUAL "i386")
    MESSAGE(FATAL_ERROR "32-bit architecture is not supported. Installation cannot continue.")
ENDIF()

## Check for required installed programs and versions
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})
SET(TRFDLName "trf${TRFVer}.${ARCH}")
SET(GLIB2_LEGACY_VERSION "2.14")
INCLUDE("FindDeps.cmake")
GCC_REQ_VER(4.1.2)
GLIBC_REQ_VER(2.14)
MYSQL_REQ_VER(5.0.95)
PERL_REQ_VER(5.8.8)

MESSAGE(STATUS "Your processor is ${CMAKE_SYSTEM_PROCESSOR} and you are running ${CMAKE_SYSTEM_NAME}. This means we'll download ${TRFDLName}")
# May need to change if other platforms are supported
#SET(TRFBin "${TRFDLName}.exe")
SET(TRFBin "trf${TRFVer}-ngs.linux.exe")

## Check for required Perl modules
FIND_PACKAGE(PerlModules COMPONENTS DBI DBD::mysql)
IF(NOT PERLMODULES_FOUND)
    MESSAGE(FATAL_ERROR "Not all required perl modules were found on your system")
ENDIF(NOT PERLMODULES_FOUND)

# Build files in src
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(perl)
ADD_SUBDIRECTORY(scripts)
#ADD_SUBDIRECTORY(etc)

# Post install message
ADD_SUBDIRECTORY(postinst)

# Install config and other files
INSTALL(CODE "
    if (EXISTS \"${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vs.cnf\")
        MESSAGE(STATUS \"vs.cnf file exists. Will install new version as 'vs.cnf.new'\")
        file(INSTALL \"${CMAKE_SOURCE_DIR}/vs.cnf\"
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}
            RENAME \"vs.cnf.new\"
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
        )
    else()
        file(INSTALL \"${CMAKE_SOURCE_DIR}/vs.cnf\"
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}
            PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ
        )
    endif(EXISTS \"${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vs.cnf\")
")
SET(MISC_FILES
    HISTORY
    COPYING CLUSTER
    README RUNNING
    sqlite_schema.sql
    )
INSTALL(FILES ${MISC_FILES}
    DESTINATION ${InstallSuffix}
    )

# Download and install trf
FILE(DOWNLOAD http://tandem.bu.edu/trf/downloads/${TRFDLName}
    "${CMAKE_CURRENT_BINARY_DIR}/${TRFBin}"
    TIMEOUT 30
    STATUS TRFdlStat
    SHOW_PROGRESS)
LIST(GET TRFdlStat 0 TRFdlErrNo)
IF (TRFdlErrNo)
    MESSAGE(WARNING "Error encountered downloading TRF. Download it manually from: http://tandem.bu.edu/trf/trf.download.html")
ELSE (TRFdlErrNo)
    INSTALL(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/${TRFBin}" DESTINATION ${InstallSuffix})
ENDIF (TRFdlErrNo)
