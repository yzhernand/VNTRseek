FILE(GLOB perlscripts RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.pl")

# Install perl scripts.
INSTALL(PROGRAMS ${perlscripts}
    DESTINATION ${InstallSuffix}
    )

# Install vntr perl modules
INSTALL(DIRECTORY vntr/
    DESTINATION ${InstallSuffix}/vntr
    )

# Create symlink to vntrseek.pl in bin
if(CMAKE_HOST_UNIX)
    INSTALL(CODE "
    EXECUTE_PROCESS(COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl vntrseek
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
       )
    ")
    INSTALL(FILES "vntrseek"
        DESTINATION bin
        OPTIONAL
        )
endif(CMAKE_HOST_UNIX)
