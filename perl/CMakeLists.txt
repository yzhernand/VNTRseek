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
# if(CMAKE_HOST_UNIX)
#     INSTALL(CODE "
#     EXECUTE_PROCESS(COMMAND ln -sf ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl vntrseek
#         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#        )
#     ")
#     INSTALL(FILES "vntrseek"
#         DESTINATION bin
#         OPTIONAL
#         )
#     # INSTALL(CODE "execute_process( \
#     #         COMMAND ${CMAKE_COMMAND} -E create_symlink \
#     #         ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl \
#     #         ${CMAKE_INSTALL_PREFIX}/bin/vntrseek   \
#     #     )"
#     #     ERROR_VARIABLE VNTRSEEK_SYMLINK_ERROR
#     # )
#     # IF(VNTRSEEK_SYMLINK_ERROR)
#     #     MESSAGE(WARNING "Problem creating symlink for vntrseek command. \
#     #         If you'd like the command 'vntrseek' to be in your PATH either \
#     #         add the directory '${CMAKE_INSTALL_PREFIX}/${InstallSuffix}' to \
#     #         your PATH or create a symbolic link to \
#     #         ${CMAKE_INSTALL_PREFIX}/${InstallSuffix}/vntrseek.pl named \
#     #         'vntrseek' in a directory in your path (eg, /usr/local/bin).")
#     # ENDIF()
# endif(CMAKE_HOST_UNIX)
