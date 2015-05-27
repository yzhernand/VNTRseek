FILE(GLOB shellscripts RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.sh")

# Install shell scripts.
INSTALL(PROGRAMS ${shellscripts}
    DESTINATION ${InstallSuffix}
    )
