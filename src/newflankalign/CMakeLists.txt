SET(FLANKALIGN_SRCS
    flankalign.c
    "bitwise edit distance alignment multiple word no end penalty.c"
    )
ADD_EXECUTABLE(flankalign.exe ${FLANKALIGN_SRCS})
TARGET_LINK_LIBRARIES(flankalign.exe PUBLIC m)

INSTALL(TARGETS flankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
