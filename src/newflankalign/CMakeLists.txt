SET(FLANKALIGN_SRCS
    flankalign.c    easylife.c
    "bitwise edit distance alignment multiple word no end penalty.c"
    )
ADD_EXECUTABLE(flankalign.exe ${FLANKALIGN_SRCS})
TARGET_LINK_LIBRARIES(flankalign.exe m)

INSTALL(TARGETS flankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
