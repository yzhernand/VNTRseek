SET(REFFLANKALIGN_SRCS
    refflankalign.c narrowbandDistanceAlignment.c
    easylife.c
    )
ADD_EXECUTABLE(refflankalign.exe ${REFFLANKALIGN_SRCS})
TARGET_LINK_LIBRARIES(refflankalign.exe m)

INSTALL(TARGETS refflankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
