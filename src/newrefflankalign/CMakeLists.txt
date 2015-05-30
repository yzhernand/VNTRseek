SET(REFFLANKALIGN_SRCS
    refflankalign.c narrowbandDistanceAlignment.c
    )
ADD_EXECUTABLE(refflankalign.exe ${REFFLANKALIGN_SRCS})
TARGET_LINK_LIBRARIES(refflankalign.exe PUBLIC m)

INSTALL(TARGETS refflankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
