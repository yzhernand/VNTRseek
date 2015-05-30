SET(JOINC_SRCS
    joinc.c
    )
ADD_EXECUTABLE(join_clusters.exe ${JOINC_SRCS})
TARGET_LINK_LIBRARIES(join_clusters.exe PUBLIC m)

INSTALL(TARGETS join_clusters.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
