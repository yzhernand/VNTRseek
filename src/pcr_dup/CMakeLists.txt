SET(PCRDUP_SRCS
    main.c   LinkedList.c
    "bitwise LCS single word.c" "bitwise LCS multiple word.c"
    )
ADD_EXECUTABLE(pcr_dup.exe ${PCRDUP_SRCS})
TARGET_LINK_LIBRARIES(pcr_dup.exe m)

INSTALL(TARGETS pcr_dup.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
