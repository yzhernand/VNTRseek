# source files for redund.exe
SET(REDUND_SRCS
    redund.c
    )
ADD_EXECUTABLE(redund.exe ${REDUND_SRCS})
TARGET_LINK_LIBRARIES(redund.exe LINK_PUBLIC m)

INSTALL(TARGETS redund.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
