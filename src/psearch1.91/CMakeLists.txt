# source files for psearch
SET(PSEARCH_SRCS
    psearch.c            bitwise\ edit\ distance\ alignment\ multiple\ word\ no\ end\ penalty.c
    bitwise\ LCS\ single\ word.c
    bitwise\ LCS\ multiple\ word.c
    )
ADD_EXECUTABLE(psearch.exe ${PSEARCH_SRCS})
TARGET_LINK_LIBRARIES(psearch.exe easylife m)

INSTALL(TARGETS psearch.exe
    RUNTIME DESTINATION ${InstallSuffix}
    )
