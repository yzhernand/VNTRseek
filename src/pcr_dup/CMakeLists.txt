add_executable(pcr_dup.exe)
target_link_libraries(pcr_dup.exe m)
target_sources(pcr_dup.exe
    PRIVATE main.c
    PRIVATE LinkedList.c
    PRIVATE "bitwise LCS single word.c"
    PRIVATE "bitwise LCS multiple word.c"
)

install(TARGETS pcr_dup.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
