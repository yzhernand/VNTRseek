# source files for psearch
add_executable(psearch.exe)
target_link_libraries(psearch.exe easylife m)
target_sources(psearch.exe
    PRIVATE psearch.c
    PRIVATE bitwise\ edit\ distance\ alignment\ multiple\ word\ no\ end\ penalty.c
    PRIVATE bitwise\ LCS\ single\ word.c
    PRIVATE bitwise\ LCS\ multiple\ word.c
)

install(TARGETS psearch.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
