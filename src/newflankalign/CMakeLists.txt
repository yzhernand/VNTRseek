add_executable(flankalign.exe)
target_sources(flankalign.exe
    PRIVATE flankalign.c
    PRIVATE "bitwise edit distance alignment multiple word no end penalty.c"
)
target_link_libraries(flankalign.exe easylife m)

install(TARGETS flankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
