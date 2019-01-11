add_executable(trf2proclu-ngs.exe)
target_link_libraries(trf2proclu-ngs.exe easylife m)
target_sources(trf2proclu-ngs.exe
    PRIVATE trf2proclu-ngs.c
)
install(TARGETS trf2proclu-ngs.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
