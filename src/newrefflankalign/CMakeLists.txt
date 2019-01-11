add_executable(refflankalign.exe)
target_link_libraries(refflankalign.exe easylife m)
target_sources(refflankalign.exe
    PRIVATE refflankalign.c
    PRIVATE narrowbandDistanceAlignment.c
)
install(TARGETS refflankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
