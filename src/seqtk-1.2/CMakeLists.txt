set(SEQTK_SRCS
    seqtk.c	khash.h	kseq.h
    )
add_executable(seqtk ${SEQTK_SRCS})
target_link_libraries(seqtk z m)
target_compile_options(seqtk PRIVATE -Wno-unused-function)

install(TARGETS seqtk
    RUNTIME DESTINATION ${InstallSuffix}
    )
