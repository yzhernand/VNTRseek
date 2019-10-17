set(FLANKALIGN_SRCS
	flankalign.c
    "bitwise edit distance alignment multiple word no end penalty.c"
    )

add_executable(flankalign.exe ${FLANKALIGN_SRCS})
target_sources(flankalign.exe
    PRIVATE ${FLANKALIGN_SRCS}
)
target_link_libraries(flankalign.exe easylife m)

install(TARGETS flankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
