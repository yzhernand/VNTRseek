set(PCR_DUP_SRCS
	main.c
	LinkedList.c
	"bitwise LCS single word.c"
	"bitwise LCS multiple word.c"
	)

add_executable(pcr_dup.exe ${PCR_DUP_SRCS})
target_link_libraries(pcr_dup.exe m)
target_sources(pcr_dup.exe
    PRIVATE ${PCR_DUP_SRCS}
)

install(TARGETS pcr_dup.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
