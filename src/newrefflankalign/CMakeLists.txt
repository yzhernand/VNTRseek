set(REFFLANKALIGN_SRCS
	refflankalign.c
	narrowbandDistanceAlignment.c
	)

add_executable(refflankalign.exe ${REFFLANKALIGN_SRCS})
target_link_libraries(refflankalign.exe easylife m)
target_sources(refflankalign.exe
    PRIVATE ${REFFLANKALIGN_SRCS}
)
install(TARGETS refflankalign.exe
    RUNTIME DESTINATION ${InstallSuffix}
)
