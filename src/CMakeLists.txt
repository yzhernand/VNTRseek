# Install euclidean distance binary file
INSTALL(FILES eucledian.dst
    DESTINATION ${InstallSuffix}
    )

ADD_SUBDIRECTORY(pcr_dup)
ADD_SUBDIRECTORY(redund_code)
ADD_SUBDIRECTORY(newflankalign)
#vntr
ADD_SUBDIRECTORY(join_clusters)
ADD_SUBDIRECTORY(trf2proclu-ngs)
ADD_SUBDIRECTORY(psearch1.91)
ADD_SUBDIRECTORY(newrefflankalign)
