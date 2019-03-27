AUTOMAKE_OPTIONS = subdir-objects

# the program to build
bin_PROGRAMS = 	pourRNA
				
SUFFIXES = _cmdl.c _cmdl.h .ggo

GENGETOPT_CMDL =  pourRNA_cmdl.c pourRNA_cmdl.h
GENGETOPT_FILES =  pourRNA.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL} PriorityQueue.cpp

# remove all generated content
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
.ggo_cmdl.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdl.c_cmdl.h :

# pourRNA sources
pourRNA_SOURCES =  ${GENGETOPT_CMDL} \
                   SpookyHash/SpookyV2.h \
                   SpookyHash/SpookyV2.cpp \
                   SpookyHash/SpookyHashMap.h \
                   BIUlibPart/LimitedHash.hh \
                   BIUlibPart/LimitedHash.icc \
                   BIUlibPart/assertbiu.hh \
                   BIUlibPart/MatrixSparse.hh \
                   BIUlibPart/MatrixSparse.icc \
                   RNAkinetics/RateMatrixUtil.h \
                   RNAkinetics/RateMatrixUtil.cpp \
                   RNAkinetics/RNA_NeighMinFilter.h \
                   RNAkinetics/RNA_NeighMinFilter.cpp \
                   TypeID.h \
                   TypeID.cpp \
                   Concurrent_Queue.h \
                   MyState.h \
                   MyState.cpp \
                   PriorityQueue.h \
                   Flooder.h \
                   Flooder.cpp \
                   WalkGradientHashed.h \
                   WalkGradientHashed.cpp \
                   StructureUtils.h \
                   StructureUtils.cpp \
                   PairHashTable.h \
                   Concurrent_Pair_Hash_Map.h \
                   PairHashMap.h \
                   SC_PartitionFunction.h \
                   SC_PartitionFunction.cpp \
                   SC_DotPlot.h \
                   SC_DotPlot.cpp \
                   StateCollector.h \
                   StateCollector.cpp \
                   StatePairCollector.h \
                   StatePairCollector.cpp \
                   pourRNA_cmdl.h \
                   pourRNA.cpp 




if HAVE_HELP2MAN
  man1_MANS = $(ax_help2man_MANS)
  @AX_HELP2MAN_RULES@
endif



