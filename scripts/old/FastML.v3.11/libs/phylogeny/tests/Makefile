CPPFLAGS=      -g -Wall -Wno-sign-compare -I.. -DLOG -ftemplate-depth-25 -O0
CPPFLAGSDEBUG= -g -Wall -Wno-sign-compare -I.. -DLOG -ftemplate-depth-32

# -O3
LDFLAGS= -L.. 
LDLIBS= -lEvolTree

CC=g++
CXX=g++

.PHONY: test tests run all clean

TESTS= split_test splitMap_test bootstrap_test				\
	given2seqTheDistBetweenThem-3					\
	given2seqEstimateTheDistBetweenThem computeNJtreeJCmodel	\
	checkTreeLikelihoodGivenBranches optimizeBranchesJC_EM		\
	given2seqEstimateTheDistBetweenThemGamma			\
	optimizeBranchesJC_EM_gamma					\
	optimizeBranchesJC_EM_gam_estimate_alp				\
	given2seqEstimateTheDistanceK2P optimize_HKY_param		\
	doubleRep readTreeWithComments
#ludouble_test 	DistanceTableFromTree 
OTHER_TESTS= exhaustiveSearch

all: test
run:  $(addsuffix .out.tmp,$(TESTS)) 

libEvolDebug=../libEvolTreeDebug.a

DEBUGEXEC = $(addsuffix .debug,$(TESTS))
$(DEBUGEXEC): $(libEvolDebug)

$(TESTS): ../libEvolTree.a

test: $(addsuffix .test,$(TESTS)) 

alltest: $(addsuffix .test,$(TESTS) $(OTHER_TESTS)) 

%.test: %.out.tmp %.out.standard
	diff $^

readTreeWithComments.out.tmp: readTreeWithComments treeWithComments.tree
	./$^>$@

%.out.tmp: %
	$(*) > $@

tests: $(TESTS)


%.debug.o: %.c
	$(CC) -c $(CPPFLAGSDEBUG) $(CFLAGS) $< -o $@

%.debug.o: %.cpp
	$(CXX) -c $(CPPFLAGSDEBUG) $(CXXFLAGS) $< -o $@

debug: $(DEBUGEXEC)


clean:
	-rm -f $(TESTS) *.out.tmp *.o
