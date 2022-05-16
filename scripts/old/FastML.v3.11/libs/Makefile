#	$Id: Makefile 942 2006-10-18 12:28:12Z ninio $	

# There might be need for a split (as done in programs/Makefile) becouse of a bug in make 3.80.1 - see 
# http://www.cygwin.com/ml/cygwin/2004-09/msg01659.html

LIBS= phylogeny

# all has to be the FIRST task!
TASKS= all clean test depend debug All install doubleRep
.PHONY: $(TASKS) $(LIBS)

 define TASKS_template
    $(1): $$(addsuffix .$(1),$(LIBS))
 endef

$(foreach task,$(TASKS),$(eval $(call TASKS_template,$(task))))

define LIB_template
   $(1).%:
	+cd $(1) && make $$(*)
endef

$(foreach lib,$(LIBS),$(eval $(call LIB_template,$(lib))))



$(LIBS):
	+cd $@ && make 

