#	$Id: Makefile 11987 2014-01-30 10:23:04Z haim $	

# this split is vital becouse of a bug in make 3.80.1 - see 
# http://www.cygwin.com/ml/cygwin/2004-09/msg01659.html

PROGRAMS1= fastml gainLoss
PROGRAMS2= indelCoder
PROGRAMS = $(PROGRAMS1) $(PROGRAMS2)

# all has to be the FIRST task!
TASKS= all clean test depend debug All install doubleRep
.PHONY: $(TASKS) $(PROGRAMS)

 define TASKS_template1
    $(1): $$(addsuffix .$(1),$(PROGRAMS1))
 endef

 define TASKS_template2
    $(1): $$(addsuffix .$(1),$(PROGRAMS2))
 endef

$(foreach task,$(TASKS),$(eval $(call TASKS_template1,$(task))))
$(foreach task,$(TASKS),$(eval $(call TASKS_template2,$(task))))

define PROGRAM_template
   $(1).%:
	+cd $(1) && make $$(*)
endef

$(foreach prog,$(PROGRAMS),$(eval $(call PROGRAM_template,$(prog))))



$(PROGRAMS):
	+cd $@ && make 

