.PHONY: all libs semphy programs clean install

all: libs programs

debug: libs.debug

%: libs.% programs.%
	echo $@

libs: libs.all

programs: programs.all

programs.all: libs
programs.debug: libs.debug

semphy: programs.semphy

install: programs.install

programs.install programs.all semphy: libs

clean: libs.clean programs.clean

libs.%:
	+cd libs;make $(*)

programs.%:
	+cd programs;make $(*)

tags: libs/*/*.cpp libs/*/*.h programs/*/*.h programs/*/*.cpp
	etags --members --language=c++ $^ 

