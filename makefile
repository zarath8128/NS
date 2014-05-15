override OPT_FLAGS+=
override FLAGS += -Wall -Wextra ${OPT_FLAGS}
override CFLAGS+=${FLAGS}
override CXXFLAGS+=${FLAGS} -std=c++1y
ALL=test

.PHONY:all clean syntax debug release profile optimize reset

all:${ALL} 
clean:
	rm -rf ${ALL} *.o *.out *.gcda *.profile *.optimized
syntax:
	@${MAKE} FLAGS="-fsyntax-only"
debug:
	@${MAKE} FLAGS="-g"
release:
	@${MAKE} FLAGS="-DNDEBUG" OPT_FLAGS="-O2 -flto"
profile:
	@${MAKE} ${foreach target, ${ALL}, ${target}.profile}
optimize:
	@${MAKE} ${foreach target, ${ALL}, ${target}.optimized}

reset:
	reset && ${MAKE} -B

%.profile:%.cpp
	@${MAKE} FLAGS="-pg -DNDEBUG" -B
	@./${*F} >& /dev/null
	@gprof ./${*F} > $@
	@rm -f ${*F} gmon.out

%.profile:%.c
	@${MAKE} FLAGS="-pg" -B
	@./${*F} >& /dev/null
	@gprof ./${*F} > $@
	@rm -f ${*F} gmon.out

%.optimize:%.cpp
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F}
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda

%.optimized:%.c
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F}
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda
	@mv ${*F} $@

test:Index.h Grid.h Iterator.h

