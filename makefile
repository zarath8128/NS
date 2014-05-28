override OPT_FLAGS+=-march=native
override FLAGS += -Wall -Wextra ${OPT_FLAGS}
override CFLAGS+=${FLAGS}
override CXXFLAGS+=${FLAGS} -std=c++1y

LDLIBS=-lglsc -lX11
ALL=test

.PHONY:all clean syntax debug release profile optimize full reset

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
full:profile optimize ${ALL}

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

%.optimized:%.cpp
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F} >& /dev/null
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda
	@mv ${*F} $@

%.optimized:%.c
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F} >& /dev/null
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda
	@mv ${*F} $@

test:Grid.h Iterator.h Staggered.h Diffusion.h Coordinate.h Monitor.o Divergence.h Gradient.h
Monitor.o:Monitor.h
