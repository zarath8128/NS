override OPT_FLAGS+=
override FLAGS += -Wall -Wextra ${OPT_FLAGS}
override CFLAGS+=${FLAGS}
override CXXFLAGS+=${FLAGS} -std=c++1y
ALL=test

.PHONY:all clean syntax debug release profile optimize reset

all:${ALL} 
clean:
	rm -rf ${ALL} *.o *.out *.gcda
syntax:
	@${MAKE} FLAGS="-fsyntax-only"
debug:
	@${MAKE} FLAGS="-g"
release:
	@${MAKE} FLAGS="-DNDEBUG" OPT_FLAGS="-O2 -flto"
profile:
	@${MAKE} FLAGS="-pg"
optimize:
	@${MAKE} ${foreach target, ${ALL}, ${target}.optimize}

reset:
	reset && ${MAKE} -B

test:Index.h Grid.h

%.optimize:%.cpp
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F}
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda

%.optimize:%.c
	@${MAKE} FLAGS="-fprofile-generate" 
	@./${*F}
	@${MAKE} FLAGS="-fprofile-use" -B
	@rm -f ${*F}.gcda
