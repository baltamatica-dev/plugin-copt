# SPDX-License-Identifier: MIT

BEX_SDK_DIR:=/s/bex-dev/bex/baltamatica/bex
COPT_DIR:=../copt50

ALL_SRC := coptbex.cpp coptinit.c copt_computeiis.c copt_feasrelax.c copt_read.c copt_solve.c copt_write.c coptmex.c

CXXFLAGS:=-Wall -Wno-unused-function -std=c++17 -shared -O3 -fPIC -I. \
	-I$(BEX_SDK_DIR)/include/  \
	-L$(BEX_SDK_DIR)/lib -lbex  \
	-I$(COPT_DIR)/include  \
	-L$(COPT_DIR)/lib -lcopt

BEX_EXT:=dll


all: copt/main.$(BEX_EXT)

copt/main.$(BEX_EXT): $(ALL_SRC)
	$(CXX) $(CXXFLAGS) -o copt/main.$(BEX_EXT)  $(ALL_SRC) $(COPT_DIR)/lib/copt.lib

clean:
	rm -f *.o *.obj *.a *.$(BEX_EXT)
	rm -f copt/main.$(BEX_EXT)
