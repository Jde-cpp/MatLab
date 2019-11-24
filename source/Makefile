CCACHE_EXISTS := $(shell ccache -V)
ifdef CCACHE_EXISTS
    CC := ccache $(CC)
    CXX := ccache g++-8
endif
DEBUG ?= 1
CXXFLAGS_BASE = -c -g -pthread -fPIC -std=c++17 -Wall -Wno-unknown-pragmas -Wno-ignored-attributes

BIN_DIR_BASE:=../../bin/
OUT_DIR_BASE := .obj/
LIBS = 
ifeq ($(DEBUG), 1)
	OUT_DIR=$(OUT_DIR_BASE)debug
	CXXFLAGS= $(CXXFLAGS_BASE) -O0 -I$(OUT_DIR) -fsanitize=address -fno-omit-frame-pointer
	BIN_DIR=$(BIN_DIR_BASE)debug
else
	OUT_DIR=$(OUT_DIR_BASE)release
	CXXFLAGS= $(CXXFLAGS_BASE) -DNDEBUG -march=native -O3 -I$(OUT_DIR)
	BIN_DIR=$(BIN_DIR_BASE)release
endif
OUTPUT_FILE = libMat.so
OUTPUT = $(BIN_DIR)/$(OUTPUT_FILE)

all: $(OUTPUT)
	$(NOECHO) $(NOOP)

OBJECTS = $(OUT_DIR)/MatLabFile.o $(OUT_DIR)/MatLabVariable.o
LDFLAGS = -pthread -shared -Wl,-z,origin -Wl,-rpath='$$ORIGIN'

$(OUTPUT): $(OUT_DIR)/stdafx.h.gch $(OBJECTS)
	$(CXX) -L$(BIN_DIR) $(LDFLAGS) $(OBJECTS) -o$(OUTPUT) $(LIBS)

$(OUT_DIR)/stdafx.h.gch: pc.h
	$(CXX) $(CXXFLAGS) pc.h -o $(OUT_DIR)/stdafx.h.gch -I/home/duffyj/code/libraries/spdlog/include -I$$BOOST_ROOT -I../../Framework/source
#-I/home/duffyj/code/libraries/json/include  
$(OUT_DIR)/%.o:%.cpp %.h pc.h
	$(CXX) $(CXXFLAGS)  ./$< -o$@ -c -I/home/duffyj/code/libraries/eigen

clean:
	rm -rf -d $(OUT_DIR)/*.*
	rm -rf -d $(OUTPUT)
