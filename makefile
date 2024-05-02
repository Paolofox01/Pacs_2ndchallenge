CXX      ?= g++
CXXFLAGS ?= -std=c++20

#change this path
PACS_ROOT?=../pacs-examples/Examples

export DOXYFILE=./Doxyfile

#change these paths
CPPFLAGS ?= -O3 -Wall -Wno-conversion-null -Wno-deprecated-declarations -I. -I$(PACS_ROOT)/src/Utilities/ -I$(PACS_ROOT)/include/


EXEC     = main

#change this path
LDFLAGS ?= -L$(PACS_ROOT)/src/Utilities/
LIBS  ?= -lpacs

# Get all files *.cpp
SRCS=$(wildcard *.cpp)
# Get the corresponding object file
OBJS = $(SRCS:.cpp=.o)
# Get all headers in the working directory
HEADERS=$(wildcard *.hpp)
#
exe_sources=$(filter main%.cpp,$(SRCS))
EXEC=$(exe_sources:.cpp=)

#========================== ORA LA DEFINIZIONE DEGLI OBIETTIVI
.PHONY: all clean distclean doc cleandoc

cleandoc:
	-\rm -rf ./doc

.DEFAULT_GOAL = all

all: $(EXEC)

clean:
	-\rm -f $(EXEC) $(OBJS)

distclean: clean cleandoc
	-\rm -f ./doc $(DEPEND)
	-\rm -f *.out *.bak *~

doc:
	doxygen $(DOXYFILE)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

# Compile each source file into an object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
