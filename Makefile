# Warnings frequently signal eventual errors:
CXXFLAGS= -std=c++11 -g -O3 #-W -Wall -Weffc++ -Wextra -pedantic -O3

#Inludes
INCLUDES=#-I/Users/emmanueljohn/libs/NetworKit/include #-I/home/emmanuj/root/include
LFLAGS= #-L/Users/emmanueljohn/libs/NetworKit

# Linker flags for both OS X and Linux
LDFLAGS=#-fopenmp #-lexpat #-ltcmalloc -lprofiler -fopenmp #-lNetworKit

# Generates list of object files from all the
#   source files in directory
OBJS = $(addsuffix .o, $(basename $(shell find *.cpp)))
DEPS = $(OBJS:.o=.d)

# Set executable name
EXEC = msparse

# Declare the phony targets
.PHONY: echo clean r clang gcc setclang setgcc vg

# Phony targets to run dependencies in order
gcc: | setgcc $(EXEC)
clang: | setclang $(EXEC)

# For use with the clang static analyzer by
#  using the environment variable for CXX
sb: | $(clean) $(EXEC)

vg: $(EXEC)
	valgrind ./$(EXEC)

# Run the executable
r:
	./$(EXEC)

clean:
	rm -rf $(OBJS)
	rm -rf $(DEPS)
	rm -rf $(EXEC)

# Phony target to use clang for compile and linking
setclang:
	@echo "Setting clang"
	$(eval CXX = clang++)
	$(eval CXX_LINK = clang++)

# Phony target to use g++ for compile and linking
setgcc:
	@echo "Setting g++-5"
	$(eval CXX = g++-5)
	$(eval CXX_LINK = g++-5)

# $< refers to the first dependency
# Uses static pattern rule to keep from compiling all
#   objects every time.
$(OBJS): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LFLAGS) -c $< -o $@ $(LDFLAGS)

$(DEPS): %.d: %.cpp
	@echo "Generating "$@
	@set -e; rm -f $@; \
      $(CXX) -MM $(CXXFLAGS) $(INCLUDES) $< > $@.$$$$; \
      sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
      rm -f $@.$$$$
#g++ -o LibDemo -std=c++11 -I/home/emmanuj/libs/networkit/include -L/home/emmanuj/libs/networkit/ LibDemo.cpp -lNetworKit -fopenmp
# $@ refers to the target
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(INCLUDES) $(LFLAGS) $(LDFLAGS)

include $(DEPS)
