CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -O0 -p -g
IFLAGS = -Iinclude

OBJ = ferrari.o Index.o IntervalList.o Graph.o Timer.o

ferrari: $(OBJ)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o ferrari $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf $(BIN) $(OBJ) ferrari