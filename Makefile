CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

SRCS = src/tsp.cpp
OBJS = $(SRCS:.cpp=.o)

TARGETS = src/christofides_edmonds src/christofides_gabow src/mst2approx src/heldkarp src/hylos

all: $(TARGETS)

src/christofides_edmonds: src/christofides_edmonds.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

src/christofides_gabow: src/christofides_gabow.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

src/mst2approx: src/mst2approx.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

src/heldkarp: src/heldkarp.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

src/hylos: src/hylos.cpp src/hylos_tsp.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGETS) 