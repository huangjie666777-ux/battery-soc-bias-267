CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic -Iinclude -isystem third_party/eigen3

BUILD_DIR := build
LIB := $(BUILD_DIR)/libbattery_soc267.a
DEMO := $(BUILD_DIR)/battery_soc267_demo

SOURCES := src/battery_model.cpp src/augmented_ekf.cpp src/soc_estimator.cpp
OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

.PHONY: all clean test

all: $(LIB) $(DEMO)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB): $(OBJECTS)
	ar rcs $@ $^

$(DEMO): test/demo.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $< $(LIB) -o $@

test: $(DEMO)
	./$(DEMO)

clean:
	rm -rf $(BUILD_DIR)
