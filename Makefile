CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude -Ithird_party/eigen3

BUILD := build
LIB := $(BUILD)/libbattery_soc267.a
SRCS := src/model.cpp src/ekf.cpp src/estimator.cpp
OBJS := $(SRCS:%.cpp=$(BUILD)/%.o)

.PHONY: all test demo run clean

all: $(LIB) $(BUILD)/selftest $(BUILD)/demo

$(BUILD)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

$(BUILD)/selftest: tests/selftest.cpp $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -lbattery_soc267 -o $@

$(BUILD)/demo: examples/demo.cpp $(LIB)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -L$(BUILD) -lbattery_soc267 -o $@

test: $(BUILD)/selftest
	./$(BUILD)/selftest

demo run: $(BUILD)/demo
	./$(BUILD)/demo

clean:
	rm -rf $(BUILD)

-include $(OBJS:.o=.d)

