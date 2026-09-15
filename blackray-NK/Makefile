CXX ?= g++
BUILD_DIR ?= build
PREFIX ?= /usr/local

CXXFLAGS ?= -O3 -march=native -std=c++17
CXXFLAGS += -Wall -Wextra -pedantic -fopenmp -fPIC -Iinclude
LDFLAGS ?=
LDLIBS ?= -fopenmp

CORE_SOURCES := src/metric.cpp src/diffeqs.cpp src/isco.cpp src/raytracer.cpp src/convolver.cpp src/blackray_c_api.cpp
CORE_OBJECTS := $(CORE_SOURCES:src/%.cpp=$(BUILD_DIR)/%.o)
CLI_OBJECT := $(BUILD_DIR)/blackray_cli.o
CORE_LIBRARY := $(BUILD_DIR)/libblackray_core.so
CLI := $(BUILD_DIR)/blackray_cli

.PHONY: all clean install run

all: $(CORE_LIBRARY) $(CLI)

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CORE_LIBRARY): $(CORE_OBJECTS)
	$(CXX) -shared $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(CLI_OBJECT): src/blackray_cli.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CLI): $(CLI_OBJECT) $(CORE_LIBRARY)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(CLI_OBJECT) -L$(BUILD_DIR) -lblackray_core $(LDLIBS) \
		-Wl,-rpath,'$$ORIGIN'

run: $(CLI)
	LD_LIBRARY_PATH=$(BUILD_DIR) $(CLI)

install: all
	install -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/include
	install -m 755 $(CORE_LIBRARY) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(CLI) $(DESTDIR)$(PREFIX)/bin/
	install -m 644 include/*.hpp $(DESTDIR)$(PREFIX)/include/

clean:
	rm -rf $(BUILD_DIR)
