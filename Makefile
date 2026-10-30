CXX      = g++
CXXFLAGS = -O2 -fopenmp -Wall -Wextra
LDFLAGS  = -lm

TARGET = sph
SRC    = sph.cpp

.PHONY: all run fresh animate rocks format clean
# prerequisites of `all` must run in order: simulate, then plot
.NOTPARALLEL:

# `clean` is deliberately NOT part of `all`: it deletes every .dat file, and
# losing a finished run to an automatic cleanup is not worth the tidiness.
# Run `make clean` yourself when you actually want the data gone.
all: run animate rocks

$(TARGET): $(SRC) geometry.h kernel.h rocks.h
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

run: $(TARGET)
	./$(TARGET)

# Start a run with the old frames cleared. Frame numbering restarts at 001 every
# run, so without this an aborted run leaves the previous run's tail behind and
# the frames become a mix of unrelated states.
fresh: $(TARGET)
	rm -f *.dat
	./$(TARGET)

# Render the frames in DIR to an mp4 with PyVista (real z-buffer, so occlusion is
# correct when the camera is rotated). Override on the command line, e.g.
#   make animate DIR=part_1 OUT=part1.mp4
#   make animate ARGS="--opacity 0.5 --rotate"
DIR    = .
OUT    = animation.mp4
ARGS   =
PYTHON = .venv/bin/python

animate:
	$(PYTHON) render_pyvista.py $(DIR) --out $(OUT) $(ARGS)

# One still picture of the rocks only (no fluid), with the box outline. The rocks
# come from DIR/rocks.txt, written by every run, so run the solver once after
# editing rocks.csv. Sizes follow DIR's last frame (full size after the ramp).
#   make rocks DIR=part_1 ROCKS_OUT=part1_rocks.png
#   make rocks ARGS="--azimuth 0 --elevation 0 --window-size 1280 960"
ROCKS_OUT = rocks.png

rocks:
	$(PYTHON) render_pyvista.py $(DIR) --no-fluid --out $(ROCKS_OUT) $(ARGS)

# Apply the layout in .clang-format to the C++ sources, in place.
format:
	clang-format -i $(SRC) geometry.h kernel.h rocks.h

clean:
	rm -f $(TARGET) *.dat
