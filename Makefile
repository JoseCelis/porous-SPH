CXX      = g++
CXXFLAGS = -O2 -fopenmp
LDFLAGS  = -lm

TARGET = sph
SRC    = sph.cpp

.PHONY: all run fresh animate clean
# prerequisites of `all` must run in order: simulate, then plot
.NOTPARALLEL:

# `clean` is deliberately NOT part of `all`: it deletes every .dat file, and
# losing a finished run to an automatic cleanup is not worth the tidiness.
# Run `make clean` yourself when you actually want the data gone.
all: run animate

$(TARGET): $(SRC) kernel.h
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

clean:
	rm -f $(TARGET) *.dat
