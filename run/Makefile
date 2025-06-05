CC = gcc
CFLAGS = -I/../src/
SRC = $(wildcard ../src/*.c)
OBJ = $(SRC:.c=.o)

# The simulation file can be changed by overwriting this variable, e.g.
# make SIM=kelvin_helmholtz_instability.c
SIM = viscosity_investigation.c 

TARGET = simulate
LIBS = -lfftw3f -lm

all: $(TARGET) 

# For creating simulation target
$(TARGET): $(OBJ) $(SIM)
	$(CC) $(CFLAGS) $(SIM) -o $@ $(OBJ) $(LIBS)

# For compiling object files in /src
../src/%.o: ../src/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

# A phony so that we always recompile if we change SIM
.PHONY: force
$(TARGET): force

# Cleans target in /run and .o files in /src
clean:
	rm -f $(OBJ) $(TARGET) 

# Cleans output directory of gifs, binary files and pngs
clean_outputs:
	rm -f outputs/*.gif outputs/*.bin outputs/*.png