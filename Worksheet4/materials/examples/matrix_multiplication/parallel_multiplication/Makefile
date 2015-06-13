# Include files
SOURCES=algo_mul.c master.c main.c read_file.c slave.c write_file.c 

# Compiler
# --------
CC=mpicc

CFLAGS= -Wall

# Linker flags
# ------------
LDFLAGS= 

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=matrix

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ 

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
