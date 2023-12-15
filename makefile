
#used chatgpt for this


# Define the compiler
CC = gcc

# Compiler flags
CFLAGS = -Wall -g -pg -lm

# Source files
SOURCES = GA.c functions.c OF.c

# Include directories
INCLUDES = -I.

# Libraries
LIBS = -lm

# Output executable
TARGET = GA

# Default rule
all: $(TARGET)

# Compile the source files and link them to create the executable
$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(SOURCES) $(LIBS)

# Clean up object files and the executable
clean:
	rm -f $(TARGET)

# PHONY rule to avoid conflicts with filenames
.PHONY: clean