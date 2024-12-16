CC = gcc
CCFLAGS = -Wall -O3
TARGET = pso
SRCS = main.c PSO.c OF.c

all: $(TARGET)

# Compile and link source files to create the executable
$(TARGET): $(SRCS)
	$(CC) $(CCFLAGS) -o $(TARGET) $(SRCS) -lm

clean:
	rm -f $(TARGET)
