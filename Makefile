CC = cc
LD_FLAGS = -lm
CC_FLAGS = --std=c17

# Nombres de archivos:

EXEC = io
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.c=.o)

# Target principal:

all: $(EXEC) removeobj

$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LD_FLAGS)

# Obtener los .o
%.o: %.c
	$(CC) $(CC_FLAGS) -c $< -o $@

# Borrar los ejecutables generados:
removeobj:
	rm -f $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(EXEC)
