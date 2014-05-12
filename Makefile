CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=RNASeqFold.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=RNASeqFold

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

clean:
	rm -rf *o RNASeqFold
