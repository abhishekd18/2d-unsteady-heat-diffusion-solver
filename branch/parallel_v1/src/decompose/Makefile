IDIR =../src
CC=g++
CFLAGS=-I$(IDIR) -lm -O3
BIN =../../bin
ODIR=./obj
LDIR =./
SRCDIR=../src

_DEPS = decomposePar.h tri.h settings.h
DEPS = $(patsubst %,%,$(_DEPS))

_OBJ = decomposePar.o tri.o settings.o decompose.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(BIN)/decomposePar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
	@echo DONE!

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
