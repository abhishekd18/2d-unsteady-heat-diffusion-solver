CC = mpic++
SOURCE = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCE))
EXECUTABLE = 2d_Unsteady_Diffusion
VTK_CPPFLAGS=-I/usr/include/vtk-5.8
CFLAGS =-O3 -Wno-deprecated -Wall $(VTK_CPPFLAGS)
VTK_LDFLAGS=-L/usr/lib
LDFLAGS = $(VTK_LDFLAGS)
LIBS = -lvtkCommon -lvtkFiltering -lvtkGraphics -lvtkIO -lvtkRendering -lvtkWidgets -lvtkHybrid

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o  $(EXECUTABLE) $(LIBS)
	@echo DONE!

-include $(OBJECTS:.o=.d)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $*.cpp -o $*.o
	@$(CC) -MM $(CFLAGS) $*.cpp > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -rf ../run/proc* *.o *.d *.vtk $(EXECUTABLE)
	@echo ALL CLEANED UP!

rebuild:
	make clean
	make
