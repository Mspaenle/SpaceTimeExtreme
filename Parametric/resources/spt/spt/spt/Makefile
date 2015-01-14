CC=gcc

INCLUDE=-I/Library/Frameworks/R.framework/Resources/include -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include

CFLAGS=-DNDEBUG -fPIC  -Wall -mtune=core2 -g -O2 -fopenmp

LDFLAGS=-dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation

LD_LIBRARY_PATH=-L/usr/local/lib -lgomp

LIB=sptlib.so
EXEC=exe

all: $(LIB)

$(LIB): timer.o spt.o main.o
	$(CC) $(LD_LIBRARY_PATH)  -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $< $(CFLAGS)

$(EXEC): timer.o spt.o main.o
	$(CC) $(LD_LIBRARY_PATH) -o $@ $^

clean:
	rm -rf *.o $(LIB) $(EXEC)
