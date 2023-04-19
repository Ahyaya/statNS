CC = gcc
CLINK = -lm -lpthread
CFLAGS = -O3

DYNLIB = lib/libstatNS.so
STATICLIB = lib/libstatNS.a

HEADERS = $(wildcard include/*.h)
OBJECTS = src/statNS.o

TESTOBJ = test.o

%.o: %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@

lib: $(DYNLIB) $(STATICLIB)
	@echo "lib: native build success!"

$(DYNLIB): $(OBJECTS)
	$(CC) -shared -fPIC $(CFLAGS) $(CLINK) $(OBJECTS) -o $(DYNLIB)

$(STATICLIB): $(OBJECTS)
	ar rcs $(STATICLIB) $(OBJECTS)

test: test.out
	@echo "Test start"
	./test.out

test.out: $(DYNLIB) $(TESTOBJ)
	$(CC) $(CFLAGS) $(TESTOBJ) $(DYNLIB) -o test.out

clean:
	rm -f *.o
	rm -f src/*.o
	rm -f *.out

cleanall: clean
	rm -f lib/*.so
	rm -f lib/*.a
	rm -f EoS_lib/test*.txt