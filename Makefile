CC = gcc
CLINK = -lm -lpthread
CFLAGS = -O3 -fPIC

DYNLIB = lib/libstatNS.so
STATICLIB = lib/libstatNS.a

HEADERS = $(wildcard include/*.h)
OBJECTS = src/statNS.o src/nvdata.o src/ptdata.o

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
	@echo "Simple Test"
	./test.out

testtov: test_tovec.out
	@echo "Testing for adaptive Runge-Kutta method"
	./test_tovec.out

test.out: $(DYNLIB) $(TESTOBJ)
	$(CC) $(CFLAGS) $(TESTOBJ) $(DYNLIB) -o test.out

test_tovec.out: $(DYNLIB) test_tovec.c
	$(CC) $(CFLAGS) test_tovec.c $(DYNLIB) -lm -o test_tovec.out

clean:
	rm -f *.o
	rm -f src/*.o
	rm -f *.out
	rm -f *.log

cleanall: clean
	rm -f lib/*.so
	rm -f lib/*.a
	rm -f EoS_lib/test*.txt
