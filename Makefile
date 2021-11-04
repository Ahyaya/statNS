[build]:
	gcc src/statNS.c -lm -lpthread -shared -fPIC -o lib/libstatNS.so

test:
	gcc test.c lib/libstatNS.so -lm -lpthread -o test.out
	./test.out
clean:
	rm -f *.out
