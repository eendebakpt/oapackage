ifdef DEBUG
	FLAGS=-DDEBUG=1 --debug
endif

all:
	g++ -o bit_array_test $(FLAGS) -Wall bit_array_test.cpp bit_array.cpp

clean:
	if test -e bit_array_test; then rm bit_array_test; fi
	for file in $(wildcard *.dSYM); do rm -r $$file; done
	for file in $(wildcard *.greg); do rm $$file; done
