CXXFLAGS = -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format -fsanitize=address
OPT = -O3
all: a.out

solve.o: solve.cpp solve.h
	g++ -c $(CXXFLAGS) $(OPT) solve.cpp

matrix.o: matrix.cpp matrix.h
	g++ -c $(CXXFLAGS) $(OPT) matrix.cpp

main.o: main.cpp
	g++ -c $(CXXFLAGS) $(OPT) main.cpp

my_errors.o: my_errors.cpp my_errors.h
	g++ -c $(CXXFLAGS) $(OPT) my_errors.cpp

a.out: main.o matrix.o solve.o my_errors.o
	g++ $^ $(CXXFLAGS) $(OPT) -o $@

clean:
	rm -f *.o a.out