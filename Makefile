

GGG = g++ -O3 -Wall -march=native -fomit-frame-pointer -fexpensive-optimizations
# GGG = clang++ -O3 -Wall -march=native -fomit-frame-pointer

OBJ = cipolla_primality_main.o \
      cipolla_primality.o \
      cipolla_primality_alloc.o \
      cipolla_primality_precompute.o \
      expression_parser.a

cipolla: $(OBJ)
	$(GGG) -static -o cipolla $(OBJ) -lgmp -lpthread -lm

cipolla_primality_main.o: cipolla_primality_main.cpp cipolla_primality.h cipolla_primality_alloc.h bison.gmp_expr.tab.h
	$(GGG) -c -o cipolla_primality_main.o cipolla_primality_main.cpp

cipolla_primality_alloc.o: cipolla_primality_alloc.cpp cipolla_primality_alloc.h
	$(GGG) -c -o cipolla_primality_alloc.o cipolla_primality_alloc.cpp

cipolla_primality.o: cipolla_primality.cpp cipolla_primality.h cipolla_primality_precompute.h
	$(GGG) -c -o cipolla_primality.o cipolla_primality.cpp

cipolla_primality_precompute.o: cipolla_primality_precompute.cpp cipolla_primality_precompute.h
	$(GGG) -c -o cipolla_primality_precompute.o cipolla_primality_precompute.cpp

expression_parser.a : bison.gmp_expr.o lex.gmp_expr.o bison.gmp_expr.tab.h
	ar vr expression_parser.a bison.gmp_expr.o lex.gmp_expr.o

bison.gmp_expr.o : bison.gmp_expr.tab.c bison.gmp_expr.h
	$(GGG) -c -o bison.gmp_expr.o bison.gmp_expr.tab.c

bison.gmp_expr.tab.c bison.gmp_expr.tab.h : parser.y
	bison -d parser.y

lex.gmp_expr.o : lex.gmp_expr.c
	$(GGG) -Wno-unused-function -DYY_BUF_SIZE=65540 -DYYLMAX=65540 -c -o lex.gmp_expr.o lex.gmp_expr.c

lex.gmp_expr.c : parser.l bison.gmp_expr.tab.h
	flex parser.l

check: cipolla
	./cipolla -st

clean:
	rm -f ./cipolla $(OBJ) bison.gmp_expr.o bison.gmp_expr.tab.c bison.gmp_expr.tab.h lex.gmp_expr.o lex.gmp_expr.c


