FLAGS =-Wall -Werror -Wextra -std=c++17 -pedantic
TEST_FLAGS= -lgtest -lstdc++ --coverage
SOURCES=s21_*.cc
OBJECTS=s21_*.o

ifeq ($(shell uname -s), Linux)
	TEST_FLAGS= -lgtest -lpthread -fprofile-arcs -ftest-coverage
endif

all: clean test

test: s21_matrix_oop.a
	g++ $(TEST_FLAGS) s21_matrix_oop.a test.cc -o test -lgtest
	./test

s21_matrix_oop.a:
	g++ -c $(FLAGS) $(SOURCES)
	ar rc s21_matrix_oop.a $(OBJECTS)
	ranlib s21_matrix_oop.a	

leaks: test
	leaks --atExit -- ./test
	
style: 
	clang-format --style=google -i *.cc *.h
	clang-format --style=google -n *.cc *.h

clean:
	rm -rf *.o gcov_report *.a *.gcda *.gcno test.info report/ test