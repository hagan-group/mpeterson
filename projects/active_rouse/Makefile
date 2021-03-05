CC = g++ -std=c++17
CFLAGS = -Iinclude -Wall -Wextra -pedantic
LFLAGS = -lstdc++fs
RELEASE = -O3 -march=x86-64 -DNDEBUG
DEBUG = -g -Og 

all: test debug release

debug:
	mkdir -p build/debug
	$(CC) $(CFLAGS) $(DEBUG) source/rouse.cpp -o build/debug/rouse $(LFLAGS)
	$(CC) $(CFLAGS) $(DEBUG) source/semiflexible.cpp -o build/debug/semiflexible $(LFLAGS)
	$(CC) $(CFLAGS) $(DEBUG) source/semiflexible_small_dt.cpp -o build/debug/semiflexible_small_dt $(LFLAGS)

release:
	mkdir -p build/release
	$(CC) $(CFLAGS) $(RELEASE) source/rouse.cpp -o build/release/rouse $(LFLAGS)
	$(CC) $(CFLAGS) $(RELEASE) source/semiflexible.cpp -o build/release/semiflexible $(LFLAGS)
	$(CC) $(CFLAGS) $(RELEASE) source/semiflexible_small_dt.cpp -o build/release/semiflexible_small_dt $(LFLAGS)

test:
	mkdir -p build/test
	$(CC) $(CFLAGS) $(DEBUG) source/test.cpp -o build/test/test $(LFLAGS)

clean:
	rm -rf build/

.PHONY : all
