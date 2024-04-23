all: main

main: main.c
	gcc -o $@ $^ -lm -g

run:
	./main

clean:
	rm main