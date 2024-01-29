target=minixy_e

$(target): main.cpp -lm
		g++  $^ -O2 -o $@

.PHONY: clean
clean:
	-rm -f *.o ${target}