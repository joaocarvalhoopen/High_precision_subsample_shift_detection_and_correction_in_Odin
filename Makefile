all:
	odin build . -out:shift_1d_detection.exe -o:speed -no-bounds-check -microarch:native

clean:
	rm -f ./shift_1d_detection.exe

run:
	./shift_1d_detection.exe
