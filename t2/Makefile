all: 
	python2.7 t2_datagen.py > data.txt
	make -C mat
	make -C sim
	make -C doc
	cp doc/Capa.pdf .

clean:
	make -C mat clean
	make -C sim clean
	make -C doc clean

.PHONY: all clean 
