install:
	python setup.py install --prefix ~/software/python-modules
build:
	python setup.py build --prefix ~/software/python-modules
clean:
	rm -rf build
