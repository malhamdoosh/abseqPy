# Created on Fri 27 Jul 2018 13:47:37 AEST

.PHONY:
build:
	python setup.py sdist bdist_wheel

.PHONY:
install:
	python setup.py install

.PHONY:
pypi:
	twine upload dist/*

.PHONY:
test:
	pytest

.PHONY:
clean:
	rm -rf dist/
