# Created on Fri 27 Jul 2018 13:47:37 AEST

.PHONY: build install test clean pypi

build:
	python setup.py sdist bdist_wheel

install:
	python setup.py install

pypi:
	# not uploading any wheels
	rm dist/*.whl
	twine upload dist/*

test:
	pytest

clean:
	rm -rf dist/
