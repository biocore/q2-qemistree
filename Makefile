.PHONY: all lint test test-cov install dev clean distclean

all: ;

lint:
	flake8 --exclude ./versioneer.py

test: all
	py.test

test-cov: all
	py.test --cov=q2_qemistree

install: all
	python setup.py install

dev: all
	pip install -e .

clean: distclean

distclean: ;
