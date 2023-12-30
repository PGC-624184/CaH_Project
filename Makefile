install:
	#install
	pip install --upgrade pip && pip install -r requirements.txt
format:
	#format
	find src -name '*.py' -exec black {} +
lint:
	#lint
	- pylint --disable=R,C src/*.py
