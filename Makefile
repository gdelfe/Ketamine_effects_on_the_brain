env:
	python3 -m venv env
	source env/bin/activate && \
	pip install --upgrade pip && \
	pip install -r requirements.txt

notebook:
	source env/bin/activate && jupyter notebook

mount:
	sshfs -o allow_other lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix data

unmount:
	diskutil unmount force data

.PHONY: env notebook mount unmount
