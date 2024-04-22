env:
	python3 -m venv env
	source env/bin/activate && \
	pip install --upgrade pip && \
	pip install -r requirements.txt
	rm -rf kfx.egg-info

notebook:
	source env/bin/activate && jupyter notebook

mount:
	sshfs -o allow_other lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix data/mnt

unmount:
	diskutil unmount force data/mnt

.PHONY: env notebook mount unmount
