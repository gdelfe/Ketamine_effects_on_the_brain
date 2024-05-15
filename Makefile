env:
	python3 -m venv env
	source env/bin/activate && \
	pip install --upgrade pip && \
	pip install -r requirements.txt && \
	pip install -e . && rm -rf src.egg-info

notebook:
	source env/bin/activate && jupyter notebook

mount:
	mkdir -p data/mnt
	sshfs lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine data/mnt

unmount:
	fusermount3 -u data/mnt

.PHONY: env notebook mount unmount
