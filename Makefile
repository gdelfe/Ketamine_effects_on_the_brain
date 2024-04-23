env:
	python3 -m venv env
	source env/bin/activate && \
	pip install --upgrade pip && \
	pip install -r requirements.txt
	rm -rf kfx.egg-info

notebook:
	source env/bin/activate && jupyter notebook

mount:
    mkdir -p data/mnt
	sshfs lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine

unmount:
	fusermount3 -u data/mnt

.PHONY: env notebook mount unmount
