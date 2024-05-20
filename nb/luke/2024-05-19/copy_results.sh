#!/bin/bash
for date in \
2022-07-27 \
2022-07-28 \
2022-08-01 \
2022-08-03 \
2022-08-04 \
2022-08-06 \
2022-08-08 \
2022-08-11 \
2022-08-12 \
2022-08-13 \
2022-08-15 \
2022-09-12 \
2022-09-14 \
2022-09-16
do
    for region in hpc pfc
    do
        cp "/mnt/home/larend/ceph/results/20240519/3461462/date=${date}_region=${region}.pkl" \
        "/mnt/home/larend/kfx/data/${date}/${date}-${region}-csd.pkl"
    done
done
