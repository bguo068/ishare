#! /usr/bin/env
set -e
cd ..
# cargo run --release --bin ibdutils -- \
#     plot-ibd -g testdata/genome.toml \
#     -i testdata/ibd_hmmibd  \
#     -f hmmibd \
#     -s testdata/samples.hap.txt \
#     -I testdata/ibd_hmmibd2 \
#     -F hmmibd  \
#     -S testdata/samples.hap.txt \
#     -o plot -x 0 -X 2

# for i in {1..100}; do
# cargo run --release --bin ibdutils -- \
#     plot-ibd -g testdata/genome.toml \
#     -i testdata/ibd_hapibd  \
#     -f hapibd \
#     -s testdata/samples.d2h.txt \
#     -I testdata/ibd_hmmibd2 \
#     -F hmmibd  \
#     -S testdata/samples.hap.txt \
#     -o plot -x 0 -X $i
# done


# for i in {1..100}; do
# cargo run --release --bin ibdutils -- \
#     plot-ibd -g testdata/genome.toml \
#     -i testdata/ibd_hmmibd2 \
#     -f hmmibd  \
#     -s testdata/samples.hap.txt \
#     -I testdata/ibd_tskibd  \
#     -F tskibd \
#     -S testdata/samples.hap.txt \
#     -o plot -x 0 -X $i
# done

for i in {1..100}; do
cargo run --release --bin ibdutils -- \
    plot-ibd -g testdata/genome.toml \
    -i testdata/ibd_hapibd  \
    -f hapibd \
    -s testdata/samples.dip.txt \
    -I testdata/ibd_tskibd  \
    -F tskibd \
    -S testdata/samples.h2d.txt \
    -o plot -x 0 -X $i
done

cd -