gcc -o3 sim-poly-blend.c -lm -o sim-poly-blend.ce -g

./sim-poly-blend.ce > tmp0 &
./sim-poly-blend.ce --outputall 1 > tmp1 &
./sim-poly-blend.ce --directed 0 --outputall 1 > tmp2 &
./sim-poly-blend.ce --directed -1 --outputall 1 > tmp3 &
./sim-poly-blend.ce --npar 100 > tmp4 &
./sim-poly-blend.ce --mut 0.001 > tmp5 &
./sim-poly-blend.ce --mut 0.01 > tmp6 &
./sim-poly-blend.ce --mut 1 > tmp7 &
./sim-poly-blend.ce --targetsize 8 > tmp8 &
./sim-poly-blend.ce --targetsize 15 > tmp9 &
./sim-poly-blend.ce --confound 1 > tmp10 &
./sim-poly-blend.ce --confound 1 --outputall 1 > tmp11 &
./sim-poly-blend.ce --confound 1 --directed 0 --outputall 1 > tmp12 &
./sim-poly-blend.ce --confound 2 > tmp13 &
./sim-poly-blend.ce --confound 2 --outputall 1 > tmp14 &
./sim-poly-blend.ce --confound 2 --directed 0 --outputall 1 > tmp15 &
./sim-poly-blend.ce --confound 3 > tmp16 &
./sim-poly-blend.ce --confound 3 --outputall 1 > tmp17 &
./sim-poly-blend.ce --confound 3 --directed 0 --outputall 1 > tmp18 &

./sim-poly-blend.ce --numr 0 --nsamp 10000000 > tmp19 &
./sim-poly-blend.ce --numr 0 --nsamp 100000000 > tmp20 &
./sim-poly-blend.ce --numr 0 --nsamp 1000000000 > tmp21 &
./sim-poly-blend.ce --targetsize 8 --numr 0 > tmp22 &
./sim-poly-blend.ce --targetsize 15 --numr 0 > tmp23 &
./sim-poly-blend.ce --numr 0 --ntile 2 --ncol 8 --nbitcol 3 > tmp24 &
./sim-poly-blend.ce --numr 0 --ntile 4 --ncol 16 --nbitcol 4 > tmp25 &
./sim-poly-blend.ce --numr 0 --ntile 8 --ncol 32 --nbitcol 5 > tmp26 &

./sim-poly-blend.ce --targetsize 32 --numr 10 --outputall 1 > tmp27 &
./sim-poly-blend.ce --targetsize 35 --numr 10 --outputall 1 > tmp28 &
./sim-poly-blend.ce --targetsize 40 --numr 10 --outputall 1 > tmp29 &

Rscript plot-structs.R out-blend-1-0.100-10-16-16-64-0-0-1-5000-1e+08
Rscript plot-structs.R out-blend-1-0.100-10-8-16-64-0-0-1-5000-1e+08
Rscript plot-structs.R out-blend-1-0.100-10-15-16-64-0-0-1-5000-1e+08
