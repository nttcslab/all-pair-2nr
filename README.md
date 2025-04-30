# Efficient and Exact Algorithm for All Pair-Wise Network Reliability and Its Applications to Enhance Unreliable Pairs

This repository includes the codes to reproduce the results of experiments in the paper "Efficient and Exact Algorithm for All Pair-Wise Network Reliability and Its Applications to Enhance Unreliable Pairs."

## Build

If your environment has CMake version >=3.8, you can build all codes with the following commands:

```shell
(moving to src/ directory)
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running this, all binaries are generated in `release/` directory. Instead of `make`, you can build individual binary by the following commands:
```shell
make main          // Building the proposed method for solving 2-NR++
make maincsrel     // Building the SOTA method for solving 2-NR++
make augment       // Building the proposed method for link augmentation
make augmentcsrel  // Building the SOTA method for link augmentation
```

## Limitation

In our experiment, we include the algorithm for determining the link ordering by [the beam-search based heuristics of Inoue and Minato](https://www-alg.ist.hokudai.ac.jp/~thomas/TCSTR/TR/ABSTRACTS/abstr_tcstr16_80.html). This is due to the following reasons: (i) we included the elapsed time for computing the link ordering, and (ii) we computed the link ordering after performing some preprocessings such as removing degree 1 nodes. However, since the code for it is currently not publicly available, we excluded the codes for determining the link ordering. Instead, we applied this algorithm for the original graph (before performing preprocessings) and wrote the computed link ordering as the edgelist. This may produce similar experimental results, but not exactly the same.

## How to reproduce experimental results

All data used in our experiments are in `data.tar.gz`. After extracting, `*.txt` describes the graph files (as an edgelist), and `*.txt.prob` specifies each link's working probability.

### Section V-A

The proposed method can be executed by the following command:

```shell
./main [graph_file] [probability_file] [order_file]
```

`[graph_file]`, `[probability_file]` specify the path to the file describing the edgelist of the graph and each link's availability, respectively. When $T^\prime=\emptyset$, `[order_file]` specifies the path to the file describing the ordering among links. Note that we already computed the link order by the path-decomposition heuristics and wrote it on `*.txt`, so `[order_file]` can be specified with the same path as `[graph_file]`. We can also specify `!` for the input of `[order_file]`; it will run the link ordering determination adaptively. However, due to the limitation described above, it currently reproduces the same result as when we specify the same argument as `[graph_file]` for `[order_file]`.

After execution, the program computes $R(u,v)$ for every node pair $u,v\in V$, i.e., solves 2-NR++. The output consists of the following lines:

```shell
(u)-(v) : (the value of R(u,v))
```

_Running example:_ To solve 2-NR++ on Interoute topology, run:

```shell
./main ../data/0146-real-Interoute.edgelist.txt ../data/0146-real-Interoute.edgelist.txt.prob ../data/0186-real-UsCarrier.edgelist.txt
```

The SOTA method can also be executed by the following command:

```shell
./maincsrel [graph_file] [probability_file] [order_file]
```

### Section VI-B

Link augmentation problem can also be solved by the following command:

```shell
./augment [graph_file] [probability_file] [order_file]
```

After execution: the program computes the least reliable pair and its reliability value for every augmented link. Finally, the program outputs the best link to augment and the least reliable pair when the best link is augmented.

```shell
(The endpoints of the link to augment) : (The least reliable node pair) (The reliability value)
...
BEST:
(The endpoints of the best link) : (The least relabile node pair) (The reliability value)
```

_Running example:_ To solve link augmentation problem on RF-1221 topology, run:

```shell
./augment ../data/1221.dat.txt ../data/1221.dat.txt.prob ../data/1221.dat.txt
```

The link augmentation problem with the SOTA method can also be executed by the following command:

```shell
./augmentcsrel [graph_file] [probability_file] [order_file]
```

## License

This software is released under the NTT license, see `LICENSE.pdf`.
