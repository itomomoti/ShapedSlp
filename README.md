ShapedSlp:
===============
Author: Tomohiro I

### Download

To download all the necessary source codes:
```sh
git clone --recursive https://github.com/itomomoti/ShapedSlp.git
```

### Compile

```sh
mkdir build
cd build
cmake ..
make
```

### Usage

Executables for benchmarks. See help by running without options.

```sh
./SlpEncBuild
./SubstrBenchmark
```

### Example

```sh
# build encoding for method SelfShapedSlp_SdSd_Sd
./SlpEncBuild -i path_to_data/chr19x1000.fa -o ds_fname -e SelfShapedSlp_SdSd_Sd -f Bigrepair
# With "-e All", build encodings for all methods. The file name for each method is prefixed by the string given by -o option and suffixed by the method name
./SlpEncBuild -i path_to_data/chr19x1000.fa -o chr_ -e All -f Bigrepair
```

```sh
# substring
./SubstrBenchmark -i ds_fname -e SelfShapedSlp_SdSd_Sd -n 10000 -l 10 -j 1000000
# With "-e All", test all methods. The base name should be given by -i option.
./SubstrBenchmark -i chr_ -e All -n 10000 -l 10 --dummy_flag 1 -j 1000000
```
