ShapedSlp:
===============

Implementation for the experiments in the paper:
```
Practical Random Access to SLP-Compressed Texts
Travis Gagie, Tomohiro I, Giovanni Manzini, Gonzalo Navarro, Hiroshi Sakamoto, Louisa Seelbach Benkner, Yoshimasa Takabatake
Proc. 27th International Symposium on String Processing and Information Retrieval (SPIRE) 2020, pp. 221--231, 2020
https://dx.doi.org/10.1007/978-3-030-59212-7_16
```

The following external libraries are used:
- simongog/sdsl-lite
- vigna/sux


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

The following executables for benchmarks are created. See help by running without options.

```sh
./SlpEncBuild
./SubstrBenchmark
```


### Examples
`./SlpEncBuild` builds SLP encoding and outputs it to file. Input SLPs under the format of the output of NavarroRepair [https://users.dcc.uchile.cl/∼gnavarro/software/repair.tgz] or Bigrepair [https://gitlab.com/manzai/bigrepair] are supported.
```sh
# build encoding for method SelfShapedSlp_SdSd_Sd
./SlpEncBuild -i path_to_data/chr19x1000.fa -o ds_fname -e SelfShapedSlp_SdSd_Sd -f Bigrepair
# With "-e All", build encodings for all methods. The file name for each method is prefixed by the string given by -o option and suffixed by the method name
./SlpEncBuild -i path_to_data/chr19x1000.fa -o chr_ -e All -f Bigrepair
```

`./SubstrBenchmark` runs benchmarks to test substring extraction speed of SLP encodings created by `./SlpEncBuild`
```sh
# test substring extraction speed
./SubstrBenchmark -i ds_fname -e SelfShapedSlp_SdSd_Sd -n 10000 -l 10 -j 1000000
# With "-e All", test all methods. The base name should be given by -i option.
./SubstrBenchmark -i chr_ -e All -n 10000 -l 10 --dummy_flag 1 -j 1000000
```


### Methods

The following encoding methods of SLPs with random access are implemented.
- PlainSlp: An input SLP is converted into a binary SLP so that its derivation tree becomes a binary tree. Then the left/right child of variables are stored plainly in FixedBitLenCode or IncBitLenCode. Expansion lengths of variables are stored independently.
- PoSlp: The encoding method proposed in [Maruyama et al., Fully-Online Grammar Compression, SPIRE 2013] based on post order SLPs.
- SelfShapedSlp: The encoding method proposed in [Gagie et al., Practical Random Access to SLP-Compressed Texts, SPIRE 2020] that utilizes the expansion length of a variable for shaving off the information to identify the variable.
- ShapedSlp: A prototype of SelfShapedSlp.


### Components

Basic data structures used as components of encoding methods can be switched by template arguments.
- FixedBitLenCode: A vector that uses fixed-length bits for each element. By default, the bit-length is set to be the smallest one to store the maximum element in the vector.
- DirectAccessibleGammaCode: Direct accessible Gamma Code that uses a succinct select data structure to encode the bit-length part of Elias Gamma Code. Select data structures are switchable by template arguments.
- IncBitLenCode: A vector that uses lg(i + c) bits for the i-th element under the assuption that the i-th element of a given vector has a value at most i + c for a predefined constant c.
