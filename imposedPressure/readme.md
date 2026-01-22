# Instructions to Build the case

```bash
mkdir build && cd build
```

```bash
cmake .. -DCMAKE_PREFIX_PATH=`spack location -i tfel`/share/tfel/cmake -DCMAKE_INSTALL_PREFIX=../install
```

```bash
make -j 2 install
```