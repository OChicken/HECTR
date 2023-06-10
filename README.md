# HECTR

(Leveled) Homomorphic Encrypted (Model Predictive) ConTRoller

## How to use

```sh
# step 1: get HECTR and build
git clone https://github.com/OChicken/HECTR.git
cd HECTR
git submodule init
git submodule update
mkdir -p lib
make

# step 2: run tests
cd tests
LD_LIBRARY_PATH=$PWD/../lib:$LD_LIBRARY_PATH
make test-hectr
./test-hectr cstr-mpc
./test-hectr cstr-hempc
```
