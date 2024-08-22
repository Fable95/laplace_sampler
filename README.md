# laplace_sampler

This repository contains the code for the publication: "Practical two-party computational differential privacy with active security"

Find it on [eprint](https://eprint.iacr.org/2024/4)
## Installation of dependencies

Setup repo by cloning the submodule and patching the updates.

```bash ./update_changes.sh```

Install the requirements for [MP-SPDZ](https://github.com/data61/MP-SPDZ):
```
sudo apt-get install automake build-essential clang cmake git libboost-dev libboost-thread-dev libgmp-dev libntl-dev libsodium-dev libssl-dev libtool python3
```
Run the following commands to setup the library:
```
make clean-deps boost libote
```

## Compile and run
To compile the Laplace Sampler, go to ```MP-SPDZ``` directory and run: (num threads dictates the number of threads used when executing make)
```
make sampler -j <number threads>
```
Go to ```MP-SPDZ/bin``` and find the executable ```sampler```. Execute the program in two terminals:

```./sampler <index of party> <number of parties> <security param>```

On two terminals with a $\kappa = 40$ this would look like:

```
./sampler 0 2 40
./sampler 1 2 40
```
