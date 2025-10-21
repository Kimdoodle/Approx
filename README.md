# Approx

## Env
- Ubuntu 22.04 LTS

## Install
git clone https://github.com/Kimdoodle/Approx.git
cd Approx
git clone https://github.com/CPET-lab/CPET_SEAL.git
cd CPET_SEAL
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=install
cmake --build build --config Release
cmake --install build
