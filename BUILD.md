# Windows
Windows => Use precompiled binary, or compile it with VS2008/2010 (Express or pro, Pro will allow you to enable Opemp in CMVS)
        => Use CMake GUI in order to generate the Visual Studio project file (in ./program you will find the main CMakeLists.txt).

# Linux compilations (Ubuntu used as example)
```
#Prepare and empty machine for building:
sudo apt-get update -qq && sudo apt-get install -qq
sudo apt-get -y install git jpeg boost boost-graph
git clone https://github.com/pmoulon/CMVS-PMVS
mkdir CMVS-PMVS_build && cd CMVS-PMVS_build
cmake ../CMVS-PMVS/program
make
sudo make install
```
