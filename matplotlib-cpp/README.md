# matplotlib from c++

Download the code

```shell
git clone https://github.com/lava/matplotlib-cpp
```

Try to compile an example.

## MACOS using homebrew and miniforge

```shell
g++ basic.cpp -std=c++11 -DWITHOUT_NUMPY \
    -I /opt/homebrew/Caskroom/miniforge/base/include/python3.12 \
    -L /opt/homebrew/Caskroom/miniforge/base/lib -lpython3.12
```

or

```shell
g++ basic.cpp -std=c++11 -DWITHOUT_NUMPY \
    -I$(python3-config --includes) \
    -L$(python3-config --prefix)/lib -lpython3.12
```
