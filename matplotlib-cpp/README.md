# matplotlib from c++

Download the code

```shell
git clone https://github.com/lava/matplotlib-cpp
```

Try to compile an example.

## MACOS using homebrew and miniforge

```shell
g++ basic.cpp -std=c++11 -DWITHOUT_NUMPY \
    -I /usr/local/Caskroom/miniforge/base/include/python3.9 \
    -L /usr/local/Caskroom/miniforge/base/lib -lpython3.9
```
