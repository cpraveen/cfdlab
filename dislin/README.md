# DISLIN

Set where to install DISLIN

```shell
export DISLIN=$HOME/Applications/disin
```

Download DISLIN from 

https://www.dislin.de/distributions.html

and extract it. Run the installer

```shell
sh ./INSTALL
```

We need openmotif which can be installed with brew (on MACOS)

```shell
brew install openmotif
```

Compile your code

```shell
g++ foo.cpp -I$DISLIN -L$DISLIN -L/usr/local/lib -ldiscpp -lXm
```

Add path to DISLIN library, e.g., on MACOS

```shell
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:$DISLIN
```

Run the code

```shell
./a.out
```
