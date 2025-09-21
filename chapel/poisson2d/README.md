# 2d Poisson equation

## Single locale version
See https://chapel-lang.org/blog/posts/bns2/

```shell
chpl nsPoisson.chpl --fast
./nsPoisson --sourceMag=500 --makePlots=true
open initial.png solution.png
```

## Distributed version

See https://chapel-lang.org/blog/posts/bns3/
