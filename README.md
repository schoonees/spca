# An R Package for Sparse Principal Components Analysis

This package implements the penalized matrix decomposition of Witten, Hastie and Tibshirani (2009), and 
applies it to sparse principal components analysis. Only the special case where lasso-type penalties are used is
currently implemented. A [vignette](vignettes/) is available with some examples.

The package can be installed using:

```
devtools::install_github("schoonees/spca")
```

Note that the package was created for educational purposes and has not been tested extensively.

## References

Witten, D. M., Tibshirani, R., & Hastie, T. (2009). A penalized matrix decomposition, with applications to 
sparse principal components and canonical correlation analysis. *Biostatistics*, 10(3), 515-534.