# IterFS
This repository contains a Matlab implementation of the algorithm described in Ordozgoiti et al, ICDM '16. 

The algorithm provides an approximate solution to the Column Subset Selection Problem.

# Usage
```
   iter_fs(data, k, MAX_ITER, ZERO)
       data: The input data matrix
       k: The number of columns to choose
       MAX_ITER: (Optional) The maximum number of iterations
       ZERO: (Optional) The threshold below which scores are set to zero, to avoid numerical issues
   
     Output:
        A set of column indices indicating the ones that were chosen

```
Check out `example.m` for an example.

For better performance, the use of an optimized implementation of the SVD is recommended. https://es.mathworks.com/matlabcentral/fileexchange/47132-fast-svd-and-pca

Author: Bruno Ordozgoiti

If you use this code, please cite the following paper:

Bruno Ordozgoiti, Sandra Gomez Canaval, Alberto Mozo, "A Fast Iterative 
Algorithm for Improved Unsupervised Feature Selection". IEEE International 
Conference on Data Mining, 2016.
