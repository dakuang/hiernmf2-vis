# Hierarchical rank-2 NMF with visualization

**Hierarchical rank-2 nonnegative matrix factorization (HierNMF2)** is an unsupervised algorithm for large-scale document clustering and topic modeling. It is about 20 times faster than LDA with comparable quality. HierNMF2 has also been successfully applied in the area of bioinformatics.

This package provides code for the algorithm used in this paper:
```
Da Kuang, P. Jeffrey Brantingham, Andrea L. Bertozzi,
Crime Topic Modeling,
Crime Science, 6(12), 2017.
```
and is based on the following paper:
```
Da Kuang, Haesun Park,
Fast rank-2 nonnegative matrix factorization for hierarchical document clustering,
The 19th ACM SIGKDD International Conference on Knowledge, Discovery, and Data Mining (KDD '13), pp. 739-747, 2013.
```
Please cite the above papers if you find the code useful.

## Prerequisites
The code was developed and tested on Ubuntu Linux 14.04 LTS. 
* Matlab R2011a or above
* Python 2.7.x
* Latex with TikZ and tikz-qtree packages

## Testing
The entry point of the code is `workflow.py`. It starts from a text file to the generation of the tree structure visualization. To test the code:
1. Create a text file `test.txt`, where each line contains a document.
2. Run `python workflow.py test.txt 20` to generate 20 leaf clusters/topics. Adjust the number of clusters/topics according to your needs.

## Example
Below is a visualization of 20 leaf topics generated on ~30,000 image captions from the [Yelp Challenge dataset](https://www.yelp.com/dataset/challenge).

![](https://dakuang.github.io/images/hiernmf2-vis-example.png)

Reference: https://github.com/dakuang/hiernmf2
