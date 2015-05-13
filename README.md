# Discriminant Non-Negative Matrix Factorization

[![Travis-CI Build Status](https://travis-ci.org/zhilongjia/DNMF.png?branch=master)](https://travis-ci.org/zhilongjia/DNMF)

Discriminant Non-Negative Matrix Factorization is to extend the Non-negative Matrix Factorization algorithm in order to extract features that enforce not only the spatial locality, but also the separability between classes in a discriminant manner. Two kinds of Discriminant Non-Negative Matrix Factorization were implemented so far.

![Type 1 motiflogo](figure/DNMF.png)

Reference: 
+ Zafeiriou, Stefanos, et al. [*Exploiting discriminant information in nonnegative matrix factorization with application to frontal face verification.*](http://www.ncbi.nlm.nih.gov/pubmed/16722172) Neural Networks, IEEE Transactions on 17.3 (2006): 683-695.
+ Kim, Bo-Kyeong, and Soo-Young Lee. [*Spectral Feature Extraction Using dNMF for Emotion Recognition in Vowel Sounds.*](http://link.springer.com/chapter/10.1007%2F978-3-642-42051-1_59) Neural Information Processing. Springer Berlin Heidelberg, 2013.
+ Lee, Soo-Young, Hyun-Ah Song, and Shun-ichi Amari. [*A new discriminant NMF algorithm and its application to the extraction of subtle emotional differences in speech.*](http://link.springer.com/article/10.1007%2Fs11571-012-9213-1#page-1) Cognitive neurodynamics 6.6 (2012): 525-535.

This package has been submitted to [CRAN](http://cran.r-project.org/web/packages/DNMF).

Installation:

	devtools::install_github("zhilongjia/DNMF")
