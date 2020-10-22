# Kendall-Correlation-for-Large-Spike-Trains
A method for quickly calculating the Kendall correlation for large spike trains based on the method presented in the paper "An O(n) Method for Calculating Kendall Correlations of Spike Trains" William T. Redman 2018. 

This method, unlike Knight's method which is O(n ln n) (Knight 1966), is O(n). We have shown it to be around 50x faster for various neurally plausible sparseness regimes and, for spike trains in the tens of thousands of elements long, tens of thousands of times of faster than MATLAB's O(n^2) implementation. 

If you use this method in your research, please consider citing the paper.
