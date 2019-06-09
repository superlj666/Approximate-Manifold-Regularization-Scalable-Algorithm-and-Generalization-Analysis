# Approximate Manifold Regularization: Scalable Algorithm and Generalization Analysis
Experiment part in Approximate Manifold Regularization: Scalable Algorithm and Generalization Analysis.
The paper has been accepted by IJCAI-19.

# Usage
The codes are implemented in MATLAB.
## Structure
- ./datasets: All datasets are available in https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/.
- ./data: Store processed data including kernel matrix and graph Laplacian.
- ./result: Store final results used in the paper.
- ./core_functions: Implementation of compared algorithms.
- ./parameter_tune: Tune parameters.
- ./utils: Some utils including constructing kernel matrix and graph Laplacian, drawing curves and optimal parameters setting.
## Steps
1. Download data sets into ./datasets
2. Run Exp1_*.m for experiment 1.
3. Run Exp2_test_labeled_curve.m for experiment 2.
4. Run Exp3_test_sample_curve.m for experiment.