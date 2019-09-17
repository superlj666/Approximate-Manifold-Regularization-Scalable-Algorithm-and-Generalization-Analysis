# Approximate Manifold Regularization: Scalable Algorithm and Generalization Analysis
This repository provides the code used to run the experiments of the paper "Approximate Manifold Regularization: Scalable Algorithm and Generalization Analysis".
The paper has been published in [IJCAI-19](https://www.ijcai.org/proceedings/2019/0400.pdf).
The paper applied Nystom and PCG to LapRLS, borrowing the idea from [FALKON](https://arxiv.org/abs/1705.10958).
The implementation also use tricks provided in the repository (https://github.com/LCSL/FALKON_paper).


# Usage
The codes are implemented in MATLAB.
## Structure
- ./datasets: All datasets are available in [Libsvm Data](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/).
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