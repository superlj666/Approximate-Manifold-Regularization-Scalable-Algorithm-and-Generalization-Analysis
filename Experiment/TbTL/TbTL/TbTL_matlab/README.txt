Codes for Large-Scale Graph-based Semi-Supervised Learning via Tree Laplacian Solver
Copyright 2015 Yan-Ming Zhang (ymzhang@nlpr.ia.ac.cn)


******************************************************************
This software is currently released under the GNU General Public License 
http://www.gnu.org/copyleft/gpl.html

If you use this code in research for publication, please cite our AAAI2015
paper "Large-Scale Graph-based Semi-Supervised Learning via Tree Laplacian Solver". 

Please provide feedback if you have any questions or suggestions via 
ymzhang@nlpr.ia.ac.cn.
*******************************************************************

 
0. Installation
==================

 mex TbTL_matlab.cpp ../TbTL_exe/detect_disconnected_components.cpp ../TbTL_exe/llgc_on_tree.cpp ../TbTL_exe/mini_spanning_tree.cpp ../TbTL_exe/shortest_path_tree.cpp ../TbTL_exe/random_spanning_tree.cpp ../TbTL_exe/build_graph.cpp


1. Usage
==================


f = TbTL_matlab(edges, n, y, k, tree_type, n_tree, alpha);


f : a n-by-1 vector that contains the predicted labels for all nodes. 
edges : a m-by-3 matrix in which each row [idx1,idx2,weight] is one edge in the graph. The node index must be continuous integers and start from 1.
n : the number of nodes in the graph.
y : a n-by-1 vector in which  y_i is the label of node i for labeled nodes, otherwise y_i=-1.
k : the number of classes.
tree_type : the type of spanning trees (0 -- minimum spanning tree, 1 -- shortest path tree, 2 -- random spanning tree).
n_tree : the number of spanning trees.
alpha : the regularization factor





