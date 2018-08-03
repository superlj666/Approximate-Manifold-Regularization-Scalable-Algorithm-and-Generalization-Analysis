#include "mex.h"
#include "..\\TbTL_exe\\llgc_on_tree.h"

/*f = llgc_on_tree_matlab(graph, n_nodes, y, n_class, tree_type, n_tree, lambda)
@graph: a n_nodes-by-3 matrix, each row is an edge
@n_nodes: 
@y: y[i] \in {0,1,....,n_classes-1} if i is labeled; y[i] = -1 if i is unlabeled
@n_classes:
@tree_type: 0, MST; 1: SPT; 2: RST
@n_tree:
@lambda: 
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    int n_nodes, n_edges, n_classes, n_trees, tree_type;    
    double *ptr, alpha;
    
    ptr = mxGetPr(prhs[1]);
    n_nodes = (int)(ptr[0]);
    
    ptr = mxGetPr(prhs[3]);
    n_classes = (int)(ptr[0]);    
    
    ptr = mxGetPr(prhs[4]);
    tree_type = (int)(ptr[0]);    

	ptr = mxGetPr(prhs[5]);
    n_trees = (int)(ptr[0]);    

	ptr = mxGetPr(prhs[6]);
    alpha = (int)(ptr[0]);    

    ptr = mxGetPr(prhs[0]);
    n_edges = mxGetM(prhs[0]);
    struct Edge *graph = (struct Edge *)malloc(n_edges*sizeof(struct Edge));
    for(int i=0; i<n_edges; i++ ){
        graph[i].x = (int)(ptr[i])-1;
        graph[i].y = (int)(ptr[n_edges+i])-1;
        graph[i].w = ptr[2*n_edges+i];
    }
    ptr = mxGetPr(prhs[2]);
    int *y = (int *)malloc(n_nodes*sizeof(int));
    for(int i=0; i<n_nodes; i++) y[i] = ptr[i];
    /* 
    mexPrintf("# of nodes: %d, # of edges: %d\n", n_nodes, n_edges);
    mexPrintf("# of classses: %d\n", n_classes);
    if( tree_type==0 )
        mexPrintf("tree type: minimum spannining tree\n");
    else if( tree_type==1 )
        mexPrintf("tree type: shortest path tree\n");
    else if( tree_type==2 )
        mexPrintf("tree type: random spanning tree\n");                
    mexPrintf("# of tree: %d\n", n_trees);    */
    	              
    llgc_on_tree(graph, n_nodes, n_edges, y, n_classes, tree_type, n_trees, alpha);
    
    plhs[0] = mxCreateDoubleMatrix(n_nodes, 1, mxREAL);
    ptr = mxGetPr(plhs[0]);
	for(int i=0; i<n_nodes; i++) ptr[i]=y[i];

    free(graph);
    free(y);        
    return;
}
