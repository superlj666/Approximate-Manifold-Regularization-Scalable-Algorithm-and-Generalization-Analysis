#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "llgc_on_tree.h"

//sort the edges
void compute_edge_order_by_preorder_traversal(struct Edge *tree, int n_nodes, int total, struct Edge *edge_order)
{
	struct Edge_ *pool = (struct Edge_ *)malloc(2*(n_nodes-1)*sizeof(struct Edge_));
	struct Edge_ **g = (struct Edge_ **)malloc(total*sizeof(struct Edge_ *));
	build_undirected_graph(tree, total, n_nodes-1, g, pool);

	int cnt=0;
	int stack_sz=0;
	int *stack, *parent_stack;
	stack = (int *)malloc(2*n_nodes*sizeof(int));
	parent_stack = stack + n_nodes;

	stack[stack_sz] = tree[0].x;//root
	parent_stack[stack_sz++] = -1;

	while(stack_sz){
		int cur = stack[stack_sz-1];
		if(g[cur]==NULL){
			stack_sz--;
			continue;
		}

		int chi = g[cur]->idx;
		if(stack_sz!=1 && chi==parent_stack[stack_sz-1]){
			g[cur] = g[cur]->next;
			continue;
		}
		
		edge_order[cnt].x = chi;
		edge_order[cnt].y = cur;
		edge_order[cnt++].w = g[cur]->w;

		stack[stack_sz] = chi;
		parent_stack[stack_sz++] = cur;
		g[cur] = g[cur]->next;
	}
	
	free(stack);
	free(pool);
	free(g);
}



//min_f f'Lf + (f-y)'C(f-y)
//<=>  (L+C)f = Cy
void tree_laplacian_solver(struct Edge *tree, int *nodes, int n_nodes, int total, REAL *c, int n_col, REAL *y)
{

	int i,j,k,idx;

	//1. re-order the nodes
	struct Edge *edge_order = (struct Edge *)malloc((n_nodes-1)*sizeof(struct Edge));
	compute_edge_order_by_preorder_traversal(tree, n_nodes, total, edge_order);

	//2. diagnolize the tree: (C+L) = S'*D*S
	REAL *diag = (REAL *)malloc(total*sizeof(REAL));
	
	for( i=0; i<n_nodes; i++ ) diag[nodes[i]] = c[nodes[i]];
	for( i=0; i<n_nodes-1; i++ ){//compute degrees
		diag[tree[i].x] += tree[i].w;
		diag[tree[i].y] += tree[i].w;
	}

	for( k = n_nodes-2; k >= 0; k-- ){
		i = edge_order[k].x;
		j = edge_order[k].y;

		REAL v = -edge_order[k].w/diag[i];
		diag[j] = diag[j] + v*edge_order[k].w;
		edge_order[k].w = v;
	}

	struct Edge_ *pool = (struct Edge_ *)malloc(2*(n_nodes-1)*sizeof(struct Edge_));
	struct Edge_ **S = (struct Edge_ **)malloc(total*sizeof(struct Edge_ *));
	build_directed_graph(edge_order, total, n_nodes-1, S, pool);
	struct Edge_ **S_T = (struct Edge_ **)malloc(total*sizeof(struct Edge_ *));
	reverse_directed_graph(edge_order, total, n_nodes-1, S_T, pool+n_nodes-1);

	int *node_order = (int *)malloc(n_nodes*sizeof(int));
	node_order[0] = edge_order[0].y;//root
	for( i=1; i<n_nodes; i++ ) node_order[i] = edge_order[i-1].x;
	free(edge_order);

	//3. solve S'*D*Sx = b
	REAL *b = (REAL *)malloc(total*sizeof(REAL));
	for( k = 0; k < n_col; k++ ){
		REAL *rlts = y + k*total;
		for( i = 0; i < n_nodes; i++ ){
			idx = nodes[i];	
			b[idx] = rlts[idx]*c[idx];
		}	
		for( i = n_nodes-1; i >= 0  ; i-- ){
			int idx = node_order[i];
			REAL r = b[idx];
			for( struct Edge_ *ptr = S_T[idx]; ptr; ptr = ptr->next){
				r -= rlts[ptr->idx]*ptr->w;
			}
			rlts[idx] = r;
		}
		for( i = 0; i < n_nodes; i++ ){
			idx = nodes[i];
			rlts[idx] = rlts[idx]/diag[idx];
		}
		for( i = 0; i < n_nodes; i++ ){
			int idx = node_order[i];
			REAL r = rlts[idx];
			for( struct Edge_ *ptr = S[idx]; ptr; ptr = ptr->next){
				r -= rlts[ptr->idx]*ptr->w;
			}
			rlts[idx] = r;
		}
	}

	free(pool);
	free(S);
	free(S_T);
	free(node_order);
	free(diag);/**/
}

void _llgc_on_tree(struct Edge *tree, int *nodes, int n_nodes, int total, REAL *c, int n_class, int *y)
{
	int i, k;
	//build b by the one-vs-all strategy. b is a total*n_class vector.
	REAL *b = (REAL *)malloc(n_class*total*sizeof(REAL));
	for( int k = 0; k < n_class; k++ ){
		for( i = 0; i < n_nodes; i++ ){
			int idx = nodes[i];
			if( y[idx]==UNLABELED )
				b[k*total+idx] = 0; //unlabeled data
			else if( y[idx]==k )
				b[k*total+idx] = 1;
			else
				b[k*total+idx] = -1;			
		}
	}

	tree_laplacian_solver(tree, nodes, n_nodes, total, c, n_class, b);

	for( i=0; i<n_nodes; i++ ) y[nodes[i]] = 0;
	for( k=1; k<n_class; k++){	
		REAL *ptr = b + k*total;
		for( i=0; i<n_nodes; i++){
			int idx = nodes[i];
			if( ptr[idx]>b[idx] ){
				y[idx] = k;
				b[idx] = ptr[idx];
			}
		}
	}

	free(b);
}



void majority_vote(int n_node, int n_class, int n_tree, int **rlts, int *y)
{
	int i, j, cmax;

	bool *is_labeled = (bool *)malloc(n_node*sizeof(bool));
	for( i=0; i<n_node; i++ ) is_labeled[i] = (rlts[0][i]==-1) ? false:true;
	int *hist = (int *)malloc(n_class*sizeof(int));
	for( i=0; i<n_node; i++ ){
		if( !is_labeled[i] ){
			y[i] = -1;
			continue;
		}
		memset(hist, 0, n_class*sizeof(int));
		for( j=0; j<n_tree; j++ ){
			hist[rlts[j][i]]++;
		}
		cmax = 0;
		for( j=1; j<n_class; j++ ){
			if( hist[j]>hist[cmax] ) cmax = j;			
		}
		y[i] = cmax;
	}

	free(hist);
	free(is_labeled);
}

void llgc_on_tree(struct Edge *graph, int n_nodes, int n_edges, int *y, int n_class, int tree_type, int n_trees, REAL alpha)
{	
	int i, j;
	//build the regularization matrix C = diag(c)
	REAL *c = (REAL *)malloc(n_nodes*sizeof(REAL));
	for( i=0; i<n_nodes; i++){	
		if( y[i]==UNLABELED )
			c[i] = 0; // no regularization on unlabeled data
		else
			c[i] = alpha;
	}

	int *buf, **rlts;
	rlts = (int **)malloc(n_trees*sizeof(int *));
	buf = (int *)malloc(n_trees*n_nodes*sizeof(int));
	if( !rlts || !buf ){	
		printf("error in allocating memory.\n");
		exit(1);
	}

	for( i=0; i<n_trees; i++ ){
		rlts[i] = buf + i*n_nodes;
		memcpy(rlts[i], y, n_nodes*sizeof(int));
	}

	struct Edge *e_com;
	e_com = (struct Edge *)malloc(n_edges*sizeof(struct Edge));
	int n_components, *e_com_sz, *n_com, *n_com_sz;	
	e_com_sz = (int *)malloc(3*n_nodes*sizeof(int));
	n_com = e_com_sz + n_nodes;
	n_com_sz = n_com + n_nodes;
	if( !e_com || !e_com_sz ){	
		printf("error in allocating memory.\n");
		exit(1);
	}

	n_components = detect_disconnected_components(graph, n_nodes, n_edges, e_com, e_com_sz, n_com, n_com_sz);
	if( tree_type==0 || tree_type==1 ) 
		for(i=0; i<n_edges; i++) e_com[i].w = -e_com[i].w;

	bool *is_labeled = (bool *)malloc(n_components*sizeof(bool));
	for(i=0; i<n_components; i++) is_labeled[i] = false;
	int s = 0;
	for(i=0; i<n_components; i++){
		for(j=0; j<n_com_sz[i]; j++){
			if( y[n_com[s+j]] != -1 ){
				is_labeled[i] = true;
				break;
			} 
		}
		s += n_com_sz[i];
	}
	
	int ii, jj, pos1, pos2, r;
	for( jj=0; jj<n_trees; jj++ ){
	
		pos1 = 0;
		pos2 = 0;	

		for( ii=0; ii<n_components; ii++ ){
			if( !is_labeled[ii] || n_com_sz[ii]==1){
				pos1 += e_com_sz[ii];
				pos2 += n_com_sz[ii];
				continue;
			} 
//clock_t t1, t2, t3;
//t1 = clock();
			//a. generate a spanning tree for component ii
			struct Edge *tree_graph = (struct Edge *)malloc((n_nodes-1)*sizeof(struct Edge));
			switch(tree_type){
				case 0:
					mini_spanning_tree(e_com+pos1, n_com+pos2, n_com_sz[ii], e_com_sz[ii], tree_graph, n_nodes);
					for(i=0; i<n_com_sz[ii]-1; i++) tree_graph[i].w = -tree_graph[i].w;
					break;
				case 1:
					r = floor(((double)rand())/((double)RAND_MAX)*n_com_sz[ii]);
					shortest_path_tree(e_com+pos1, n_com+pos2, n_com_sz[ii], e_com_sz[ii], tree_graph, n_nodes, n_com[pos2+r]);
					for(i=0; i<n_com_sz[ii]-1; i++) tree_graph[i].w = -tree_graph[i].w;
					break;
				case 2:
					random_spanning_tree(e_com+pos1, n_com+pos2, n_com_sz[ii], e_com_sz[ii], tree_graph, n_nodes);
					break;
				default:
					printf("unknown tree type\n");
					exit(1);
					break;	
			}
//t2 = clock();			
			_llgc_on_tree(tree_graph, n_com+pos2, n_com_sz[ii], n_nodes, c, n_class, rlts[jj]);
			pos1 += e_com_sz[ii];
			pos2 += n_com_sz[ii];
			free(tree_graph);
//t3 = clock();
//printf("tree construction time: %.5f; transductive learning time: %.5f\n", ((double)(t2-t1))/CLOCKS_PER_SEC, ((double)(t3-t2))/CLOCKS_PER_SEC);
		}	
			
	}

	majority_vote(n_nodes, n_class, n_trees, rlts, y);

	free(c);
	free(is_labeled);
	free(e_com);
	free(e_com_sz);
	free(buf);
	free(rlts);
}