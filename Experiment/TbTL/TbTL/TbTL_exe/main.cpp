#include "llgc_on_tree.h"

void load_graph(char *graph_file, int *n_node_p, int *n_edge_p, struct Edge **edges);
void load_labels(char *train_file, int n, int *n_classes, int **y, int **label_idx);
void output_rlts(char *rlts_file, int n, int *y, int *label_idx);

/**/
void exit_with_help()
{
	printf(
		"Useage: mtc [options] graph_file train_file output_file\n"
		"options:\n"
		"-t: select the type of spanning tree (default: 0)\n"
		"	0 -- minimum spanning tree\n"
		"	1 -- shortest path tree\n"
		"	2 -- random spanning tree\n"
		"-n: set the number of spanning trees (default: 1)\n"
		"-a: set the value of regulariztor (default: 1)\n"
		);
	exit(1);
}

void parse_command_line(int argc, char **argv, int *tree_type, int *n_tree, float *alpha, char *graph_file, char *train_file, char *output_file)
{
	int i;

	*tree_type = 0;
	*n_tree = 1;
	*alpha = 1;

	for( i=1; i<argc; i++ ){
	
		if(argv[i][0] != '-') break;
		if(++i >= argc) exit_with_help();

		switch(argv[i-1][1]){
		
			case 't':
				*tree_type = atoi(argv[i]);
				break;
			case 'n':
				*n_tree = atoi(argv[i]);
				break;

			case 'a':
				*alpha = atof(argv[i]);
				break;

			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
	
	if( i != argc-3 )
		exit_with_help();
	else{	
		strcpy(graph_file, argv[i]);
		strcpy(train_file, argv[i+1]);
		strcpy(output_file, argv[i+2]);
	}
}

int main(int argc, char **argv)
{
	int tree_type, n_tree, n_node, n_edge, n_class;
	float alpha;
	struct Edge *graph;
	int *y, *label_idx;
	char graph_file[LINE_LEN], train_file[LINE_LEN], output_file[LINE_LEN];

	parse_command_line(argc, argv, &tree_type, &n_tree, &alpha, graph_file, train_file, output_file);

	load_graph(graph_file, &n_node, &n_edge, &graph);
	load_labels(train_file, n_node, &n_class, &y, &label_idx);

	//alpha = alpha/n_node;

	clock_t tstart, tend;
	tstart = clock();	
	llgc_on_tree(graph, n_node, n_edge, y, n_class, tree_type, n_tree, alpha);
	tend = clock();
	printf("The running time: %f seconds\n", ((double)(tend-tstart))/CLOCKS_PER_SEC);

	output_rlts(output_file, n_node, y, label_idx);
	
	free(y);
	free(graph);

	return 0;
}


//假设node idx是从0开始的连续整数
void load_graph(char *graph_file, int *n_node_p, int *n_edge_p, struct Edge **edges)
{
	FILE *fd = fopen(graph_file, "r");
	if(!fd){
		printf("error in reading file: %s\n", graph_file);
		exit(1);
	}

	int n_node = -1;
	int n_edge = 0;
	int buf_sz = 0;
	int block_sz = 1000000;
	struct Edge *buf = NULL;
	char line[LINE_LEN];
	while( fgets(line, LINE_LEN, fd) )
	{
		if( n_edge == buf_sz ){
			buf_sz += block_sz;
			buf = (struct Edge *)realloc( buf, buf_sz*sizeof(struct Edge) );
			if( !buf ){
				printf("error in (re)allocating memory when load graph\n");
				exit(1);
			}
		}
		int x, y;
		float w;
		sscanf(line, "%d,%d,%f", &x, &y, &w);
		buf[n_edge].x = x;
		buf[n_edge].y = y;
		buf[n_edge++].w = w;

		if( n_node<x ) n_node = x;
		if( n_node<y ) n_node = y;
	}

	fclose(fd);
	*n_node_p = n_node+1;
	*n_edge_p = n_edge;
	*edges = buf;

	return;
}

//假设label是从0开始的连续整数
void load_labels(char *train_file, int n, int *n_classes, int **y, int **label_idx)
{
	FILE *fd = fopen(train_file, "r");
	if(!fd){
		printf("error in reading label file.\n");
		exit(1);
	}

	int *labels = (int *)malloc(2*n*sizeof(int));	
	if( !labels ){
		printf("error in allocating memory when load training file\n");
		exit(1);
	}
	for( int i=0; i<2*n; i++ ) labels[i] = -1;

	int K = -1;
	int idx, c;
	while( fscanf(fd, "%d:%d\n", &idx, &c)==2 ){
		labels[idx] = c;
		labels[n+idx] = 1;
		if( c>K ) K = c;
	}

	*y = labels;
	*label_idx = labels+n;
	*n_classes = K+1;
	fclose(fd);
	return;
}

void output_rlts(char *rlts_file, int n, int *y, int *label_idx)
{
	FILE *fd = fopen( rlts_file, "w" );
	for( int i=0; i<n; i++ )
		if(label_idx[i]==-1) fprintf(fd, "%d:%d\n", i, y[i]);

	fclose(fd);
}
