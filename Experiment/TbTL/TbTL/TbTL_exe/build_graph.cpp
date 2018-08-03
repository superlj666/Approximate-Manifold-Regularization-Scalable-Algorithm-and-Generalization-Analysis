#include "llgc_on_tree.h"

//build an undirected graph
void build_undirected_graph(struct Edge *graph, int n_nodes, int n_edges, struct Edge_ **g, struct Edge_ *pool)
{
	memset(g, 0, n_nodes*sizeof(struct Edge_ *));
	for(int i=0; i<n_edges; i++){
		int j = 2*i;
		pool[j].idx = graph[i].x;
		pool[j].w = graph[i].w;
		pool[j].next = g[graph[i].y];
		g[graph[i].y] = pool + j;
		j++;
		pool[j].idx = graph[i].y;
		pool[j].w = graph[i].w;
		pool[j].next = g[graph[i].x];
		g[graph[i].x] = pool + j;
	}

}

void build_directed_graph(struct Edge *graph, int n_nodes, int n_edges, struct Edge_ **g, struct Edge_ *pool)
{
	memset(g, 0, n_nodes*sizeof(struct Edge_ *));
	for(int i=0; i<n_edges; i++){
		pool[i].idx = graph[i].y;
		pool[i].w = graph[i].w;
		pool[i].next = g[graph[i].x];
		g[graph[i].x] = pool + i;
	}
}

void reverse_directed_graph(struct Edge *graph, int n_nodes, int n_edges, struct Edge_ **g, struct Edge_ *pool)
{
	memset(g, 0, n_nodes*sizeof(struct Edge_ *));
	for(int i=0; i<n_edges; i++){
		pool[i].idx = graph[i].x;
		pool[i].w = graph[i].w;
		pool[i].next = g[graph[i].y];
		g[graph[i].y] = pool + i;
	}
}