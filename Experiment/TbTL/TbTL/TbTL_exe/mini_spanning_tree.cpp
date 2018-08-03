#include "llgc_on_tree.h"
#include <vector>
#include <algorithm>

int cmp(const void *i, const void *j)
{
	if(((struct Edge *)i)->w < ((struct Edge *)j)->w) return -1;
	if(((struct Edge *)i)->w == ((struct Edge *)j)->w) return 0;
	
	return 1;
}

struct edge_less {
    inline bool operator ()(struct Edge const& a, struct Edge const& b) const 
	{
		return a.w<b.w;
    }
};


int *s;//union-find set
int *rank;

int find(int x)
{
	while( x!=s[x] ) x = s[x];		
	return x;
}

void Union(int x, int y)
{
	if(rank[x]>rank[y])
		s[y] = x;
	else{
		s[x] = y;
		if(rank[x]==rank[y]) rank[y]++;
	}
}


//input: total is the max num of node idx.
void mini_spanning_tree( struct Edge *graph, int *nodes, int n_nodes, int n_edges, struct Edge *tree, int total )
{
	int i, cnt;	

	std::vector<Edge> g(graph, graph+n_edges);
	std::sort(g.begin(), g.end(), edge_less()); 
	for( i = 0; i < n_edges; i++ ) graph[i] = g[i];

	s = (int *)malloc(2*total*sizeof(int));
	rank = s + total;
	for(i=0; i<total; i++){
		s[i] = i;
		rank[i] = 0;
	}
	i = 0;
	cnt = 0;
	while(cnt<n_nodes-1){
		int x = find(graph[i].x);
		int y = find(graph[i].y); 
		if( x != y ){
			tree[cnt++] = graph[i];	
			Union(x, y);
		}
		i++;
	}

	free(s);
}

