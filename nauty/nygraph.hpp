// nygraph.hpp
// a C++ wrapper for nauty graphs

// Cytonaut - Python interface to Nauty
// Written in 2015 by Gregory J. Puleo, puleo@illinois.edu

// To the extent possible under law, the author(s) have dedicated all
// copyright and related and neighboring rights to this software to the
// public domain worldwide. This software is distributed without any
// warranty.

#define MAXN WORDSIZE
#include <nauty.h>
#include <naugroup.h>
#include <gtools.h>
#include <vector>
#include <string>
#include <cstdio>

#define WORKSPACE_FACTOR 66

/*
typedef std::vector<int> permu;
std::vector<permu> *automorphisms_target=0;

void my_write_autom(permutation *p, int n)
{
    permu p2(n);
    for(int i=0; i < n; i++)
    {
	p2[i] = p[i];
    }
    automorphisms_target->push_back(p2);
}
*/
extern "C" {
  
void testx();
}

class NyGraph
{
public:
    int num_nodes;
    int setwords_in_set;
    std::vector<setword> graph_data;
    std::vector<setword> canong;    
    std::vector<int> canon_labeling;
    std::vector<int> ptn;
    std::vector<int> orbits;

    //std::vector<permu> automorphisms;

    NyGraph(int _num_nodes) :
	num_nodes( _num_nodes ),
	setwords_in_set( (num_nodes + WORDSIZE - 1) / WORDSIZE ),
	graph_data( _num_nodes ),
	canon_labeling( _num_nodes ),
	ptn( _num_nodes )
    {};

    inline void add_edge(int i, int j)
    {
	ADDELEMENT(GRAPHROW(&(graph_data)[0],i,setwords_in_set), j);
	ADDELEMENT(GRAPHROW(&(graph_data)[0],j,setwords_in_set), i);
    }

    void do_nauty(bool autom)
    {
	DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = TRUE;
	options.defaultptn = FALSE; // we now use coloring data given by client
	if(autom)
	{
	    options.userautomproc = groupautomproc;
	    options.userlevelproc = grouplevelproc;
	}
	statsblk stats;
	int worksize = num_nodes * setwords_in_set * WORKSPACE_FACTOR;
	canon_labeling.resize(num_nodes);
	orbits.resize(num_nodes);
	canong.resize(num_nodes);
	std::vector<setword> workspace(worksize);

	nauty(&(graph_data)[0], &(canon_labeling)[0], &(ptn)[0], NULL, &(orbits)[0], &options, &stats,
	      &(workspace)[0], worksize, setwords_in_set, num_nodes, &(canong)[0]);
    }

    /*
    void get_automorphisms()
    {
	grouprec *group;
	group = groupptr(FALSE);
	makecosetreps(group);
	automorphisms_target = &automorphisms;
	allgroup(group, my_write_autom);
	automorphisms_target = 0;
    }
*/
    std::string canon_string()
    {
	// assumes canong is already known
	return std::string(ntog6(&(canong[0]), setwords_in_set, num_nodes));
    }
};
