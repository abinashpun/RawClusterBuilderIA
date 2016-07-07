
// Requirements:
//
// the class type Hit needs to provide:
// 
//   an operator< that sorts by ix and then by iz
//   a function is_adjacent() which returns true if the argumment hit is adjacent
//   a function ...() that returns true if the argument Hit is far enough away
//      from this one to allow breaking out of the inner loop early
//

#include <vector>
#include <map>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

template<class Hit> 
int PHMakeGroups(std::vector<Hit>& hits, std::multimap<int, Hit>& groups) {

    using namespace boost;
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // Sort in eta; if same eta then sort in phi.
    std::sort(hits.begin(), hits.end());

    for(unsigned i = 0; i < hits.size(); i++) {
        for(unsigned j = i + 1; j < hits.size(); j++) {
            if (hits[i].is_adjacent(hits[j])) add_edge(i,j,G);
        }
        add_edge(i,i,G); // edge to itself?...
    }

    // Find the connections between the vertices of the graph (vertices are the rawhits, 
    // connections are made when they are adjacent to one another)
    std::vector<int> component(num_vertices(G));
    int num = connected_components(G, &component[0]);
    std::cout << "Found " << num << " groups of hits" << std::endl;
    std::cout << "\tsizeof(seeds) = " << hits.size() << std::endl;

    // Loop over the components(vertices) compiling a list of the unique
    // connections (ie clusters).
    std::set<int> comps; // Number of unique components
    for(unsigned i=0; i<component.size(); i++) {
        comps.insert(component[i]);
        groups.insert(std::make_pair(component[i], hits[i]));
    } 
    std::cout << "unique groups: " << comps.size();


    return 0; 
}
