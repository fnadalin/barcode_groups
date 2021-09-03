#ifndef GRAPH_H_
#define GRAPH_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>

#include "Reads.h"
#include "Hash.h"

class Graph
{
    // container for graph structure
public:

    static const uint32_t MAX_NEIGHBORS = 100;

    const double MIN_FRAC_IN_DEGREE = 0.8; // fraction of hub's neighbors inside the group (i.e. in degree / total degree)

    enum color_type { WHITE, GRAY, BLACK };

    static bool verbose;
    inline static void set_verbose () { verbose = true; }

    Graph (uint32_t max_distance, double max_read_frac_neighbor);
    ~Graph () {}

    /* Build the graph starting from the sequences
     * - use the hamming distance between sequences
     * - use fingerprints to quickly find neighbors
     * A Graph is a vector of Nodes.
     * Each Node contains a list of indexes of elements of the Graph (neighbors) and a list of distances (associated to the neighbors).
     */
    void                          build (Hash & H);

    // number of nodes
    inline size_t                 size () { return _graph.size(); }

    // set all nodes to white
    void                          reset_color ();

    // set all nodes to black
    void                          set_color ();

    // check whether all nodes are black
    bool                          all_black ();

    class Node
    {
    public:

        Node (uint32_t idx);
        Node (uint32_t idx, uint32_t count);
        ~Node();

        std::vector<uint32_t>     neighbors;
        std::vector<uint32_t>     distances; // hamming distance
        uint32_t                  node_idx; // node index
        uint32_t                  parent;
        uint32_t                  count;
        Graph::color_type         color = WHITE;
    };

    // return a vector with the degree of each node in _graph
    std::vector <uint32_t>        compute_degree ();
    
    // return a vector v of size max_degree+1 with v[i] = vector of nodes ID with degree i
    std::vector < std::vector <uint32_t> > degree_sets ();

    // return a vector v of size max_degree+1 with v[i] = map of nodes ID with degree i
    std::vector < std::map <uint32_t,uint32_t> > degree_maps ();

    // print stats on the graph
    void                          print_graph_attributes (const Hash & H, const std::string & out);

    // compute the stars in the network, iteratively
    // hub of the star = node with max degree
    // the neighbors of the neighbors of the hub are also neighbors of the hub 
//    void                          count_stars (const Hash & H);

    // compute the stars in the network, iteratively
    // hub of the star = node with max degree, parent of all the other nodes in the star
    // print the number of stars and the number of edges linking two distinct stars 
    // the aim would be to minimize this number 
    void                          create_stars ();

    // return a map containing the read count (first) and the node indexes with that read count (= read ID, second)
    std::map < uint32_t, std::vector <uint32_t> >  count_map (const Hash & H);

    // compute the stars in the network, iteratively
    // hub of the star = node with max count, parent of all the other nodes in the star
    // the aim would be to minimize the links between high count nodes
    void                          create_stars_with_counts (const Hash & H);

    // collect the groups into sets labelled with the hub node 
    std::map < uint32_t, std::map <uint32_t, uint32_t> > & collect_stars ();

    // sort the groups by total count
    std::map < uint32_t, std::vector <uint32_t> > & sort_stars (const std::map < uint32_t, std::map <uint32_t, uint32_t> > & groups);

    // compute the minimum number of reads in a group to guarantee that all groups satisfy hub_in_degree >= MIN_FRAC_IN_DEGREE * hub_total_degree
    uint32_t min_group_count_sel (const std::map < uint32_t, std::vector <uint32_t> > & groups_sorted, const std::map < uint32_t, std::map <uint32_t, uint32_t> > & groups);

    // print degree and count for each node
//    void                          print_node_attributes (const Hash & H, const std::string & out);

    // print node count for each edge (i.e. node pair)
//    void                          print_edge_attributes (const Hash & H, const std::string & out);

    // print read count for the top and the 2nd-to-top neighbor of each node 
//    void                          print_neighbor_attributes (const Hash & H, const std::string & out);

    // print a vector v of size max_degree+1 with v[i] = number of nodes with degree i
//    void                          print_degree_frequency (const std::string & out);

    // print the barcode groups
    void                          print_stars (const Hash & H, const std::string & out) const;

    // print the barcode groups with count info
//    void                          print_stars_with_counts (const Hash & H, const std::string & out);

    // print node and parent information (parent = hub of the group, node = element of the group)
//    void                          print_stars_info (const Hash & H, const std::string & out);

    // print the links between barcodes belonging to different groups
//    void                          print_out_of_stars (const Hash & H, const std::string & out);

    // print a vector v of size max_star_size+1 with v[i] = number of stars with i elements
//    void                          print_group_size_frequency (const std::string & out);

private:

    std::vector <Node>            _graph;
    uint32_t                      _max_distance;
    double                        _max_read_frac_neighbor;

};

#endif /* GRAPH_H_ */
