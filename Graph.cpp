
#include "Graph.h"

bool Graph::verbose = false;

Graph::Graph (uint32_t max_distance, double max_read_frac_neighbor)
{
    _max_distance = max_distance;
    _max_read_frac_neighbor = max_read_frac_neighbor;
}

void
Graph::build (Hash & H)
{
    _graph.reserve(H.numReads());

    // Assign node properties
    for (auto i = 0; i < H.numReads(); ++i)
    {
         auto c = H.getRead(i).count();
        _graph.emplace_back(i, c); // NEW: add read count
    }

    // Compute distances and add a link only if it is
    // bounded by _max_distance
    for (auto iter = H._HASHcounter.begin(); iter != H._HASHcounter.end(); ++iter)
    {
        auto start = *iter; // first pos. in _HASHvalues with the current fingerprint
        auto jter = iter; ++jter;
        auto end = (jter != H._HASHcounter.end()) ? *jter : H._HASHvalues.size(); // past-to-last pos. in _HASHvalues with the current fingerprint

        for (auto i = start; i < end; ++i)
        {
            if (i >= H._HASHvalues.size()) std::cerr << "i >= H._HASHvalues.size() [" << i << "," << H._HASHvalues.size() << "]" << std::endl;
            uint32_t r1 = H._HASHvalues.at(i).read_id;
            uint32_t f1 = H._HASHvalues.at(i).block_id;
            for (auto j = i+1; j < end; ++j)
            { 
                if (j >= H._HASHvalues.size()) std::cerr << "j >= H._HASHvalues.size() [" << j << "," << H._HASHvalues.size() << "]" << std::endl;
                uint32_t f2 = H._HASHvalues.at(j).block_id;
                if (f1 != f2) continue; // only consider the sequences witn the same fingerprint on the same block

                uint32_t r2 = H._HASHvalues.at(j).read_id;
                uint32_t dist = H.ComputeHammingDistance (r1, r2);
                if (dist > _max_distance) continue;

                // if j no not exist among the neighbors of i yet, then add it and add the distance as well
                // if j exists already, then check the distance
                // if the distance is larger, then update it, otherwise do nothing
                if (r1 >= _graph.size()) std::cerr << "r1 >= _graph.size() [" << r1 << "," << _graph.size() << "]" << std::endl;
                auto n1 = _graph.at(r1).neighbors.begin();
                auto d1 = _graph.at(r1).distances.begin();
                while (n1 != _graph.at(r1).neighbors.end() and *n1 != r2)
                {
                    ++n1;
                    ++d1;
                }
                if (r2 >= _graph.size()) std::cerr << "r2 >= _graph.size() [" << r2 << "," << _graph.size() << "]" << std::endl;
                if (n1 == _graph.at(r1).neighbors.end())
                {
                    _graph.at(r1).neighbors.push_back(r2);
                    _graph.at(r1).distances.push_back(dist);
                    _graph.at(r2).neighbors.push_back(r1);
                    _graph.at(r2).distances.push_back(dist);
                } 
                else
                {
                    if (*d1 <= dist) continue;
                    
                    *d1 = dist;
                    auto n2 = _graph.at(r1).neighbors.begin();
                    auto d2 = _graph.at(r1).distances.begin();
                    while (n2 != _graph.at(r2).neighbors.end() and *n2 != r1)
                    {
                        ++n2;
                        ++d2;
                    }
                    if (n2 == _graph.at(r2).neighbors.end()) 
                    {
                        std::cerr << "the graph is not undirected!" << std::endl;
                        return;
                    }
                    *d2 = dist;
                }
            }
        }
    }
}

Graph::Node::Node (uint32_t idx)
{
    neighbors.reserve(MAX_NEIGHBORS);
    distances.reserve(MAX_NEIGHBORS);
    node_idx = idx;
    parent = idx;
    count = 1;
    color = WHITE;
}

Graph::Node::Node (uint32_t idx, uint32_t c)
{
    neighbors.reserve(MAX_NEIGHBORS);
    distances.reserve(MAX_NEIGHBORS);
    node_idx = idx;
    parent = idx;
    count = c;
    color = WHITE;
}

Graph::Node::~Node ()
{
    neighbors.clear();
    distances.clear();
}

void
Graph::reset_color ()
{
    for (std::vector <Node>::iterator N = _graph.begin(); N != _graph.end(); ++N)
        N->color = WHITE;
}

void
Graph::set_color ()
{
    for (std::vector <Node>::iterator N = _graph.begin(); N != _graph.end(); ++N)
        N->color = BLACK;
}

bool
Graph::all_black ()
{
    bool all_black = true;
    for (std::vector <Node>::iterator N = _graph.begin(); N != _graph.end(); ++N)
        all_black = (all_black and (N->color == BLACK));
    return all_black;
}

std::vector <uint32_t>        
Graph::compute_degree ()
{
    std::vector <uint32_t> v;
    v.reserve(_graph.size());
    for (auto iter = _graph.begin(); iter != _graph.end(); ++iter)
        v.push_back(iter->neighbors.size());

    return v;
}

std::vector < std::vector <uint32_t> >        
Graph::degree_sets ()
{
    std::vector <uint32_t> v = compute_degree();
    std::vector < std::vector <uint32_t> > w;
    auto max_iter = std::max_element(v.begin(),v.end());
    if (max_iter == v.end()) return w;

    size_t max = *max_iter; 
    w.resize(max+1);
    auto i = 0;
    for (auto iter = v.begin(); iter != v.end(); ++iter) 
    {
        if (w.size() == w.capacity())
            w.reserve(w.size() + MAX_NEIGHBORS);
        w.at(*iter).push_back(i);
        ++i;
    }
    
    return w;
}

std::vector < std::map <uint32_t,uint32_t> >        
Graph::degree_maps ()
{
    std::vector <uint32_t> v = compute_degree();
    std::vector < std::map <uint32_t, uint32_t> > w;
    auto max_iter = std::max_element(v.begin(),v.end());
    if (max_iter == v.end()) return w;
 
    size_t max = *max_iter; 
    w.resize(max+1);
    auto i = 0;
    for (auto iter = v.begin(); iter != v.end(); ++iter) 
    {
        w.at(*iter).insert(std::make_pair(i,i));
        ++i;
    }
    
    return w;
}

void                        
Graph::print_graph_attributes (const Hash & H, const std::string & out_dir)
{
    // node attributes
    std::string outname(out_dir); outname.append("/node_attr.tsv");
    std::ofstream outfile; outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "id	num_neighbors	counts" << std::endl;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        outfile << H.getRead(ii->node_idx).id() << "	"
                << ii->neighbors.size() << "	"
                << H.getRead(ii->node_idx).count() << std::endl;
    }
    outfile.close();

    // edge attributes
    outname = out_dir; outname.append("/edge_attr.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "id1	id2	count1	count2" << std::endl;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        for (auto jj = ii->neighbors.begin(); jj != ii->neighbors.end(); ++jj)
        {
            if (ii->node_idx < _graph.at(*jj).node_idx)
            {
                auto idx1 = (H.getRead(ii->node_idx).count() >= H.getRead(_graph.at(*jj).node_idx).count()) ? ii->node_idx : _graph.at(*jj).node_idx;
                auto idx2 = (H.getRead(ii->node_idx).count() >= H.getRead(_graph.at(*jj).node_idx).count()) ? _graph.at(*jj).node_idx : ii->node_idx;
                outfile << H.getRead(idx1).id() << "	"
                        << H.getRead(idx2).id() << "	"
                        << H.getRead(idx1).count() << "	"
                        << H.getRead(idx2).count() << std::endl;
            }
        }
    }
    outfile.close();

    // neighbor attributes
    outname = out_dir; outname.append("/neighbor_attr.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "id	top_neigh_count	top2_neigh_count" << std::endl;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        if (ii->neighbors.size() > 1) 
        {
            // store the read count for all neighbors of ii
            std::vector<uint32_t> neigh_counts;
            neigh_counts.reserve(ii->neighbors.size());
            for (auto jj = ii->neighbors.begin(); jj != ii->neighbors.end(); ++jj)
                neigh_counts.push_back(H.getRead(*jj).count());

            // find the top and the 2nd top read count neighbors
            auto it_max_count = std::max_element(neigh_counts.begin(), neigh_counts.end());
            uint32_t max_count = *it_max_count;
            neigh_counts.erase(it_max_count);
            auto it_sec_max_count = std::max_element(neigh_counts.begin(), neigh_counts.end());
            uint32_t sec_max_count = *it_sec_max_count;

            // check that the two neighbors of ii are not mutual neighbors
            auto max_pos = it_max_count - neigh_counts.begin();
            auto sec_max_pos = it_sec_max_count - neigh_counts.begin();
            if (max_pos <= sec_max_pos) ++sec_max_pos;
            auto jj = ii->neighbors.at(max_pos);
            auto kk = ii->neighbors.at(sec_max_pos);

            // print the non-mutual top count neighbors of ii
            std::vector<uint32_t> nn = _graph.at(jj).neighbors;
            if (std::find(nn.begin(), nn.end(), kk) == nn.end())
                outfile << H.getRead(ii->node_idx).id() << "	"
                        << max_count << "	"
                        << sec_max_count << std::endl;
        }
    }
    outfile.close();

    // degree frequency
    outname = out_dir; outname.append("/degree_freq.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "degree	num_nodes" << std::endl;
    std::vector < std::vector<uint32_t> > w = degree_sets();
    auto i = 0;
    for (auto iter = w.begin(); iter != w.end(); ++iter) 
        outfile << i++ << "	" << iter->size() << std::endl;
    outfile.close();

}

void
Graph::create_stars ()
{
    // set all nodes to WHITE = not visited
    reset_color();

    auto deg_nodes = compute_degree(); // idx = node idx, value = degree
    auto degrees = degree_maps(); // idx = degree, value = map of nodes with the same degree (specified by idx)

    // scan nodes by non-increasing degree order
    uint32_t stars = 0;
    uint32_t out_of_stars = 0; // initialize the count of the edges linking two distict stars  
    for (auto nodes = degrees.rbegin(); nodes != degrees.rend(); ++nodes) 
    {
        while (not nodes->empty())
        {
            auto ii = nodes->begin()->first; // first node in the map
            if (ii >= _graph.size()) 
            {
                std::cerr << "Node " << ii << " is outside of node idx span" << std::endl;
                continue;
            }
            if (_graph.at(ii).color == BLACK) continue; 

            // visit node ii
            _graph.at(ii).color = GRAY;

            // visit ii's neighbors jj
            std::vector<uint32_t> ne;
            ne.reserve(_graph.at(ii).neighbors.size());
            for (auto jj = _graph.at(ii).neighbors.begin(); jj != _graph.at(ii).neighbors.end(); ++jj) 
            {
                if (_graph.at(*jj).color == BLACK) continue; // the edge (ii,jj) has already been counted in out_of_stars
                _graph.at(*jj).color = GRAY;
                _graph.at(*jj).parent = ii; // THIS CONTAINS THE STAR/GROUP INFORMATION!!!
                ne.push_back(*jj);
            }

            // scan the neighbors of jj
            std::vector <uint32_t> nne;
            nne.reserve(MAX_NEIGHBORS);
            for (auto jj = ne.begin(); jj != ne.end(); ++jj) 
            {
                auto v = _graph.at(*jj).neighbors;
                for (auto kk = v.begin(); kk != v.end(); ++kk)
                { 
                    // new neighbor
                    if (_graph.at(*kk).color == WHITE and std::find(nne.begin(), nne.end(), *kk) == nne.end())
                        nne.push_back(*kk);
                }
            }
            out_of_stars += nne.size();

            // update the degree of the neighbors of the neighbors
            for (auto kk = nne.begin(); kk != nne.end(); ++kk)
            {
                if (*kk >= deg_nodes.size()) 
                {
                    std::cerr << "Node " << *kk << " is outside of deg nodes span" << std::endl;
                    continue;
                }
                auto current_degree = deg_nodes.at(*kk);
                uint32_t edges_to_remove = 0;
                auto v = _graph.at(*kk).neighbors;
                for (auto hh = v.begin(); hh != v.end(); ++hh)
                    edges_to_remove += (_graph.at(*hh).color == GRAY);
                if (edges_to_remove > current_degree)
                {
                    std::cerr << "The number of edges to remove from node " << *kk << " is greater than its degree (" 
                              << edges_to_remove << " > " << current_degree << ")" << std::endl;
                    continue;
                }
                auto new_degree = current_degree - edges_to_remove;
                // find the current slot and remove kk
                auto nodes_n = degrees.begin() + current_degree;
                auto iter = nodes_n->find(*kk);
                if (iter == nodes_n->end()) 
                {
                    std::cerr << "Node " << *kk << " was not found on degree slot " << current_degree << std::endl;
                    continue;
                }
                nodes_n->erase(iter);
                // find the new slot and add kk to it
                nodes_n = degrees.begin() + new_degree;
                nodes_n->insert(std::make_pair(*kk,*kk));
                deg_nodes.at(*kk) = new_degree; // IMPORTANT!! otherwise the node won't be found in the correct degrees slot
            }

            // remove ii from the degree slot
            nodes->erase(nodes->begin()); 

            // remove ii's neighbors from their degree slot
            for (auto jj = ne.begin(); jj != ne.end(); ++jj)
            {
                if (*jj >= deg_nodes.size()) 
                {
                    std::cerr << "Node " << *jj << " is outside of deg nodes span" << std::endl;
                    continue;
                }
                if (deg_nodes.at(*jj) >= degrees.size())
                {
                    std::cerr << "Node " << *jj << " degree is greater than maximum degree (" 
                              << deg_nodes.at(*jj) << " >= " << degrees.size() << ")" << std::endl;
                }
                auto nodes_n = degrees.begin() + deg_nodes.at(*jj); // iterator to the map containing the nodes with degree = degree of j
                auto iter = nodes_n->find(*jj);
                if (iter == nodes_n->end()) 
                {
                    std::cerr << "Node " << *jj << " was not found on degree slot " << deg_nodes.at(*jj) << std::endl;
                    continue;
                }
                nodes_n->erase(iter);
            }

            // end the visit of ii and jj
            // N.B.: kk are NOT marked as visited yet
            _graph.at(ii).color = BLACK;
            for (auto jj = ne.begin(); jj != ne.end(); ++jj) 
                _graph.at(*jj).color = BLACK;

            ++stars;
        }
    }

    std::cout << "number of groups: " << stars << std::endl;
    std::cout << "number of edges connecting two groups: " << out_of_stars << std::endl;

}

std::map < uint32_t, std::vector <uint32_t> >   
Graph::count_map (const Hash & H)
{
    std::map < uint32_t, std::vector <uint32_t> > m;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        auto cc = H.getRead(ii->node_idx).count();
        if (m.find(cc) == m.end()) 
        {
            std::vector <uint32_t> v;
            m.insert(std::make_pair(cc, v));
        }
        m.find(cc)->second.push_back(ii->node_idx);
    }

    return m;
}

void
Graph::create_stars_with_counts (const Hash & H)
{
    // set all nodes to WHITE = not visited
    reset_color();

    auto counts = count_map(H); // first = count, second = vector of node indexes

    // scan nodes by non-increasing degree order
    uint32_t stars = 0;
    uint32_t out_of_stars = 0; // initialize the count of the edges linking two distict stars  
    for (auto cc = counts.rbegin(); cc != counts.rend(); ++cc) 
    {
//        std::cout << "read count of current node: " << cc->first << std::endl;
        // scan node indexes
        for (auto ii = cc->second.begin(); ii != cc->second.end(); ++ii) 
        { 
            if (*ii >= _graph.size()) 
            {
                std::cerr << "Node " << *ii << " is outside of node idx span" << std::endl;
                continue;
            } 
            if (_graph.at(*ii).color == BLACK) continue; 

            // visit node ii
            _graph.at(*ii).color = GRAY;

            // visit ii's neighbors jj
            std::vector<uint32_t> ne;
            ne.reserve(_graph.at(*ii).neighbors.size());
            for (auto jj = _graph.at(*ii).neighbors.begin(); jj != _graph.at(*ii).neighbors.end(); ++jj) 
            {
                if (_graph.at(*jj).color == BLACK) continue; // the edge (ii,jj) has already been counted in out_of_stars
                // NEW: check that the count of the neighbor is below _max_read_frac_neighbor fraction 
                if (_max_read_frac_neighbor * _graph.at(*ii).count >= _graph.at(*jj).count)
                {
                    _graph.at(*jj).color = GRAY;
                    _graph.at(*jj).parent = *ii; // THIS CONTAINS THE STAR/GROUP INFORMATION!!!
                    ne.push_back(*jj);
                } // else {
//                    std::cout << _max_read_frac_neighbor << "*" << _graph.at(*ii).count << " = " << _max_read_frac_neighbor * _graph.at(*ii).count << " < " << _graph.at(*jj).count << std::endl;
//                }
            }

            // scan the neighbors of jj
            std::vector <uint32_t> nne;
            nne.reserve(MAX_NEIGHBORS);
            for (auto jj = ne.begin(); jj != ne.end(); ++jj) 
            {
                auto v = _graph.at(*jj).neighbors;
                for (auto kk = v.begin(); kk != v.end(); ++kk)
                { 
                    // new neighbor
                    if (_graph.at(*kk).color == WHITE and std::find(nne.begin(), nne.end(), *kk) == nne.end())
                        nne.push_back(*kk);
                }
            }
            out_of_stars += nne.size();

            // end the visit of ii and jj
            // N.B.: kk are NOT marked as visited yet
            _graph.at(*ii).color = BLACK;
            for (auto jj = ne.begin(); jj != ne.end(); ++jj) 
                _graph.at(*jj).color = BLACK;

            ++stars;
        }
    }

    std::cout << "number of groups: " << stars << std::endl;
    std::cout << "number of edges connecting two groups: " << out_of_stars << std::endl;

}

void
Graph::print_stars (const Hash & H, const std::string & out_dir) const
{
    // collect stars
    std::map < uint32_t, std::map <uint32_t, uint32_t> > groups; // key = parent idx, value = (idx,idx) 
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        if (groups.find(ii->parent) == groups.end())
        {
            std::map <uint32_t, uint32_t> m;
            groups.insert(std::make_pair(ii->parent, m));
        }
        auto map_iter = groups.find(ii->parent);
        map_iter->second.insert(std::make_pair(ii->node_idx,ii->node_idx));
    }

    // sort stars
    std::map < uint32_t, std::vector <uint32_t> > groups_sorted; // key = group count, value = list of parents identifying the groups 
    for (auto ii = groups.begin(); ii != groups.end(); ++ii)
    {
        auto c = 0;
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
            c += H.getRead(jj->first).count();
        if (groups_sorted.find(c) == groups_sorted.end())
        {
            std::vector <uint32_t> v;
            groups_sorted.insert(std::make_pair(c,v));
        }
        groups_sorted.find(c)->second.push_back(ii->first);
    }

    // select min group count
    uint32_t c;
    bool found = false;
    for (auto ii = groups_sorted.rbegin(); ii != groups_sorted.rend(); ++ii)
    {
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
        {
            double tot_degree = _graph.at(*jj).neighbors.size();
            double in_degree = groups.find(*jj)->second.size() - 1;
            found = (in_degree < MIN_FRAC_IN_DEGREE * tot_degree);
            if (found) break;
        }
        if (found) break;
        c = ii->first;
    }
    std::cout << "Min group count: " << c << std::endl;

    // compute group size frequency
    uint32_t max_size_all = 0;
    for (auto ii = groups.begin(); ii != groups.end(); ++ii)
    {
//        if (max_size_all < _graph.at(ii->first).neighbors.size())
//            max_size_all = _graph.at(ii->first).neighbors.size();
        if (max_size_all < ii->second.size())
            max_size_all = ii->second.size();
    } 
    std::vector <uint32_t> group_sizes_all; group_sizes_all.resize(max_size_all, 0);
    for (auto ii = groups.begin(); ii != groups.end(); ++ii)
    {
        if (ii->second.size()-1 >= group_sizes_all.size())
        { 
            std::cerr << "ii->second.size()-1 >= group_sizes_all.size() (" << (ii->second.size()-1) << " > " << group_sizes_all.size() << ")" << std::endl;
            return;
        }
        group_sizes_all.at(ii->second.size()-1)++;
    }

    // compute selected group size frequency
    uint32_t max_size = 0;
    for (auto ii = groups_sorted.rbegin(); ii != groups_sorted.rend(); ++ii)
    {
        if (ii->first < c) break;
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
        {
            auto kk = groups.find(*jj);
            if (max_size < kk->second.size())
                max_size = kk->second.size();
        }
    }
    std::vector <uint32_t> group_sizes; group_sizes.resize(max_size, 0);
    for (auto ii = groups_sorted.rbegin(); ii != groups_sorted.rend(); ++ii)
    {
        if (ii->first < c) break;
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
        {
            auto kk = groups.find(*jj);
            if (kk != groups.end())
            {
                if (kk->second.size()-1 >= group_sizes.size())
                { 
                    std::cerr << "kk->second.size()-1 >= group_sizes.size() (" << kk->second.size()-1 << " > " << group_sizes.size() << ")" << std::endl;
                    return;
                }
                group_sizes.at(kk->second.size()-1)++;
            }
        }
    }

    // print groups: all and selected
    std::string outname(out_dir); outname.append("/groups_all.tsv");
    std::ofstream outfile_all; outfile_all.open(outname.c_str());
    if (not outfile_all.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outname = out_dir; outname.append("/groups.tsv");
    std::ofstream outfile; outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile_all << "group	size	hub	hub_count	group_count	out_hub_degree" << std::endl;
    outfile << "group	size	hub	hub_count	group_count	out_hub_degree" << std::endl;
    size_t group_id = 1;
    for (auto ii = groups_sorted.rbegin(); ii != groups_sorted.rend(); ++ii)
    {
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
        {
            std::ostringstream ss;
            ss << group_id++ << "	"
               << groups.find(*jj)->second.size() << "	"
               << H.getRead(*jj).id() << "	"
               << H.getRead(*jj).count() << "	" 
               << ii->first << "	"
               << (_graph.at(*jj).neighbors.size() - (groups.find(*jj)->second.size() - 1));
            outfile_all << ss.str() << std::endl;
            if (ii->first >= c) outfile << ss.str() << std::endl;
        }
    }
    outfile.close();

    // print node and parent information (parent = hub of the group, node = element of the group)
    outname = out_dir; outname.append("/groups_info.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "node	hub	dist	count" << std::endl;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        outfile << H.getRead(ii->node_idx).id() << "	" 
                << H.getRead(ii->parent).id() << "	"
                << H.ComputeHammingDistance(ii->node_idx, ii->parent) << "	"
                << H.getRead(ii->node_idx).count() << std::endl;
    }
    outfile.close();
    
    // collect the edges connecting nodes belonging to different groups
    outname = out_dir; outname.append("/out_of_groups.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "node1	node2	hub1	hub2" << std::endl;
    for (auto ii = _graph.begin(); ii != _graph.end(); ++ii)
    {
        for (auto jj = ii->neighbors.begin(); jj != ii->neighbors.end(); ++jj)
        {
            if (ii->node_idx < *jj) continue;
            if (ii->parent != _graph.at(*jj).parent) 
                outfile << H.getRead(ii->node_idx).id() << "	" 
                        << H.getRead(*jj).id() << "	"
                        << H.getRead(ii->parent).id() << "	"
                        << H.getRead(_graph.at(*jj).parent).id() << std::endl;
        }
    }
    outfile.close();
 
    // all groups size frequency
    outname = out_dir; outname.append("/group_all_size_freq.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "size	num_groups" << std::endl;
    size_t i = 1;
    for (auto ii = group_sizes_all.begin(); ii != group_sizes_all.end(); ++ii)
        outfile << i++ << "	" << *ii << std::endl;
    outfile.close();
 
    // selected groups size frequency
    outname = out_dir; outname.append("/group_size_freq.tsv");
    outfile.open(outname.c_str());
    if (not outfile.is_open()) { std::cerr << "Error while opening file " << outname << " for writing" << std::endl; return; }
    outfile << "size	num_groups" << std::endl;
    i = 1;
    for (auto ii = group_sizes.begin(); ii != group_sizes.end(); ++ii)
        outfile << i++ << "	" << *ii << std::endl;
    outfile.close();
     
}
    

/*
void                          
Graph::print_stars (const Hash & H, const std::string & out)
{
    std::ofstream outfile;
    outfile.open(out.c_str());
    if (not outfile.is_open())
    {
        std::cerr << "Error while opening file " << out << " for writing" << std::endl;
        return;
    }

    outfile << "group	size	hub	elements	dist_hub_elements" << std::endl;

    // collect the groups into sets labelled with the hub node 
    std::map < uint32_t, std::map <uint32_t, uint32_t> > groups = collect_stars ();

    size_t group_id = 0;
    for (auto ii = groups.begin(); ii != groups.end(); ++ii)
    {
        auto hub = H.getRead(ii->first).id();
        std::string el("");
        std::string d("");
        for (auto jj = ii->second.begin(); jj != ii->second.end(); ++jj)
        {
            if (el.length() > 0)
            {
                el += ',';
                d += ',';
            }
            el.append(H.getRead(jj->first).id());
            d.append(std::to_string(jj->second));
        }
        outfile << group_id << "	"
                << ii->second.size() << "	"
                << hub << "	"
                << el << "	" 
                << d << std::endl;
        ++group_id;
    }
    
}
*/




