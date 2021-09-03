/*
 ============================================================================
 Name        : barcode_groups.cpp
 Author      : Francesca Nadalin
 Version     : 
 Description : Compute the groups of sequences that share the same barcode
 ============================================================================
 */

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "Hash.h"
#include "Graph.h"

#define RESERVE_BLOCK 100

int main(int argc, char * argv[]) {

    if (argc < 5)
    {
        std::cout << std::endl;
        std::cout << "Usage: " << argv[0] << " <sequences> <mismatches> <frac_count> <outdir>" << std::endl;
        std::cout << std::endl;
        std::cout << "    <sequences>   ID and sequence and (optionally) read count, tab/space separated, one barcode per line" << std::endl;
        std::cout << "    <mismatches>  max number of mismatches between sequences in the same group" << std::endl;
        std::cout << "    <frac_count>  max fraction of read count for a sequence to be considered a neighbor of the hub" << std::endl;
        std::cout << "    <outdir>      output directory where results and stats are saved" << std::endl;
        std::cout << std::endl;
        return 0;
    }

    std::string seq_file(argv[1]);
    uint32_t mismatches = atoi(argv[2]);
    double frac_count = atof(argv[3]);
    std::string out_dir(argv[4]);

    // print params

    std::string params(out_dir);
    params.append("/params.txt");
    std::ofstream params_file;
    params_file.open(params.c_str());
    if (params_file.good())
    {
        params_file << "seq_file: " << seq_file << std::endl;
        params_file << "mismatches: " << mismatches << std::endl;
        params_file << "out_dir: " << out_dir << std::endl;
    }

    clock_t time1 = clock();

    // sequence indexing

    Hash H(mismatches);
    std::cout << "store reads" << std::endl;
    H.store_reads(seq_file);

    std::cout << "fill hash" << std::endl;
    H.fill_hash();

    // graph build and visit
  
    Graph G(mismatches, frac_count);
    std::cout << "build graph" << std::endl;
    G.build(H);
    G.print_graph_attributes(H, out_dir);

    std::cout << "create groups" << std::endl;
//  G.create_stars();
    G.create_stars_with_counts(H);
    G.print_stars(H, out_dir);
//    G.print_stars_with_counts(H, out_groups);
//    G.print_stars_info(H, out_groups_info);
//    G.print_out_of_stars(H, out_out_of_groups);
//    G.print_group_size_frequency(out_group_size_freq);
 
    double time2 = clock();    
    std::cout << std::endl;
    std::cout << "Wall time: " << (double) (time2 - time1) / CLOCKS_PER_SEC << " s" << std::endl;

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "RSS: " << usage.ru_maxrss << std::endl;

    return 0;
}

