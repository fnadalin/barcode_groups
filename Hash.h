#ifndef HASH_H_
#define HASH_H_

#include <iostream>
#include <fstream>
#include <vector>

#include "Reads.h"

typedef struct {
    uint32_t read_id; // 0-based position in _readsMulti
    uint32_t block_id; // 0-based ID of the block (from 0 to _numBlocks-1)
} fragment_type;

class Graph;

class Hash {

public:

    static const uint32_t INFTY_DIST = 1000;

    Hash (uint32_t mismatches);
    ~Hash() { }

    void                        store_reads (const std::string & filename);

    // compute _numFingerprints on non-overlapping read substrings
    void                        fill_hash ();

    const Reads &               getRead (uint32_t i) const { return _readsMulti.at(i); }

    uint32_t                    ComputeFingerprint (uint32_t i, size_t fingerprint_id) const;

    unsigned int                numReads () const { return _numReads; }
    uint32_t                    numBlocks () const { return _numBlocks; }

    uint32_t                    ComputeHammingDistance (uint32_t i, uint32_t j) const;

    uint32_t                    size () const { return _numReads; }

    friend Graph;

protected:

    std::vector <Reads>         _readsMulti;
    std::vector <uint32_t>      _HASHcounter;
    std::vector <fragment_type> _HASHvalues;
  
    unsigned int                _numReads;
    unsigned int                _readLength;
    uint32_t                    _numBlocks;
    uint32_t                    _blockSize;
    uint32_t                    _mismatches;
    uint32_t                    _numFingerprints;

};


#endif /* HASH_H_ */
