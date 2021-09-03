
#include "Hash.h"

Hash::Hash (uint32_t mismatches)
{
    _mismatches = mismatches;
    _readLength = 0;
    _numReads = 0;

}

void
Hash::store_reads (const std::string & filename)
{
    std::ifstream infile;
    infile.open(filename.c_str());
    if (not infile.is_open())
    {
        std::cerr << "Error while opening file " << filename << " for reading" << std::endl;
	return;
    }

    bool init = true;
    auto line_num = 1;
    while (infile.good())
    {
        std::string line;
        getline (infile, line);
        if (line.size() > 0)
        {
            std::stringstream ss(line);
            std::string id("");
            std::string s("");
            char c;
            uint32_t count = 0;
            while (std::isalnum(ss.peek()))
            {
                ss.get(c);
		id += c;
            }
            if (not std::isblank(ss.peek())) continue; // skip incorrectly formatted lines
            while (std::isblank(ss.peek()))
                ss.get(c);
            while (std::isalpha(ss.peek()))
            {
                ss.get(c);
                s += c;
            }
            if (!init and s.length() != _readLength)
            {
                std::cerr << "Sequences must all have the same length (line: " << line_num << ")" << std::endl;
                return;
            }
            while (std::isblank(ss.peek()))
                ss.get(c);
            while (std::isdigit(ss.peek()))
            {
                ss.get(c);
                count = count*10 + c - 48;
            }
            Reads r(id,s,count);
            _readsMulti.push_back(r);
            _readLength = r.length();
            init = false;
            ++line_num;
        }
    }
 
    _numReads += _readsMulti.size();

    _numBlocks = _mismatches + 1;
    _blockSize = (_readLength + _numBlocks - 1) / _numBlocks;

    _numFingerprints = (1 << 2*_blockSize);

    std::cout << "num sequences: " << _numReads << std::endl;
    std::cout << "block size: " << _blockSize << std::endl;
    std::cout << "num fingerprints: " << _numFingerprints << std::endl;

}

void
Hash::fill_hash ()
{
    _HASHcounter.assign(_numFingerprints, 0);
    _HASHvalues.resize(_numReads*_numBlocks+1);
    unsigned long int fp;
    // compute # reads with fingerprint = fp
    for (auto i = 0; i < _numReads; ++i) 
    {
        for (auto id = 0; id < _numBlocks; ++id)
        {
            fp = ComputeFingerprint(i, id); // compute the fingerprint for the i-th block
            _HASHcounter.at(fp)++;
        }
    }
    // compute # reads with fingerprint < fp
    uint32_t t1 = _HASHcounter.at(0);
    _HASHcounter.at(0) = 0;
    for (size_t i = 1; i < _numFingerprints; ++i) {
        int t2 = _HASHcounter.at(i);
        _HASHcounter.at(i) = _HASHcounter.at(i-1) + t1;
        t1 = t2;
    }

    // assign read pointer and update # reads with fingerprint <= fp
    for (size_t i = 0; i < _numReads; i++) 
    {
        for (size_t id = 0; id < _numBlocks; ++id)
        {
            fp = ComputeFingerprint(i,id);
            fragment_type f;
            f.read_id = i; f.block_id = id;
            _HASHvalues.at(_HASHcounter.at(fp)) = f;
            _HASHcounter.at(fp)++;
        }
    }

    // re-store # of reads with fingerprint < fp (so that it corresponds to the first position in _HASHvalues)
    for (uint32_t i = _numFingerprints - 1; i > 0 ; --i) 
        _HASHcounter.at(i) = _HASHcounter.at(i-1);
    _HASHcounter[0] = 0;

}

// TODO: this can be optimized, without converting to string
uint32_t                
Hash::ComputeFingerprint (uint32_t i, size_t f_id) const
{
    std::string read = _readsMulti.at(i).toString();
    uint32_t fingerprint = 0;

    auto start = f_id*_blockSize;
    auto end = ( (f_id+1)*_blockSize < read.length() ) ? (f_id+1)*_blockSize : read.length();
    for (size_t j = start; j < end; ++j) 
    {
        char c = 0;
        switch(read.at(j)) {
        case 'a' : case 'A' : c = 0; break;
        case 'c' : case 'C' : c = 1; break;
        case 'g' : case 'G' : c = 2; break;
        case 't' : case 'T' : c = 3; break;
        }
        fingerprint = (fingerprint << 2) + c;
    }

    return fingerprint;
}

// TODO: this can be optimized, without converting to string
uint32_t                
Hash::ComputeHammingDistance (uint32_t i, uint32_t j) const
{
    if (i >= _readsMulti.size()) std::cerr << "i >= _readsMulti.size() [" << i << "," << _readsMulti.size() << "]" << std::endl;
    std::string read1 = _readsMulti.at(i).toString();
    std::string read2 = _readsMulti.at(j).toString();

    uint32_t d = 0;
    for (size_t k = 0; k < _readLength; ++k) 
    {
        d += (read1.at(k) != read2.at(k));
    }

//    std::cout << "r1: " << read1 << std::endl;
//    std::cout << "r2: " << read2 << std::endl;
//    std::cout << "distance = " << d << std::endl << std::endl;

    return d;
}

