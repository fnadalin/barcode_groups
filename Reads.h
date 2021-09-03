#ifndef READS_H_
#define READS_H_

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <stdint.h>

class Reads {

public:

    virtual ~Reads () {}
    Reads () { };
    Reads (const std::string & id, const std::string & s);
    Reads (const std::string & id, const std::string & s, uint32_t count);
    Reads (const Reads &r); // copy constructor

    Reads &               operator=(const Reads & r);

    uint32_t              length() const { return _length; } // Return the length of the sequence
    std::string           toString() const; // returns char representation of string
    std::string           id () const { return _id; }
    uint32_t              count () const { return _count; }

protected:

    std::vector<uint8_t>  _read; // read as given in input
    std::string           _id; // sequence identifier
    uint32_t              _length;
    uint32_t              _count; // number of events

};

#endif /* READS_H_ */
