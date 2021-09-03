
#include "Reads.h"

Reads::Reads(const std::string & id, const std::string & s)
{
    _id = id;
    _length = s.length();

    unsigned short int bufferREAD_length = (_length + 3) / 4;
    _read.assign (bufferREAD_length, 0); // init with zeroes
  
    int actualByte = 0;
    int actualBit = 0;
    for (unsigned int i = 0; i < s.length() ; i++) {
        switch (s.at(i)) {
        case 'A' : case 'a' : _read.at(actualByte) &= ~(1 << (7-actualBit)); _read.at(actualByte) &= ~(1 << (7-(actualBit+1))); break;
        case 'C' : case 'c' : _read.at(actualByte) &= ~(1 << (7-actualBit)); _read.at(actualByte) |= (1 << (7-(actualBit+1)));  break;
        case 'G' : case 'g' : _read.at(actualByte) |= (1 << (7-actualBit));  _read.at(actualByte) &= ~(1 << (7-(actualBit+1))); break;
        case 'T' : case 't' : _read.at(actualByte) |= (1 << (7-actualBit));  _read.at(actualByte) |= (1 << (7-(actualBit+1)));  break;
        default : _read.at(actualByte) &= ~(1 << (7-actualBit)); _read.at(actualByte) &= ~(1 << (7-(actualBit+1))); break;
        }
        actualBit += 2;
        if(actualBit == 8) {
            actualByte++;
            actualBit=0;
        }
    }

}

Reads::Reads(const std::string & id, const std::string & s, uint32_t count)
{
    Reads r(id, s);
    _read = r._read;
    _id = r._id;
    _length = r._length;
    _count = count;
}

Reads::Reads(const Reads &r)
{
    _read = r._read;
    _id = r._id;
    _length = r._length;
    _count = r._count;
}

Reads &
Reads::operator=(const Reads & r)
{
    if (this != &r) { // protect against invalid self-assignment
        this->_read = r._read;
        this->_id = r._id;
        this->_length = r._length;
        this->_count = r._count;
    }
    return *this;
}

std::string
Reads::toString () const
{
    uint8_t Mask = 0;
    uint32_t block = 0;
    Mask |= 1 << 7;
    Mask |= 1 << 6;

    std::string out;
    out.reserve(_length);

    uint8_t R = _read.at(0);

    for (auto i = 0; i < length(); ++i) 
    {
        if(i%4 == 0 and i > 0) 
        {
            ++block;
            R = _read.at(block);
        }
        uint8_t c = (R & Mask) >> 6;
        R = R << 2;
        switch(c) {
            case 0: out.append("A"); break;
            case 1: out.append("C"); break;
            case 2: out.append("G"); break;
            case 3: out.append("T"); break;
            default: std::cout << "strange" << std::endl;
        }
    }

    return out;
}


