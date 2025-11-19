#include "bitvector.h"
#include <stdexcept>
#include <cctype>

BitVector::BitVector() = default;

BitVector::BitVector(int sizeElement) {
    setSizeElement(sizeElement);
}

BitVector::BitVector(int sizeElement, size_t dataSize) {
    setSizeElement(sizeElement);
    reserve(dataSize);
}

void BitVector::reserve(size_t bits) {
    if (bits == 0) return;
    size_t blocks = (bits + BLOCK_BITS - 1) / BLOCK_BITS;
    _bitVector.reserve(blocks);
}

void BitVector::push_back(bool bit) {
    size_t idx = _bitCount / BLOCK_BITS;
    size_t offset = _bitCount % BLOCK_BITS;

    if (idx >= _bitVector.size()) {
        _bitVector.emplace_back();
    }

    if (bit) {
        _bitVector[idx].set(offset);
    }
    _bitCount++;
}

size_t BitVector::size() const {
    return _bitCount;
}

void BitVector::clear() {
    _bitVector.clear();
    _bitCount = 0;
    _listBit.clear();
}

std::vector<bool> BitVector::to_vector() const {
    std::vector<bool> out;
    out.reserve(_bitCount);
    for (size_t i = 0; i < _bitCount; ++i) {
        size_t idx = i / BLOCK_BITS;
        size_t offset = i % BLOCK_BITS;
        out.push_back(_bitVector[idx].test(offset));
    }
    return out;
}

void BitVector::setSizeElement(int sizeElement) {
    if (sizeElement <= 0 || sizeElement > static_cast<int>(BLOCK_BITS)) {
        throw std::invalid_argument("sizeElement must be between 1 and BLOCK_BITS");
    }
    _sizeElement = static_cast<size_t>(sizeElement);
}

int BitVector::getSizeElement() const {
    return static_cast<int>(_sizeElement);
}

void BitVector::createListBit(const std::string& list) {
    _listBit.clear();
    int code = 0;
    for (char c : list) {
        _listBit.emplace_back(c, code);
        ++code;
    }
}

const std::vector<std::pair<char, int>>& BitVector::getListBit() const {
    return _listBit;
}

void BitVector::addCha(char cha) {
    switch (toupper(static_cast<unsigned char>(cha))) {
        case 'A':
            push_back(false); push_back(false); // 00
            break;
        case 'C':
            push_back(true); push_back(false);  // 10
            break;
        case 'G':
            push_back(false); push_back(true);  // 01
            break;
        case 'T':
            push_back(true); push_back(true);   // 11
            break;
        default:
            throw std::invalid_argument("Invalid nucleotide character: " + std::string(1, cha));
    }
}

std::string BitVector::readBitVector() const {
    std::string out;
    out.reserve(_bitCount / 2 + 1);
    for (size_t i = 0; i + 1 < _bitCount; i += 2) {
        bool b1 = _bitVector[i / BLOCK_BITS].test(i % BLOCK_BITS);
        bool b2 = _bitVector[(i+1) / BLOCK_BITS].test((i+1) % BLOCK_BITS);

        if (!b1 && !b2) out.push_back('A');
        else if (b1 && !b2) out.push_back('C');
        else if (!b1 && b2) out.push_back('G');
        else out.push_back('T');
    }
    return out;
}

bool BitVector::test(size_t idx) const {
    if (idx >= _bitCount) {
        throw std::out_of_range("BitVector::test: index out of range");
    }
    size_t block = idx / BLOCK_BITS;
    size_t offset = idx % BLOCK_BITS;
    return _bitVector[block].test(offset);
}
