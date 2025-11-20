#include "bitvector.h"
#include <stdexcept>
#include <cctype>
#include <iostream>

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
        case 'A': push_back(false); push_back(false); break;  // 00
        case 'C': push_back(true); push_back(false); break;  // 10
        case 'G': push_back(false); push_back(true); break;  // 01
        case 'T': push_back(true); push_back(true); break;  // 11
        default: throw std::invalid_argument("Invalid nucleotide character");
    }
}

std::string BitVector::readBitVector() const {
    std::string out;
    out.reserve(_bitCount / 2 + 1);
    for (size_t i = 0; i + 1 < _bitCount; i += 2) {
        bool b1 = test(i);
        bool b2 = test(i+1);

        if (!b1 && !b2) out.push_back('A');
        else if (b1 && !b2) out.push_back('C');
        else if (!b1 && b2) out.push_back('G');
        else out.push_back('T');
    }
    return out;
}

bool BitVector::test(size_t idx) const {
    if (idx >= _bitCount) {
        return false; // Fail-safe
    }
    size_t block = idx / BLOCK_BITS;
    size_t offset = idx % BLOCK_BITS;
    return _bitVector[block].test(offset);
}

void BitVector::append(const BitVector& other, size_t skip_nucleotides) {
    size_t start_bit = skip_nucleotides * 2;
    if (start_bit >= other.size()) return;

    for (size_t i = start_bit; i < other.size(); ++i) {
        push_back(other.test(i));
    }
}

uint8_t BitVector::getNucleotideAt(size_t index) const {
    size_t bit_idx = index * 2;
    bool b1 = test(bit_idx);
    bool b2 = test(bit_idx + 1);

    // Reconstitution de la valeur 2 bits
    uint8_t val = 0;
    if (b1) val |= 2; // 10
    if (b2) val |= 1; // 01
    return val;
}

void BitVector::resize(size_t num_bits) {
    if (num_bits >= _bitCount) return; // We only support shrinking here for safety
    _bitCount = num_bits;
    // We don't necessarily need to shrink the vector memory, just the count
}
