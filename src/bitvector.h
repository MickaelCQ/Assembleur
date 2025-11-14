#pragma once

#include <vector>
#include <bitset>
#include <string>
#include <utility>
#include <cstddef>

class BitVector
{
protected:
    std::vector<std::bitset<60>> _bitVector;
    size_t _sizeElement = 4;            // length of the element in the bitset (not used for basic push_back)
    size_t _bitCount = 0;               // total number of bits stored

    std::vector<std::pair<char, int>> _listBit;   // list of bit codes associated to a char (optional)

public :
    /*________Constructor__________*/
    BitVector();
    BitVector(int sizeElement);
    BitVector(int sizeElement, size_t dataSize);

    /*_________Vector-like API (minimal)___________*/
    void reserve(size_t bits);
    void push_back(bool bit);
    size_t size() const;
    void clear();

    // Convert to a std::vector<bool>
    std::vector<bool> to_vector() const;

    /*_________Getter/setter___________*/

    /**
    * @brief Size of one element in bits.
    */
    void setSizeElement(int sizeElement);
    int getSizeElement() const;

    /**
     * @brief Create the dictionary of characters associated with a bits code.
     */
    void createListBit(const std::string& list);
    const std::vector<std::pair<char, int>>& getListBit() const;

    /**
    * @brief Add a character in bits to the bitvector (interprets nucleotide -> two bits).
    */
    void addCha(char cha);

    /**
    * @brief Convert back to sequence (naive: interpret every 2 bits as nucleotide using A,C,G,T mapping).
    */
    std::string readBitVector() const;
};
