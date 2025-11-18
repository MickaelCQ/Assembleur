//
// Created by raphael on 11/13/25.
//

#ifndef ASSEMBLEUR_CONVERT_H
#define ASSEMBLEUR_CONVERT_H

#include <string>
#include <vector>
#include <cstdint>      // For uint8_t
#include "bitvector.h"

/**
 * @class Convert
 * @brief Processes FASTA files, converting DNA sequences into a compact bit vector.
 *
 * This class reads DNA sequences from a FASTA file in two passes.
 * 1. The first pass calculates the total number of reads and total
 * nucleotide count to reserve memory efficiently.
 * 2. The second pass converts each nucleotide into a 2-bit representation
 * (A=00, C=10, G=01, T=11), and stores the result in a single,
 * continuous BitVector.
 *
 * It also tracks the cumulative end position (in bits) of each read.
 */
class Convert
{
private:
    BitVector bitVector;
    std::vector<size_t> read_end_positions;  // Ending position of a read in bits

    /**
     * @brief Converts a single DNA sequence string and appends it to the bit vector.
     * @param sequence The DNA sequence string (e.g., "ACGT").
     */
    void convertSeq(const std::string& sequence);

public:
    /**
     * @brief Default constructor.
     */
    Convert() = default;

    /**
     * @brief Opens, reads, and processes a FASTA file using a two-pass method.
     * @param filename The path to the FASTA file.
     * @throws std::runtime_error if the file cannot be opened.
     */
    void processFile(const std::string& filename);

    /**
     * @brief Gets the BitVector containing all reads.
     * @return A const reference to the `BitVector`.
     */
    const BitVector& get_bitVector() const;

    /**
     * @brief Gets the vector of read end positions.
     * @return A const reference to the `std::vector<size_t>`.
     * Each element `i` stores the cumulative bit position
     * at the end of read `i`.
     */
    const std::vector<size_t>& get_read_end_positions() const;
};


#endif //ASSEMBLEUR_CONVERT_H
