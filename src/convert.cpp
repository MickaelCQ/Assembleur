//
// Created by raphael on 11/13/25.
//

#include "convert.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <algorithm> // For std::remove_if
#include <cctype>    // For ::isspace

void Convert::process_fasta_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // Clear any previous data
    bit_vector.clear();
    read_end_positions.clear();
    _cache_valid = false;

    size_t total_read_size = 0;
    size_t total_read_number = 0;
    std::string line;

    // --- FIRST PASS: Calculate total size and read count ---
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            total_read_number++; // Found a new read header
        } else {
            // This is a sequence line.
            // Remove any potential whitespace from the line before counting.
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            total_read_size += line.length();
        }
    }

    // --- MEMORY RESERVATION ---
    // Reserve space for all bits (2 bits per nucleotide)
    bit_vector.reserve(total_read_size * 2);
    // Reserve space for all read end positions
    read_end_positions.reserve(total_read_number);

    // --- SECOND PASS: Convert and store data ---
    // Reset file stream to the beginning
    file.clear();
    file.seekg(0, std::ios::beg);

    std::stringstream current_sequence_stream;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue; // Skip empty lines
        }

        if (line[0] == '>') {
            // This is a header line.
            // Process the sequence we've built up so far.
            std::string sequence = current_sequence_stream.str();
            if (!sequence.empty()) {
                convert_and_store_sequence(sequence);
            }

            // Reset the stream for the new sequence
            current_sequence_stream.str("");
            current_sequence_stream.clear();
        } else {
            // This is a sequence line. Append it to our current sequence stream.
            // We strip any potential whitespace from the line.
            line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
            current_sequence_stream << line;
        }
    }

    // After the loop, process the very last sequence in the file
    std::string last_sequence = current_sequence_stream.str();
    if (!last_sequence.empty()) {
        convert_and_store_sequence(last_sequence);
    }

    file.close();
}

void Convert::convert_and_store_sequence(const std::string& sequence) {
    if (sequence.empty()) {
        return;
    }

    for (char nucleotide : sequence) {
        bit_vector.addCha(nucleotide);
    }

    // Store the *new* cumulative size of the bit vector
    read_end_positions.push_back(bit_vector.size());
    _cache_valid = false;
}

const std::vector<bool>& Convert::get_bit_vector() const {
    if (!_cache_valid) {
        _bit_vec_cache = bit_vector.to_vector();
        _cache_valid = true;
    }
    return _bit_vec_cache;
}

const std::vector<size_t>& Convert::get_read_end_positions() const {
    return read_end_positions;
}