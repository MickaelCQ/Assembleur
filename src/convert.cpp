#include "convert.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cctype>

void Convert::processFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // Réinitialisation de l'état interne
    bitVector = BitVector();
    read_end_positions.clear();

    size_t total_read_size = 0;
    size_t total_read_number = 0;
    std::string line;
    
    std::string current_sequence;

    // --- RÉSERVATION MÉMOIRE ---
    // Réserve 2 bits par nucléotide.
    bitVector.reserve(total_read_size * 2);
    read_end_positions.reserve(total_read_number);
    
    // --- PREMIÈRE PASSE : Statistiques ---
    // Compte simplement pour optimiser l'allocation mémoire (critique pour les grands génomes).
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            total_read_number++;
          if (!current_sequence.empty()) {
            convertSeq(current_sequence);
            current_sequence.clear();
          }
        } else {
            // Supprime les espaces blancs pour obtenir le vrai compte de nucléotides
            line.erase(std::remove_if(line.begin(), line.end(), [](char c) {
                return std::isspace(static_cast<unsigned char>(c));
            }), line.end());
            total_read_size += line.length();
            current_sequence += line;
        }
    }

  // Traite la toute dernière lecture (car aucun '>' ne la déclenche après)
    if (!current_sequence.empty()) {
        convertSeq(current_sequence);
    }
  
    file.close();
}

void Convert::convertSeq(const std::string& sequence) {
    if (sequence.empty()) return;

    for (char nucleotide : sequence) {
        bitVector.addCha(nucleotide);
    }

    // Enregistre la frontière de cette lecture dans le flux de bits.
    read_end_positions.push_back(bitVector.size());
}

const BitVector& Convert::get_bitVector() const {
    return bitVector;
}

const std::vector<size_t>& Convert::get_read_end_positions() const {
    return read_end_positions;
}
