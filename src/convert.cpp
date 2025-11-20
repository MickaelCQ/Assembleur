#include "convert.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cctype>

using namespace std;

void Convert::processFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) 
    {
        throw runtime_error("Error: Could not open file " + filename);
    }

    // Réinitialisation de l'état interne
    bitVector = BitVector();
    endPos.clear();

    size_t totReadSize = 0;
    size_t totReadNum = 0;
    string line;
    
    string cSeq;

    // --- RÉSERVATION MÉMOIRE ---
    // Réserve 2 bits par nucléotide.
    bitVector.reserve(totReadSize * 2);
    endPos.reserve(totReadNum);
    
    // --- Statistiques ---
    // Compte simplement pour optimiser l'allocation mémoire (critique pour les grands génomes).
    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') 
        {
            totReadNum++;
          if (!cSeq.empty()) 
          {
            convertSeq(cSeq);
            cSeq.clear();
          }
        } 
        else 
        {
            // Supprime les espaces blancs pour obtenir le vrai compte de nucléotides
            line.erase(remove_if(line.begin(), line.end(), [](char c) {
                return isspace(static_cast<unsigned char>(c));
            }), line.end());
            totReadSize += line.length();
            cSeq += line;
        }
    }

  // Traite la toute dernière lecture (car aucun '>' ne la déclenche après)
    if (!cSeq.empty()) {
        convertSeq(cSeq);
    }
  
    file.close();
}

void Convert::convertSeq(const string& sequence) {
    if (sequence.empty()) return;

    for (char nucleotide : sequence) {
        bitVector.addCha(nucleotide);
    }

    // Enregistre la frontière de cette lecture dans le flux de bits.
    endPos.push_back(bitVector.size());
}

const BitVector& Convert::getBitVector() const {
    return bitVector;
}

const vector<size_t>& Convert::getEndPos() const {
    return endPos;
}
