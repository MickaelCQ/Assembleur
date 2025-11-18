//
// Created by raphael on 11/13/25.
//

#ifndef ASSEMBLEUR_COMPARE_H
#define ASSEMBLEUR_COMPARE_H
#include <string>
#include <vector>
#include <cstddef> // Pour size_t
#include "bitvector.h"

/**
 * @struct CompareKMers
 * @brief Compare les k-mers à partir d'un vecteur binaire de nucléotides.
 *
 * Cette classe stocke toutes les lectures d'un fichier FASTA sous forme
 * d'un unique vecteur binaire (où chaque nucléotide est codé sur 2 bits)
 * et un vecteur des positions de fin de chaque lecture (en bits).
 */
struct CompareKMers {
    BitVector bitVector; // Utilisation directe de BitVector au lieu de std::vector<bool>
    std::vector<size_t> reads;      // Vecteur des positions de fin de lecture (en bits)
    size_t kmersize = 31;

    CompareKMers(BitVector bitVector, std::vector<size_t> reads, size_t kmersize);
    CompareKMers(BitVector bitVector, std::vector<size_t> reads);

    // Setters
    void set_kmersize(size_t k);

    // Getters
    size_t get_kmersize() const;
    BitVector& get_bitVector();
    std::vector<size_t>& get_reads();
    size_t get_read_end_pos(size_t read_idx) const;
    size_t get_nKmers(size_t ref_read_idx) const;
    size_t get_all_nKmers() const;
    size_t get_nReads() const;

    // Methods
    /**
     * @brief Compare le chevauchement (k-1) de deux k-mers.
     * @param bit_idx1 Index binaire de début du premier k-mer.
     * @param bit_idx2 Index binaire de début du second k-mer.
     * @return 1 si le suffixe (k-1) de k-mer1 correspond au préfixe (k-1) de k-mer2, sinon 0.
     */
    size_t compare_line(size_t bit_idx1, size_t bit_idx2) const;

    /**
     * @brief Compare le premier k-mer d'une lecture de référence à tous les autres.
     * @param ref_read_idx Index de la lecture de référence.
     * @return Vecteur de résultats (1 pour match, 0 sinon).
     */
    std::vector<size_t> compare_lines(size_t ref_read_idx) const;

    /**
     * @brief Exécute compare_lines pour chaque lecture.
     * @return Matrice 2D des résultats de comparaison.
     */
    std::vector<std::vector<size_t>> compare_all() const;
};

#endif //ASSEMBLEUR_COMPARE_H
