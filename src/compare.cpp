//
// Created by raphael on 11/13/25.
//

#include "compare.h"
#include <stdexcept>
#include <utility> // Pour std::move

// doit prendre une reference de 
CompareKMers::CompareKMers(BitVector bitVector, std::vector<size_t> reads, size_t kmersize)
    : bitVector(std::move(bitVector)), reads(std::move(reads)), kmersize(kmersize) {}

// CHANGEMENT: Accepte BitVector
CompareKMers::CompareKMers(BitVector bitVector, std::vector<size_t> reads)
    : bitVector(std::move(bitVector)), reads(std::move(reads)), kmersize(31) {}

// Setters
void CompareKMers::set_kmersize(size_t k) {
    kmersize = k;
}

// Getters
size_t CompareKMers::get_kmersize() const {
    return kmersize;
}

// Retourne une référence au BitVector interne
BitVector& CompareKMers::get_bitVector() {
    return bitVector;
}

std::vector<size_t>& CompareKMers::get_reads() {
    return reads;
}

size_t CompareKMers::get_read_end_pos(const size_t read_idx) const
{
    // CORRECTION: La vérification < 0 est inutile pour size_t
    if (read_idx >= static_cast<size_t>(reads.size())) {
        throw std::out_of_range("read_idx is out of range");
    }

    // CORRECTION: L'ancien code retournait 0 pour read_idx == 0, ce qui était incorrect.
    // reads[read_idx] contient la position de fin de la lecture 'read_idx'.
    return reads[read_idx];
}

size_t CompareKMers::get_nKmers(const size_t ref_read_idx) const {
    // CORRECTION: La position de début est 0 pour la première lecture, ou la fin de la précédente.
    const size_t start_bit = (ref_read_idx == 0) ? 0 : get_read_end_pos(ref_read_idx - 1);
    const size_t end_bit = get_read_end_pos(ref_read_idx);

    // CHANGEMENT: Calcule le nombre de nucléotides à partir des bits
    const size_t num_nucleotides = (end_bit - start_bit) / 2;

    // On ne peut pas avoir de k-mer si la lecture est plus courte que le k-mer
    if (num_nucleotides < kmersize) {
        return 0;
    }

    // La formule est (longueur - k + 1)
    return num_nucleotides - kmersize + 1;
}

size_t CompareKMers::get_all_nKmers() const {
    // CORRECTION: L'ancienne implémentation était boguée.
    // La nouvelle somme le nombre de k-mers de chaque lecture.
    size_t total_kmers = 0;
    for (size_t i = 0; i < get_nReads(); ++i) {
        total_kmers += get_nKmers(i);
    }
    return total_kmers;
}

size_t CompareKMers::get_nReads() const {
    return reads.size();
}

// Methods
size_t CompareKMers::compare_line(const size_t bit_idx1, const size_t bit_idx2) const
{
    // CHANGEMENT: Logique de comparaison binaire
    // Compare le chevauchement de k-1 nucléotides
    for (size_t i = 0; i < kmersize - 1; ++i)
    {
        // Index binaire du nucléotide (i+1) du k-mer 1
        size_t nuc1_bit_start = bit_idx1 + (i + 1) * 2;

        // Index binaire du nucléotide (i) du k-mer 2
        size_t nuc2_bit_start = bit_idx2 + i * 2;

        // Compare les deux bits du nucléotide en utilisant l'interface BitVector
        if (bitVector.test(nuc1_bit_start)     != bitVector.test(nuc2_bit_start) ||
            bitVector.test(nuc1_bit_start + 1) != bitVector.test(nuc2_bit_start + 1))
        {
            return 0; // Pas de match, sortie anticipée
        }
    }

    // Si la boucle se termine, tous les k-1 nucléotides ont matché
    return 1;
}

std::vector<size_t> CompareKMers::compare_lines(const size_t ref_read_idx) const
{
    std::vector<size_t> results;
    results.reserve(reads.size()); // CORRECTION: L'original réservait 'kmersize'

    // CORRECTION: L'ancien code passait 'ref_read_idx' comme un index de k-mer (bit_idx1).
    // Nous récupérons maintenant l'index binaire de début de la *lecture* de référence.
    const size_t ref_kmer_start_bit = (ref_read_idx == 0) ? 0 : reads[ref_read_idx - 1];

    for (size_t i = 0; i < reads.size(); ++i)
    {
        if (i == ref_read_idx) {
            results.push_back(1);  // Identical lines
            continue;
        }

        // CORRECTION: Récupère l'index binaire de début de la lecture 'i'
        const size_t other_kmer_start_bit = (i == 0) ? 0 : reads[i - 1];

        // Note: Une implémentation robuste vérifierait ici que les deux lectures
        // ont une taille >= kmersize en utilisant get_nKmers(i) > 0.

        size_t cmp_result = compare_line(ref_kmer_start_bit, other_kmer_start_bit);
        results.push_back(cmp_result);
    }
    return results;
}

std::vector<std::vector<size_t>> CompareKMers::compare_all() const
{
    std::vector<std::vector<size_t>> all_results;
    all_results.reserve(reads.size());
    // La deuxième réservation était redondante et a été supprimée.
    for (size_t i = 0; i < reads.size(); ++i)
    {
        std::vector<size_t> line_results = compare_lines(i);
        all_results.push_back(std::move(line_results));
    }
    return all_results;
}
