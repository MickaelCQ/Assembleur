#ifndef ASSEMBLEUR_CONVERT_H
#define ASSEMBLEUR_CONVERT_H

#include <string>
#include <vector>
#include <cstdint>
#include "bitvector.h"

/**
 * @class Convert
 * @brief Traite les fichiers FASTA, convertissant les séquences ADN en un vecteur binaire compact.
 *
 * Stratégie :
 * 1. **Passe 1** : Scanne le fichier pour compter les nucléotides. Permet une réservation mémoire exacte ($O(1)$ réallocations).
 * 2. **Passe 2** : Lit les séquences, nettoie les retours à la ligne/espaces, et remplit le `BitVector`.
 *
 * Elle maintient un vecteur `read_end_positions` pour savoir où une lecture s'arrête et où la suivante commence
 * au sein du flux continu de bits.
 */
class Convert
{
private:
    BitVector bitVector;
    // Stocke l'index cumulatif (en bits) où chaque lecture se termine.
    // Exemple : Lecture 1 fait 10pb (20 bits). Lecture 2 fait 10pb.
    // Vecteur : [20, 40].
    std::vector<size_t> read_end_positions;

    /**
     * @brief Convertit chaîne et l'ajoute au bitVector.
     * @complexity : space : O(b), time : O(s).
     * b = nb de block (pire cas, sinon O(1) en espace). s = taille sequence.
     */
    void convertSeq(const std::string& sequence);

public:
    Convert() = default;

    /**
     * @brief Pilote principal : Ouvre le FASTA, exécute la logique en 2 passes, remplit les structures.
     * @complexity : space : O(n/64 + r + s), time : O(n + h).
     * n = taille totale des séquences nucleotide. h = taille header.
     * r = read_end_positions.size(). s = sequence.size().
     */
    void processFile(const std::string& filename);

    /**
     * @brief get BitVector.
     * @complexity : space : O(1), time : O(1).
     */
    const BitVector& get_bitVector() const;

    /**
     * @brief get read end positions.
     * @complexity : space : O(1), time : O(1).
     */
    const std::vector<size_t>& get_read_end_positions() const;
};

#endif //ASSEMBLEUR_CONVERT_H
