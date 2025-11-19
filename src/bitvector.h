#pragma once

#include <vector>
#include <bitset>
#include <string>
#include <utility>
#include <cstddef>
#include <cstdint>

/**
 * @class BitVector
 * @brief Conteneur spécialisé pour stocker des séquences d'ADN de manière efficace via compression binaire.
 *
 * Au lieu de stocker 1 octet (8 bits) par nucléotide (char), cette classe stocke
 * 2 bits par nucléotide (A=00, C=10, G=01, T=11). Cela réduit l'utilisation mémoire par 4.
 * Les données sont stockées par blocs de BLOCK_BITS bits.
 */
class BitVector
{
public:
    // Constante : nombre de bits par bloc (doit être une valeur compile-time).
    static constexpr size_t BLOCK_BITS = 64;

protected:
    // Conteneur de stockage : vecteur de blocs de BLOCK_BITS bits.
    std::vector<std::bitset<BLOCK_BITS>> _bitVector;

    size_t _sizeElement = 2;            // Bits par élément (2 pour les nucléotides)
    size_t _bitCount = 0;               // Nombre total de bits actuellement utilisés
    std::vector<std::pair<char, int>> _listBit;   // Dictionnaire pour l'encodage (optionnel)

public :
    /*________Constructeurs__________*/
    BitVector();

    /**
     * @brief Constructeur définissant la taille de l'élément (défaut 2 pour l'ADN).
     */
    BitVector(int sizeElement);

    /**
     * @brief Constructeur pré-allouant la mémoire.
     */
    BitVector(int sizeElement, size_t dataSize);

    /*_________API type Vector (minimaliste)___________*/

    /**
     * @brief Pré-alloue les blocs mémoire pour éviter les réallocations lors du remplissage.
     * @param bits Nombre total de bits attendus.
     */
    void reserve(size_t bits);

    /**
     * @brief Ajoute un seul bit à la fin du vecteur.
     * Gère la logique de saut au bloc suivant quand le bloc actuel est plein.
     */
    void push_back(bool bit);

    /**
     * @brief Retourne le nombre total de bits stockés.
     */
    size_t size() const;

    /**
     * @brief Réinitialise le vecteur (vide le contenu et remet le compteur à 0).
     */
    void clear();

    /**
     * @brief Convertit les bitsets internes en un vecteur booléen standard (pour débogage).
     * ne pas les prendre en compte dans l'algorithme micka
     */
    std::vector<bool> to_vector() const;

    /*_________Accesseurs / Mutateurs___________*/

    void setSizeElement(int sizeElement);
    int getSizeElement() const;

    void createListBit(const std::string& list);
    const std::vector<std::pair<char, int>>& getListBit() const;

    /**
     * @brief Encode un caractère nucléotidique en 2 bits et les ajoute.
     * Mapping : A->00, C->10, G->01, T->11.
     * @param cha Le caractère (A, C, G, T).
     */
    void addCha(char cha);

    /**
     * @brief Décode le vecteur binaire pour retrouver la chaîne ADN originale.
     * @return Représentation string de la séquence.
     */
    std::string readBitVector() const;

    // --- Accès aux bits ---
    /**
     * @brief Vérifie la valeur d'un bit spécifique à un index global.
     * Calcule quel bloc et quel décalage (offset) vérifier.
     * @param idx Index global du bit.
     * @return true (1) ou false (0).
     */
    bool test(size_t idx) const;

    bool operator[](size_t idx) const { return test(idx); }
};