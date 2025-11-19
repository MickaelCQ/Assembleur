#ifndef ASSEMBLEUR_GRAPHDBJ_H
#define ASSEMBLEUR_GRAPHDBJ_H

#include <vector>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <algorithm>
#include "convert.h"

/**
 * @struct Noeud
 * @brief Représente un (k-1)-mer dans le Graphe de De Bruijn.
 * * Dans un graphe de De Bruijn construit à partir de k-mers :
 * - Les nœuds représentent des séquences de longueur $k-1$.
 * - Les arêtes représentent la transition vers le $k$-ième nucléotide.
 */
struct Noeud {
    uint64_t p;             // La séquence du (k-1)-mer encodée en entier 64 bits.
    std::vector<Noeud*> c;  // Enfants (Arêtes sortantes / Out-edges).

    // --- Métadonnées pour le parcours de graphe ---
    std::vector<Noeud*> parents; // Parents (Arêtes entrantes / In-edges), crucial pour le parcours inverse (élagage).
    bool removed = false;        // Flag de suppression logique (évite de supprimer physiquement de la mémoire pendant la simplification).
    uint32_t coverage = 0;       // Combien de fois ce (k-1)-mer est apparu dans les lectures brutes (couverture).

    Noeud(uint64_t val) : p(val) {}
};

/**
 * @class GraphDBJ
 * @brief Construit, simplifie et parcourt le Graphe de De Bruijn.
 */
class GraphDBJ {
private:
    int k; // La taille du k-mer. Note : Les nœuds ont une taille k-1.

    // Hash map assurant que chaque (k-1)-mer est unique dans le graphe.
    // Clé : représentation uint64_t de la séquence. Valeur : Pointeur vers le Noeud.
    std::unordered_map<uint64_t, Noeud*> nodes_map;

    /**
     * @brief Calcule le complément inverse d'un k-mer encodé.
     * Exemple (si k=2): AC (00 10) -> GT (01 11).
     */
    uint64_t getReverseComplement(uint64_t val, int k) const;

    // ... suite de vos fonctions existantes ...
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;

    // --- Méthodes utilitaires ---

    /**
     * @brief Décode un uint64_t en string pour l'exportation.
     */
    std::string kmerToString(uint64_t val, int length) const;

    /**
     * @brief Supprime le lien entre parent et enfant dans les deux listes d'adjacence.
     */
    void disconnectNodes(Noeud* parent, Noeud* child);

    /**
     * @brief Heuristique pour trouver où deux branches divergentes se rejoignent (pour éclater les bulles).
     */
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);

public:
    /**
     * @brief Constructeur. Parse l'objet Convert et construit le graphe immédiatement.
     */
    GraphDBJ(const Convert& converter, int kmer_size);

    /**
     * @brief Destructeur. Nettoie tous les Noeuds alloués dynamiquement.
     */
    ~GraphDBJ();

    std::vector<Noeud*> getNodes() const;

    // --- Algorithmes d'Assemblage ---

    /**
     * @brief Élagage des pointes (Tip Clipping) : Supprime les courtes branches mortes (souvent erreurs de séquençage en fin de lecture).
     */
    int clipTips();

    /**
     * @brief Résolution des bulles (Bubble Collapsing) : Fusionne les chemins similaires causés par des SNP ou erreurs.
     */
    int resolveBubbles();

    /**
     * @brief Génération de Contigs : Parcourt le graphe simplifié pour sortir les séquences assemblées.
     * Utilise des heuristiques de couverture pour résoudre les ambiguïtés.
     */
    std::vector<std::string> generateContigs() const;

    // Ajoutez cette déclaration dans la section public de GraphDBJ
    void exportToGFA(const std::string& filename) const;

    /**
     * @brief Fusionne les contigs qui se chevauchent et supprime les doublons/inclusions.
     * @param contigs Liste brute des contigs.
     * @param min_overlap Taille minimale du chevauchement pour fusionner (ex: k).
     */
    static std::vector<std::string> mergeContigs(std::vector<std::string> contigs, int min_overlap);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H