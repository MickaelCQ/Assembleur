#ifndef ASSEMBLEUR_GRAPHDBJ_H
#define ASSEMBLEUR_GRAPHDBJ_H

#include <vector>
#include <cstdint> // Pour uint64_t
#include <unordered_map>
#include <string>       // NOUVEAU : Nécessaire pour std::string
#include <algorithm>    // NOUVEAU : Nécessaire pour std::find
#include "convert.h"

/**
 * @struct Noeud
 * @brief Représente un (k-1)-mer dans le graphe de De Bruijn.
 */
struct Noeud {
    uint64_t p;             // L'entier représentant la séquence du (k-1)-mer
    std::vector<Noeud*> c;  // Vecteur de pointeurs vers les enfants (arêtes sortantes)

    // --- NOUVEAUX CHAMPS ---
    std::vector<Noeud*> parents; // Parents (arêtes entrantes) pour remonter le graphe
    bool removed = false;        // Flag pour suppression logique (Soft delete)
    uint32_t coverage = 0;       // Compteur de couverture (nombre de fois vu)

    // Constructeur pour faciliter la création
    Noeud(uint64_t val) : p(val) {}
};

/**
 * @class GraphDBJ
 * @brief Construit et simplifie le graphe de De Bruijn.
 */
class GraphDBJ {
private:
    int k; // Taille du k-mer (les noeuds seront de taille k-1)

    // Hash map pour stocker l'unicité des noeuds
    std::unordered_map<uint64_t, Noeud*> nodes_map;

    /**
     * @brief Extrait une valeur entière sur 64 bits.
     */
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;

    // --- NOUVELLES MÉTHODES PRIVÉES (HELPERS) ---

    /**
     * @brief Convertit la valeur uint64_t en chaîne "ACGT".
     */
    std::string kmerToString(uint64_t val, int length) const;

    /**
     * @brief Supprime proprement les liens parent<->enfant dans les deux sens.
     */
    void disconnectNodes(Noeud* parent, Noeud* child);

    /**
     * @brief Cherche un point de convergence entre deux branches (pour les bulles).
     * C'est la déclaration qui manquait pour l'erreur ligne 173.
     */
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);

public:
    /**
     * @brief Constructeur. Construit immédiatement le graphe.
     */
    GraphDBJ(const Convert& converter, int kmer_size);

    /**
     * @brief Destructeur pour nettoyer la mémoire.
     */
    ~GraphDBJ();

    /**
     * @brief Récupère tous les noeuds du graphe.
     */
    std::vector<Noeud*> getNodes() const;

    // --- NOUVEAUX ALGORITHMES PUBLICS ---

    /**
     * @brief Élagage des pointes (Tip Clipping).
     */
    void removeTips(int length_threshold);

    /**
     * @brief Simplification des bulles (Bubble Collapsing).
     */
    void resolveBubbles();

    /**
     * @brief Génération des Contigs (extension gloutonne avec couverture).
     * C'est la déclaration qui manquait pour l'erreur ligne 217.
     */
    std::vector<std::string> generateContigs() const;
};

#endif //ASSEMBLEUR_GRAPHDBJ_H