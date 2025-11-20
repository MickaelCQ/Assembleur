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
 * @brief Représente un k-1 mer dans le graphe de De Bruijn.
 *
 * Chaque noeud stocke sa valeur canonique (sous forme d'entier 64 bits),
 * ses connexions (enfants/parents) et sa couverture (nombre d'occurrences).
 */
struct Noeud {
    uint64_t p;                 ///< Valeur encodée du k-1 mer (2 bits par nucléotide).
    std::vector<Noeud*> c;      ///< Liste des enfants (arcs sortants vers d'autres k-1 mers).
    std::vector<Noeud*> parents;///< Liste des parents (arcs entrants).
    bool removed = false;       ///< Flag pour marquer le noeud comme supprimé (lazy deletion).
    uint32_t coverage = 0;      ///< Profondeur de couverture (nombre de fois que ce k-1 mer a été vu).

    /**
     * @brief Constructeur simple.
     * @param val La valeur uint64 du k-mer.
     */
    Noeud(uint64_t val) : p(val) {}
};

/**
 * @struct GraphDBJConfig
 * @brief Conteneur pour tous les paramètres heuristiques de l'assembleur.
 *
 * Regroupe les seuils d'erreurs pour la fusion et les ratios pour le nettoyage du graphe
 * (suppression des pointes et résolution des bulles).
 */
struct GraphDBJConfig {
    // --- Paramètres de Fusion (Post-traitement) ---

    /** @brief Pourcentage d'erreur autorisé lors du chevauchement de contigs (ex: 0.05 = 5%). */
    double ERROR_PERCENT_OVERLAP = 0.05;

    /** @brief Pourcentage d'erreur autorisé pour considérer qu'un contig est inclus dans un autre. */
    double ERROR_PERCENT_CONTAINED = 0.02;

    // --- Paramètres de Nettoyage du Graphe (Simplification) ---

    size_t TOPO_MAX_LEN;    ///< Longueur max (en noeuds) pour supprimer une pointe basée uniquement sur la topologie.
    size_t RCTC_MAX_LEN;    ///< Longueur max (en noeuds) pour supprimer une pointe en considérant aussi la couverture.

    double RCTC_MAX_RATIO = 10;    ///< Ratio max utilisé pour le calcul interne.
    double TOPO_MAX_RATIO = 2.5;   ///< Facteur multiplicateur de k pour définir la longueur topologique max.

    /** * @brief Ratio de couverture (Parent / Pointe).
     * Si couverture(Parent) > couverture(Pointe) * RCTC_RATIO, la pointe est probablement une erreur.
     */
    double RCTC_RATIO = 2.0;

    /** * @brief Ratio de couverture pour choisir la branche principale lors de la génération de contigs.
     * Utilisé pour désambiguïser les bifurcations non résolues.
     */
    double COVERAGE_RATIO = 1.0;

    size_t MAX_CONTIG_LEN = 1000000;   ///< Sécurité : longueur maximale d'un contig généré.
    double SEARCH_DEPTH_FACTOR = 20.0; ///< Détermine la profondeur de recherche des bulles (profondeur = k * facteur).
    int MAX_PASSES = 50;               ///< Nombre maximum d'itérations pour les boucles de simplification.

    /**
     * @brief Constructeur de configuration.
     * Calcule les seuils dépendants de la taille des k-mers (k).
     * @param k_size La taille du k-mer utilisée pour le graphe.
     */
    GraphDBJConfig(int k_size) {
        // Initialisation dynamique : plus k est grand, plus on tolère de longues pointes.
        TOPO_MAX_LEN = (size_t)(k_size * TOPO_MAX_RATIO);
        RCTC_MAX_LEN = (size_t)(k_size * RCTC_MAX_RATIO);
    }
};

/**
 * @class GraphDBJ
 * @brief Classe principale gérant le graphe de De Bruijn.
 *
 * Responsabilités :
 * 1. Construction du graphe à partir du BitVector brut (séquences).
 * 2. Simplification du graphe (suppression des erreurs de séquençage : tips, bubbles).
 * 3. Génération des contigs initiaux (chemins simples).
 * 4. Fusion avancée des contigs (Overlap-Layout-Consensus simplifié).
 */
class GraphDBJ {
private:
    int k; ///< Taille du k-mer (Les noeuds représentent des k-1 mers).

    /** * @brief Table de hachage principale stockant tous les noeuds du graphe.
     * Clé : uint64_t (représentation binaire du k-1 mer).
     * Valeur : Pointeur vers l'objet Noeud.
     */
    std::unordered_map<uint64_t, Noeud*> nodes_map;

    GraphDBJConfig config; ///< Paramètres de l'algorithme.

    // --- Méthodes utilitaires internes ---

    /**
     * @brief Calcule le complément inverse d'une séquence encodée en entier.
     * Exemple : A(00) devient T(11), C(10) devient G(01), et l'ordre est inversé.
     * @complexity : space : O(1), time : O(k)
     * k = quantité nucléotides.
     */
    uint64_t getReverseComplement(uint64_t val, int k) const;

    /**
     * @brief Extrait la valeur entière d'un k-mer directement depuis le BitVector compressé.
     * @complexity : space : O(1), time : O(k)
     * k = quantité de nucléotide
     */
    uint64_t extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const;

    /**
     * @brief Convertit un entier k-mer en chaîne de caractères lisible (A, C, G, T).
     * @complexity : space : O(l), time : O(l)
     * l = longueur d'entrèe
     */
    std::string kmerToString(uint64_t val, int length) const;

    /**
     * @brief Supprime proprement les liens (arêtes) entre un parent et un enfant.
     * @complexity : space : O(1), time : O(d)
     * d = plus grand degré d'un noeud
     */
    void disconnectNodes(Noeud* parent, Noeud* child);

    /**
     * @brief Cherche un noeud de convergence entre deux branches (détection de bulle).
     * @complexity : space (), time ()
     */
    static Noeud* findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit);

    /**
     * @brief Ajoute une séquence k-mer (entière) à un BitVector (utilisé pour construire les contigs).
     */
    void addKmerToBitVector(BitVector& bv, uint64_t val, int length) const;

    // --- Méthodes statiques pour la fusion de contigs ---

    /**
     * @brief Calcule le meilleur chevauchement entre deux BitVectors.
     * @return La longueur du chevauchement (en nucléotides), ou 0 si aucun valide.
     */
    static int calculateBinaryOverlap(const BitVector& b1, const BitVector& b2, int min_len, double error_percent);

    /**
     * @brief Vérifie si le BitVector 'small' est entièrement contenu dans 'large' (avec tolérance d'erreur).
     */
    static bool isContained(const BitVector& large, const BitVector& small, double error_percent);

public:
    /**
     * @brief Constructeur : Construit le graphe à partir des données converties.
     * @param converter L'objet contenant le BitVector des lectures.
     * @param kmer_size La taille K.
     * @param config La configuration des seuils.
     */
    GraphDBJ(const Convert& converter, int kmer_size, const GraphDBJConfig& config);

    /** @brief Destructeur : Libère la mémoire des Noeuds. */
    ~GraphDBJ();

    /** @brief Récupère la liste de tous les pointeurs de noeuds (pour debug ou itération). */
    std::vector<Noeud*> getNodes() const;

    // --- Algorithmes d'Assemblage ---

    /**
     * @brief Supprime les "tips" (pointes mortes) causés par des erreurs en fin de lecture.
     * @return Nombre de pointes supprimées.
     */
    int clipTips();

    /**
     * @brief Résout les "bulles" (chemins divergents qui se rejoignent) causées par des SNPs ou erreurs.
     * @return Nombre de bulles résolues (branches supprimées).
     */
    int resolveBubbles();

    /**
     * @brief Parcourt le graphe simplifié pour extraire les contigs (séquences contiguës).
     * @return Vecteur de BitVectors représentant les contigs.
     */
    std::vector<BitVector> generateContigs() const;

    /**
     * @brief Exporte la structure du graphe au format GFA (Graph Fragment Assembly) pour visualisation (ex: Bandage).
     */
    void exportToGFA(const std::string& filename) const;

    /**
     * @brief Fusionne les contigs qui se chevauchent fortement (étape OLC post-graphe).
     * Utilise une approche de "Deep Seeding" et gère les inclusions.
     */
    static void mergeContigs(std::vector<BitVector>& contigs, int min_overlap, double overlap_error_percent, double contained_error_percent);
};

#endif //ASSEMBLEUR_GRAPHDBJ_H
