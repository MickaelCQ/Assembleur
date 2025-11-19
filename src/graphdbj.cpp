#include "graphdbj.h"
#include <iostream>
#include <stdexcept>

GraphDBJ::GraphDBJ(const Convert& converter, int kmer_size) : k(kmer_size) {
    if (k <= 1) {
        throw std::invalid_argument("K doit etre superieur a 1 pour construire un graphe.");
    }
    // Un uint64_t peut stocker jusqu'à 32 nucléotides (64 bits).
    // Donc k-1 doit être <= 32 => k <= 33.
    if (k > 33) {
        throw std::invalid_argument("Taille de K trop grande pour etre stockee sur 64 bits (max 33).");
    }

    const BitVector& bv = converter.get_bitVector();
    const std::vector<size_t>& read_ends = converter.get_read_end_positions();

    size_t current_read_start = 0;

    // 1. Parcourir chaque lecture
    for (size_t end_pos : read_ends) {
        // Calculer la longueur de la lecture en bits
        size_t read_len_bits = end_pos - current_read_start;
        size_t read_len_nuc = read_len_bits / 2;

        // Si la lecture est assez longue pour contenir au moins un k-mer
        if (read_len_nuc >= (size_t)k) {
            // Parcourir tous les k-mers de cette lecture
            // On s'arrête quand il ne reste plus assez de place pour un k-mer
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {

                // Position en bits du début du k-mer courant
                size_t bit_idx = current_read_start + (i * 2);

                // --- CONSTRUCTION NOEUD PREFIXE (k-1) ---
                // Prend k-1 nucléotides à partir de bit_idx
                uint64_t prefix_val = extractKmerValue(bv, bit_idx, k - 1);

                // --- CONSTRUCTION NOEUD SUFFIXE (k-1) ---
                // Le suffixe commence 1 nucléotide (2 bits) après le préfixe
                // et a aussi une longueur de k-1
                uint64_t suffix_val = extractKmerValue(bv, bit_idx + 2, k - 1);

                // 1. Récupérer ou créer le noeud Prefix
                Noeud* prefixNode;
                if (nodes_map.find(prefix_val) == nodes_map.end()) {
                    prefixNode = new Noeud(prefix_val);
                    nodes_map[prefix_val] = prefixNode;
                } else {
                    prefixNode = nodes_map[prefix_val];
                }
                prefixNode->coverage++;

                // 2. Récupérer ou créer le noeud Suffix
                Noeud* suffixNode;
                if (nodes_map.find(suffix_val) == nodes_map.end()) {
                    suffixNode = new Noeud(suffix_val);
                    nodes_map[suffix_val] = suffixNode;
                } else {
                    suffixNode = nodes_map[suffix_val];
                }

                // 3. Créer l'arête
                // Vérifier si l'arête existe déjà pour éviter les doublons (multigraphe vs graphe simple)
                bool edgeExists = false;
                for (auto* child : prefixNode->c) {
                    if (child == suffixNode) { edgeExists = true; break; }
                }

                if (!edgeExists) {
                    prefixNode->c.push_back(suffixNode);
                    suffixNode->parents.push_back(prefixNode); // <--- AJOUT IMPORTANT
                }
            }
        }

        // Mise à jour pour la prochaine lecture
        current_read_start = end_pos;
    }
}

GraphDBJ::~GraphDBJ() {
    // Nettoyage de la mémoire : delete de tous les pointeurs stockés dans la map
    for (auto& pair : nodes_map) {
        delete pair.second;
    }
    nodes_map.clear();
}

uint64_t GraphDBJ::extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const {
    uint64_t val = 0;

    // On parcourt nucléotide par nucléotide
    for (int i = 0; i < len_nucleotides; ++i) {
        size_t pos = start_bit_idx + (i * 2);

        // Lecture des 2 bits
        bool b1 = bv.test(pos);
        bool b2 = bv.test(pos + 1);

        // On décale la valeur actuelle de 2 bits vers la gauche pour faire de la place
        val = val << 2;

        // Construction de la valeur 2 bits (b1b2)
        // Si b1 est vrai, on ajoute 2 (10 en binaire)
        // Si b2 est vrai, on ajoute 1 (01 en binaire)
        if (b1) val |= 2ULL;
        if (b2) val |= 1ULL;
    }

    return val;
}

std::vector<Noeud*> GraphDBJ::getNodes() const {
    std::vector<Noeud*> result;
    result.reserve(nodes_map.size());

    for (const auto& pair : nodes_map) {
        result.push_back(pair.second);
    }
    return result;
}

void GraphDBJ::removeTips(int length_threshold) {
    std::cout << "--- Elagage des pointes (Tip Clipping) ---" << std::endl;
    bool changed = true;
    int tips_removed = 0;

    while (changed) {
        changed = false;
        auto all_nodes = getNodes(); // Copie des pointeurs pour itérer sûrement

        for (Noeud* n : all_nodes) {
            if (n->removed) continue;

            // Cas 1 : Fin de pointe (Dead end) -> In-degree > 0, Out-degree == 0
            if (n->c.empty() && !n->parents.empty()) {
                // On remonte en arrière
                std::vector<Noeud*> chain;
                Noeud* curr = n;
                bool keep_chain = false;

                // On construit la chaine inversée
                while (curr->parents.size() == 1 && curr->c.size() <= 1) {
                    chain.push_back(curr);
                    if (chain.size() > (size_t)length_threshold) {
                        keep_chain = true; // Trop long pour être une erreur
                        break;
                    }
                    curr = curr->parents[0];
                }
                // Ajouter le dernier noeud (le point d'embranchement) pour vérif
                if (!keep_chain && chain.size() <= (size_t)length_threshold) {
                    // C'est une petite pointe, on supprime tout sauf le point d'ancrage
                    for (Noeud* to_remove : chain) {
                        to_remove->removed = true;
                        // Déconnecter proprement
                        if (!to_remove->parents.empty()) {
                            disconnectNodes(to_remove->parents[0], to_remove);
                        }
                    }
                    tips_removed++;
                    changed = true;
                }
            }
            // Note: On pourrait aussi gérer les pointes d'entrée (In=0, Out>0) de façon symétrique
        }
    }
    std::cout << "Pointes supprimees : " << tips_removed << std::endl;
}

// Helper pour trouver le noeud de convergence
Noeud* GraphDBJ::findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit) {
    // Ceci est une version simplifiée. Pour une vraie implémentation,
    // il faudrait un double BFS. Ici, on regarde si branch1 rejoint branch2 rapidement.

    // Pour cet exemple, on suppose une bulle simple : S -> (chemin1) -> E  et S -> (chemin2) -> E
    if (branch1->c.size() == 1 && branch1->c[0] == branch2->c[0]) return branch1->c[0]; // Cas trivial
    if (branch1->c.size() == 1 && branch1->c[0] == branch2) return branch2; // Indel

    return nullptr;
}

void GraphDBJ::resolveBubbles() {
    std::cout << "--- Resolution des bulles (Avancee) ---" << std::endl;
    int bubbles_collapsed = 0;

    for (auto& pair : nodes_map) {
        Noeud* s = pair.second;
        if (s->removed || s->c.size() < 2) continue;

        // On a une divergence. Comparons les enfants deux à deux.
        // Simplification : on prend les deux premiers enfants
        Noeud* pathA = s->c[0];
        Noeud* pathB = s->c[1];

        // Cherchons s'ils se rejoignent immédiatement (SNP ou petite erreur)
        // Cas SNP : S->A->E et S->B->E.
        // A et B doivent avoir le même enfant unique E.
        if (pathA->c.size() == 1 && pathB->c.size() == 1) {
            if (pathA->c[0] == pathB->c[0]) {
                // C'est une bulle !
                // On garde celui avec la plus grosse couverture
                Noeud* to_keep = (pathA->coverage >= pathB->coverage) ? pathA : pathB;
                Noeud* to_remove = (pathA->coverage >= pathB->coverage) ? pathB : pathA;

                to_remove->removed = true;
                disconnectNodes(s, to_remove);
                disconnectNodes(to_remove, to_remove->c[0]);
                bubbles_collapsed++;
            }
        }
    }
    std::cout << "Bulles simplifiees : " << bubbles_collapsed << std::endl;
}

std::vector<std::string> GraphDBJ::generateContigs() const {
    std::vector<std::string> contigs;
    std::unordered_map<Noeud*, bool> visited;

    // Ratio pour considérer un chemin comme "clair" (ex: 10x plus fréquent que l'autre)
    const double COVERAGE_RATIO = 5.0;

    for (const auto& pair : nodes_map) {
        Noeud* startNode = pair.second;
        if (startNode->removed || visited[startNode]) continue;

        // On cherche un point de départ (bout de graphe ou embranchement non résolu)
        if (startNode->parents.empty() || startNode->parents.size() > 1) {

            std::string seq = kmerToString(startNode->p, k - 1);
            visited[startNode] = true;
            Noeud* curr = startNode;

            while (true) {
                Noeud* next = nullptr;

                // CAS 1 : Chemin simple (1 seul enfant)
                if (curr->c.size() == 1) {
                    next = curr->c[0];
                }
                // CAS 2 : Bifurcation (Plusieurs enfants) -> Heuristique de couverture
                else if (curr->c.size() > 1) {
                    Noeud* best_candidate = nullptr;
                    uint32_t max_cov = 0;
                    uint32_t second_max_cov = 0;

                    for (auto* child : curr->c) {
                        if (child->removed) continue;
                        if (child->coverage > max_cov) {
                            second_max_cov = max_cov;
                            max_cov = child->coverage;
                            best_candidate = child;
                        } else if (child->coverage > second_max_cov) {
                            second_max_cov = child->coverage;
                        }
                    }

                    // Si le meilleur candidat est beaucoup plus fort que le deuxième
                    if (best_candidate != nullptr && max_cov > (second_max_cov * COVERAGE_RATIO)) {
                        next = best_candidate;
                    } else {
                        // Ambiguïté trop forte (ex: répétition 50/50), on coupe le contig ici.
                        break;
                    }
                } else {
                    // 0 enfants, fin du chemin
                    break;
                }

                // Vérifications de sécurité sur 'next'
                if (next == nullptr || next->removed || visited[next]) break;

                // Ajouter le nucléotide
                uint64_t val = next->p;
                uint64_t last_nuc_val = val & 3ULL;
                char c = "ACGT"[last_nuc_val == 2 ? 1 : (last_nuc_val == 1 ? 2 : (last_nuc_val == 3 ? 3 : 0))]; // Petit hack de conversion rapide

                seq += c;
                visited[next] = true;
                curr = next;
            }
            contigs.push_back(seq);
        }
    }
    return contigs;
}

std::string GraphDBJ::kmerToString(uint64_t val, int length) const {
    std::string seq = "";
    // On lit les bits de poids fort vers poids faible (selon ta logique extractKmerValue)
    // Attention : dans extractKmerValue, tu fais val << 2.
    // Le premier nucléotide inséré se retrouve aux bits de poids fort.
    // Pour récupérer l'ordre correct :
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);

        if (two_bits == 0) seq += 'A';       // 00
        else if (two_bits == 2) seq += 'C';  // 10 (Attention: ton code: b1=1(2), b2=0(0) => 10 = 2)
        else if (two_bits == 1) seq += 'G';  // 01 (b1=0, b2=1 => 01 = 1)
        else seq += 'T';                     // 11 (3)
    }
    return seq;
}

void GraphDBJ::disconnectNodes(Noeud* parent, Noeud* child) {
    // Retirer child des enfants de parent
    auto it_c = std::find(parent->c.begin(), parent->c.end(), child);
    if (it_c != parent->c.end()) parent->c.erase(it_c);

    // Retirer parent des parents de child
    auto it_p = std::find(child->parents.begin(), child->parents.end(), parent);
    if (it_p != child->parents.end()) child->parents.erase(it_p);
}