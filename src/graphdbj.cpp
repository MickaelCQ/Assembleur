#include "graphdbj.h"

#include <complex>
#include <iostream>
#include <stdexcept>

GraphDBJ::GraphDBJ(const Convert& converter, int kmer_size) : k(kmer_size) {
    if (k < 2)
      {
        throw std::invalid_argument("K doit etre superieur a 2 pour construire un graphe.");    
      }
    // Un uint64_t peut stocker jusqu'à 32 nucléotides (64 bits). Donc k-1 doit être <= 32 => k < 32.
    if (k > 32) 
      {
        throw std::invalid_argument("Taille de K trop grande pour etre stockee sur 64 bits (max 33).");
      }

    const BitVector& bv = converter.get_bitVector();
    const std::vector<size_t>& read_ends = converter.get_read_end_positions();
    size_t current_read_start = 0;


    // 1. Parcourir des reats 
    for (size_t end_pos : read_ends) {
        // Calculer la longueur de la lecture en bits
      size_t read_len_nuc = (end_pos - current_read_start) / 2;

        // Si la lecture est assez longue pour contenir au moins un k-mer
        if (read_len_nuc >= (size_t)k) {
            // Sliding Window : Parcourir tous les k-mers de cette lecture
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {

		//Mickael : dégagage/simplification de la variable intermédiaire pr déterminer les suf/préfixes
		uint64_t prefix_val = extractKmerValue(bv, current_read_start + i * 2, k - 1);
		uint64_t suffix_val = extractKmerValue(bv, current_read_start + i * 2 + 2, k - 1);

        	 // 1. Récupérer ou créer le noeud Prefix
        	 //Mickael: Simplification des conditionnels en dégagant les else (inutile)
                 Noeud*& prefixNode = nodes_map[prefix_val];
                if (!prefixNode)
               {
    		    prefixNode = new Noeud(prefix_val);
                }
                 
                // 2. Récupérer ou créer les noeuds préfixe & suffixe
                Noeud*& suffixNode = nodes_map[suffix_val];
                if (!suffixNode)
                {
                    suffixNode = new Noeud(suffix_val);
                }

        	prefixNode->coverage++, suffixNode->coverage++;
        
                // 3. Créer l'arête
                // Vérifier si l'arête existe déjà pour éviter les doublons (multigraphe vs graphe simple)
                bool edgeExists = false;
                for (auto* child : prefixNode->c) {
                    if (child == suffixNode) { edgeExists = true; break;}
                }

                if (!edgeExists) {
                    prefixNode->c.push_back(suffixNode);
                    suffixNode->parents.push_back(prefixNode); // Maintien du lien parent pour le retour arrière
                }
            }
        }

        // Mise à jour du pointeur de début pour la prochaine lecture
        current_read_start = end_pos;
    }
}

GraphDBJ::~GraphDBJ() {
    // Nettoyage mémoire : suppression explicite de tous les pointeurs stockés dans la map
    for (auto& pair : nodes_map) {
        delete pair.second;
    }
    nodes_map.clear();
}

uint64_t GraphDBJ::extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const {
    uint64_t val = 0;

    // Construction de l'entier par décalage de bits (Bit Shifting)
    for (int i = 0; i < len_nucleotides; ++i) {
        size_t pos = start_bit_idx + (i * 2);

        bool b1 = bv.test(pos);
        bool b2 = bv.test(pos + 1);

        // Décale la valeur actuelle de 2 bits vers la gauche pour faire de la place aux nouveaux bits
        val = val << 2;

        // Ajoute la valeur 2 bits (b1b2) à la fin
        if (b1) val |= 2ULL; // Si b1 est 1, ajoute 10 (binaire) soit 2
        if (b2) val |= 1ULL; // Si b2 est 1, ajoute 01 (binaire) soit 1
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

    // Répète tant qu'on trouve des pointes (car supprimer une pointe peut en révéler une autre)
    while (changed) {
        changed = false;
        auto all_nodes = getNodes(); // Copie des pointeurs pour itérer sûrement

        for (Noeud* n : all_nodes) {
            if (n->removed) continue;

            // Définition d'une pointe de fin (Dead end) : Degré Entrant > 0, Degré Sortant == 0
            if (n->c.empty() && !n->parents.empty()) {
                // On remonte en arrière pour voir la longueur de la branche
                std::vector<Noeud*> chain;
                Noeud* curr = n;
                bool keep_chain = false;

                // On construit la chaine inversée tant qu'on est sur un chemin linéaire (1 parent, 1 enfant max)
                while (curr->parents.size() == 1 && curr->c.size() <= 1) {
                    chain.push_back(curr);
                    // Si la chaine est trop longue, c'est probablement une vraie séquence, pas une erreur
                    if (chain.size() > (size_t)length_threshold) {
                        keep_chain = true;
                        break;
                    }
                    curr = curr->parents[0];
                }

                // Si c'est une petite pointe (erreur probable), on supprime
                if (!keep_chain && chain.size() <= (size_t)length_threshold) {
                    for (Noeud* to_remove : chain) {
                        to_remove->removed = true; // Marquer comme supprimé
                        // Déconnecter proprement du reste du graphe
                        if (!to_remove->parents.empty()) {
                            disconnectNodes(to_remove->parents[0], to_remove);
                        }
                    }
                    tips_removed++;
                    changed = true;
                }
            }
        }
    }
    std::cout << "Pointes supprimees : " << tips_removed << std::endl;
}

Noeud* GraphDBJ::findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit) {
    // Version simplifiée : vérifie si les branches se rejoignent immédiatement.
    // Pour une implémentation complète, il faudrait un BFS (Breadth-First Search).
    if (branch1->c.size() == 1 && branch1->c[0] == branch2->c[0]) return branch1->c[0]; // Convergence en V
    if (branch1->c.size() == 1 && branch1->c[0] == branch2) return branch2; // Indel (Insertion/Deletion)
    return nullptr;
}

void GraphDBJ::resolveBubbles() {
    std::cout << "--- Resolution des bulles (Avancee) ---" << std::endl;
    int bubbles_collapsed = 0;

    for (auto& pair : nodes_map) {
        Noeud* s = pair.second;
        if (s->removed || s->c.size() < 2) continue;

        // On a une divergence (noeud avec >1 enfant).
        // Simplification : on regarde les deux premiers enfants.
        Noeud* pathA = s->c[0];
        Noeud* pathB = s->c[1];

        // Cas SNP typique : S->A->E et S->B->E (Losange simple)
        if (pathA->c.size() == 1 && pathB->c.size() == 1) {
            if (pathA->c[0] == pathB->c[0]) {
                // C'est une bulle ! On garde le chemin avec la plus forte couverture.
                Noeud* to_keep = (pathA->coverage >= pathB->coverage) ? pathA : pathB;
                Noeud* to_remove = (pathA->coverage >= pathB->coverage) ? pathB : pathA;

                to_remove->removed = true;
                // Déconnexion propre
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

    // Seuil heuristique : si un chemin est 5x plus couvert que l'autre, on le suit.
    // Je pense que ça on pourrait le mettre en option, faudrait le tester
    const double COVERAGE_RATIO = 5.0;

    for (const auto& pair : nodes_map) {
        Noeud* startNode = pair.second;
        if (startNode->removed || visited[startNode]) continue;

        // On cherche un point de départ valide pour un contig :
        // Soit un début absolu (pas de parents), soit une jonction complexe non résolue.
        if (startNode->parents.empty() || startNode->parents.size() > 1) {

            std::string seq = kmerToString(startNode->p, k - 1);
            visited[startNode] = true;
            Noeud* curr = startNode;

            // Extension gloutonne du contig
            while (true) {
                Noeud* next = nullptr;

                // CAS 1 : Chemin linéaire simple (1 seul enfant) -> On avance
                if (curr->c.size() == 1) {
                    next = curr->c[0];
                }
                // CAS 2 : Bifurcation (Plusieurs enfants) -> On utilise la couverture
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

                    // Si le meilleur candidat domine largement le deuxième (ratio > 5)
                    if (best_candidate != nullptr && max_cov > (second_max_cov * COVERAGE_RATIO)) {
                        next = best_candidate;
                    } else {
                        // Ambiguïté trop forte (ex: répétition 50/50), on arrête le contig ici.
                        break;
                    }
                } else {
                    // 0 enfants, fin du chemin
                    break;
                }

                // Vérifications de sécurité
                if (next == nullptr || next->removed || visited[next]) break;

                // Ajouter le dernier nucléotide du nœud suivant à la séquence
                uint64_t val = next->p;
                uint64_t last_nuc_val = val & 3ULL; // Masque pour récupérer les 2 derniers bits

                // Conversion 2 bits -> char (hack efficace via table de correspondance implicite)
                // 'A'=00, 'C'=10, 'G'=01, 'T'=11
                char c = "ACGT"[last_nuc_val == 2 ? 1 : (last_nuc_val == 1 ? 2 : (last_nuc_val == 3 ? 3 : 0))];

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
    // On lit les bits de poids fort vers poids faible (selon la logique d'extractKmerValue).
    // Le premier nucléotide inséré se retrouve aux bits de poids fort. Pour récupérer l'ordre correct,
    // on extrait les paires de bits de la position la plus significative vers la moins significative.
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);
std::
        // Correspondance inverse à celle définie dans addCha()
        if (two_bits == 0) seq += 'A';       // 00
        else if (two_bits == 2) seq += 'C';  // 10
        else if (two_bits == 1) seq += 'G';  // 01
        else seq += 'T';                     // 11
    }
    return seq;
}

void GraphDBJ::disconnectNodes(Noeud* parent, Noeud* child) {
    // Retirer child de la liste des enfants de parent
    auto it_c = std::find(parent->c.begin(), parent->c.end(), child);
    if (it_c != parent->c.end()) parent->c.erase(it_c);

    // Retirer parent de la liste des parents de child
    auto it_p = std::find(child->parents.begin(), child->parents.end(), parent);
    if (it_p != child->parents.end()) child->parents.erase(it_p);
}
