#include "graphdbj.h"
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>

// Fonction utilitaire pour calculer le chevauchement max entre la fin de s1 et le début de s2
int calculateOverlap(const std::string& s1, const std::string& s2, int min_len) {
    int max_ov = 0;
    // On ne teste pas tout, on limite la zone de recherche pour la performance
    int limit = std::min(s1.length(), s2.length());

    for (int len = min_len; len < limit; ++len) {
        // Vérifie si le suffixe de s1 de taille 'len' est égal au préfixe de s2
        if (s1.compare(s1.length() - len, len, s2, 0, len) == 0) {
            max_ov = len;
        }
    }
    return max_ov;
}

std::vector<std::string> GraphDBJ::mergeContigs(std::vector<std::string> contigs, int min_overlap) {
    std::cout << "--- Post-traitement : Fusion des Contigs (Consensus) ---" << std::endl;

    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<std::string> next_pass;
        std::vector<bool> merged(contigs.size(), false);

        // Tri par taille décroissante : on veut garder les gros morceaux comme base
        std::sort(contigs.begin(), contigs.end(), [](const std::string& a, const std::string& b) {
            return a.length() > b.length();
        });

        for (size_t i = 0; i < contigs.size(); ++i) {
            if (merged[i]) continue;

            std::string current = contigs[i];
            bool merged_current = false;

            for (size_t j = 0; j < contigs.size(); ++j) {
                if (i == j || merged[j]) continue;

                // 1. INCLUSION : Si contigs[j] est DANS current, on supprime contigs[j]
                if (current.find(contigs[j]) != std::string::npos) {
                    merged[j] = true; // J est absorbé par I
                    changed = true;
                    continue;
                }

                // 2. FUSION (Overlap)
                // Cas A : Fin de CURRENT chevauche Début de J
                int ov_right = calculateOverlap(current, contigs[j], min_overlap);
                if (ov_right > 0) {
                    // Fusion : Current + (J sans le chevauchement)
                    current += contigs[j].substr(ov_right);
                    merged[j] = true;
                    merged_current = true;
                    changed = true;
                    // On recommence la boucle avec le nouveau 'current' agrandi
                    i--; // Astuce pour re-traiter 'current' au prochain tour de i
                    break;
                }

                // Cas B : Fin de J chevauche Début de CURRENT
                int ov_left = calculateOverlap(contigs[j], current, min_overlap);
                if (ov_left > 0) {
                    // Fusion : J + (Current sans le chevauchement)
                    current = contigs[j] + current.substr(ov_left);
                    merged[j] = true;
                    merged_current = true;
                    changed = true;
                    i--;
                    break;
                }
            }

            if (!merged_current && !merged[i]) {
                next_pass.push_back(current);
            } else if (merged_current) {
                // Si on a fusionné, on ajoute le résultat (qui est dans 'current')
                // Note: si on a fait i--, on ne l'ajoute pas tout de suite, mais ici on simplifie la logique
                // Dans cette version simple, on ajoute le merged et on relance le while(changed)
                next_pass.push_back(current);
                merged[i] = true;
            }
        }
        contigs = next_pass;
        if(changed) std::cout << "Cycle de fusion termine. Contigs restants : " << contigs.size() << std::endl;
    }

    return contigs;
}

void GraphDBJ::exportToGFA(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return;

    out << "H\tVN:Z:1.0\n";

    // 1. Créer une carte pour traduire k-mer (uint64_t) -> ID simple (1, 2, 3...)
    std::unordered_map<uint64_t, int> id_map;
    int current_id = 1;

    // Assignation des IDs
    for (const auto& pair : nodes_map) {
        if (!pair.second->removed) {
            id_map[pair.first] = current_id++;
        }
    }

    // 2. Écrire les segments avec les nouveaux IDs
    for (const auto& pair : nodes_map) {
        Noeud* node = pair.second;
        if (node->removed) continue;

        std::string seq = kmerToString(node->p, k - 1);
        // On utilise id_map[...] pour l'ID
        out << "S\t" << id_map[node->p] << "\t" << seq << "\tDP:i:" << node->coverage << "\n";
    }

    // 3. Écrire les liens avec les nouveaux IDs
    for (const auto& pair : nodes_map) {
        Noeud* src = pair.second;
        if (src->removed) continue;

        for (Noeud* dest : src->c) {
            if (dest->removed) continue;

            // Vérification de sécurité
            if (id_map.find(src->p) != id_map.end() && id_map.find(dest->p) != id_map.end()) {
                out << "L\t" << id_map[src->p] << "\t+\t" << id_map[dest->p] << "\t+\t" << (k - 2) << "M\n";
            }
        }
    }

    out.close();
}

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
        size_t read_len_nuc = (end_pos - current_read_start) / 2;
        if (read_len_nuc >= (size_t)k) {
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {

                // 1. Récupérer les valeurs brutes (Forward)
                uint64_t u_fwd = extractKmerValue(bv, current_read_start + i * 2, k - 1);
                uint64_t v_fwd = extractKmerValue(bv, current_read_start + i * 2 + 2, k - 1);

                // 2. Calculer les valeurs inverses (Reverse Complement)
                // Attention : si U -> V, alors Rev(V) -> Rev(U)
                uint64_t v_rev = getReverseComplement(v_fwd, k - 1);
                uint64_t u_rev = getReverseComplement(u_fwd, k - 1);

                // 3. Fonction lambda pour ajouter une arête (évite de dupliquer le code)
                auto add_edge = [&](uint64_t source_val, uint64_t target_val) {
                    Noeud*& sourceNode = nodes_map[source_val];
                    if (!sourceNode) sourceNode = new Noeud(source_val);

                    Noeud*& targetNode = nodes_map[target_val];
                    if (!targetNode) targetNode = new Noeud(target_val);

                    sourceNode->coverage++;
                    targetNode->coverage++;

                    bool edgeExists = false;
                    for (auto* child : sourceNode->c) {
                        if (child == targetNode) { edgeExists = true; break; }
                    }
                    if (!edgeExists) {
                        sourceNode->c.push_back(targetNode);
                        targetNode->parents.push_back(sourceNode);
                    }
                };

                // 4. Ajouter l'arête sur le brin SENS (Forward)
                add_edge(u_fwd, v_fwd);

                // 5. Ajouter l'arête sur le brin COMPLÉMENTAIRE (Reverse)
                // C'est ça qui assure la continuité du graphe dans les deux sens !
                add_edge(v_rev, u_rev);
            }
        }
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

uint64_t GraphDBJ::getReverseComplement(uint64_t val, int length) const {
    uint64_t res = 0;

    // On parcourt chaque nucléotide du k-mer original
    for (int i = 0; i < length; ++i) {
        // Récupère le nucléotide le plus à droite (les 2 derniers bits)
        // Dans votre construction, les bits de poids faible sont les derniers ajoutés (fin de séquence)
        uint64_t nuc = val & 3ULL;

        // Calcule le complément :
        // A(00) <-> T(11) => 00 XOR 11 = 11
        // C(10) <-> G(01) => 10 XOR 11 = 01
        // L'opération XOR 3 (binaire 11) inverse les bits correctement pour votre encodage.
        uint64_t comp = nuc ^ 3ULL;

        // Placement dans le résultat :
        // Le dernier nucléotide de 'val' devient le premier de 'res' (Reverse)
        // On le décale donc vers la position de poids fort correspondant à l'index i inversé.
        res |= (comp << (2 * (length - 1 - i)));

        // On passe au nucléotide suivant en décalant val
        val >>= 2;
    }
    return res;
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

int GraphDBJ::clipTips() {
    std::cout << "--- Elagage des pointes (Algorithme Minia : Topo + RCTC) ---" << std::endl;
    int tips_removed_count = 0;
    bool changed = true;

    // Paramètres inspirés de Minia
    // 1. Critère Topologique : On coupe tout ce qui est <= k * 2.5
    const int TOPO_MAX_LEN = (int)(k * 2.5);

    // 2. Critère Relatif (RCTC) : On coupe jusqu'à k * 10 SI la couverture est faible
    const int RCTC_MAX_LEN = k * 10;

    // 3. Ratio de couverture : Le voisin (main path) doit être 2x plus couvert que le tip
    const double RCTC_RATIO = 2.0;

    while (changed) {
        changed = false;

        // Copie des clés pour itérer
        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);

        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue;
            Noeud* node = nodes_map[key];
            if (node->removed) continue;

            // On cherche les "Dead Ends" (Noeuds sans enfants)
            if (node->c.empty()) {

                // --- Phase de remontée (Backtracking) ---
                // On remonte depuis la pointe pour mesurer sa longueur et sa couverture moyenne
                // jusqu'à trouver un embranchement (le tronc principal).

                std::vector<Noeud*> tip_path;
                Noeud* curr = node;
                bool is_valid_tip = true;
                uint64_t sum_coverage = 0;

                // On remonte tant que :
                // 1. On n'est pas allé trop loin (limite max absolue)
                // 2. Le noeud a un seul parent (c'est une ligne droite)
                // 3. Le noeud parent n'a qu'un seul enfant (c'est encore la ligne du tip)
                //    Dès que le parent a >1 enfant, c'est le point d'ancrage (Junction), on s'arrête.

                while (tip_path.size() <= RCTC_MAX_LEN) {
                    tip_path.push_back(curr);
                    sum_coverage += curr->coverage;

                    if (curr->parents.empty()) {
                        // C'est un îlot isolé (pas connecté au reste), on peut le supprimer s'il est court
                        // mais ce n'est pas une "pointe" attachée à un graphe.
                        is_valid_tip = false;
                        break;
                    }

                    Noeud* parent = curr->parents[0];

                    // Si le parent a plusieurs enfants, c'est le point d'ancrage (la jonction).
                    // On a trouvé la base de la pointe.
                    if (parent->c.size() > 1) {
                        break; // Fin de la remontée
                    }

                    // Si le parent a plusieurs parents, c'est une structure complexe (bulle?), on arrête.
                    if (parent->parents.size() != 1) {
                        is_valid_tip = false;
                        break;
                    }

                    curr = parent;
                }

                if (!is_valid_tip || tip_path.empty()) continue;

                // Le noeud d'ancrage est le parent du dernier noeud visité dans la boucle
                Noeud* last_tip_node = tip_path.back();
                if (last_tip_node->parents.empty()) continue;
                Noeud* anchor = last_tip_node->parents[0];

                // --- Prise de décision (Minia Logic) ---
                int len = tip_path.size();
                double tip_avg_cov = (double)sum_coverage / len;
                double anchor_cov = (double)anchor->coverage; // Couverture du "Main Path"

                bool should_remove = false;

                // Règle 1 : Topologique (Court et inconditionnel)
                if (len <= TOPO_MAX_LEN) {
                    should_remove = true;
                }
                // Règle 2 : Relative Coverage (Plus long mais faible couverture comparé au voisin)
                else if (len <= RCTC_MAX_LEN) {
                    if (anchor_cov > (tip_avg_cov * RCTC_RATIO)) {
                        should_remove = true;
                    }
                }

                // --- Suppression ---
                if (should_remove) {
                    // 1. Couper le lien entre l'ancre et le début du tip
                    disconnectNodes(anchor, last_tip_node);

                    // 2. Marquer tout le tip comme supprimé
                    for (Noeud* n : tip_path) {
                        n->removed = true;
                    }

                    tips_removed_count++;
                    changed = true;
                }
            }
        }
    }
    std::cout << "Pointes supprimees : " << tips_removed_count << std::endl;
    return tips_removed_count;
}

// Fonction utilitaire à ajouter si elle n'existe pas encore
void GraphDBJ::disconnectNodes(Noeud* parent, Noeud* child) {
    if (!parent || !child) return;

    // Supprimer le child de la liste des enfants du parent
    parent->c.erase(std::remove(parent->c.begin(), parent->c.end(), child), parent->c.end());

    // Supprimer le parent de la liste des parents du child
    child->parents.erase(std::remove(child->parents.begin(), child->parents.end(), parent), child->parents.end());
}

Noeud* GraphDBJ::findConvergence(Noeud* branch1, Noeud* branch2, int depth_limit) {
    // Version simplifiée : vérifie si les branches se rejoignent immédiatement.
    // Pour une implémentation complète, il faudrait un BFS (Breadth-First Search).
    if (branch1->c.size() == 1 && branch1->c[0] == branch2->c[0]) return branch1->c[0]; // Convergence en V
    if (branch1->c.size() == 1 && branch1->c[0] == branch2) return branch2; // Indel (Insertion/Deletion)
    return nullptr;
}

int GraphDBJ::resolveBubbles() {
    // On augmente la profondeur de recherche pour attraper les grosses bulles
    // k + 2 (~33) était trop court. k * 2 (~62) permet de fermer les grosses boucles.
    int search_depth = k * 20;

    std::cout << "--- Resolution des bulles (Profondeur " << search_depth << ") ---" << std::endl;
    int bubbles_collapsed = 0;
    bool changed = true;

    int max_passes = 50;
    int pass = 0;

    while (changed && pass < max_passes) {
        changed = false;
        pass++;

        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);

        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue;

            Noeud* s = nodes_map[key];
            if (s->removed || s->c.size() != 2) continue;

            Noeud* branch1 = s->c[0];
            Noeud* branch2 = s->c[1];

            std::vector<Noeud*> path1, path2;
            path1.push_back(branch1);
            path2.push_back(branch2);

            Noeud* convergence = nullptr;

            // 3. CHANGEMENT ICI : On utilise search_depth au lieu de (k + 2)
            for (int step = 0; step < search_depth; ++step) {

                // --- (Le reste de la boucle for ne change pas) ---

                Noeud* last1 = path1.back();
                if (last1->c.size() == 1 && !last1->c[0]->removed) {
                    path1.push_back(last1->c[0]);
                }

                Noeud* last2 = path2.back();
                if (last2->c.size() == 1 && !last2->c[0]->removed) {
                    path2.push_back(last2->c[0]);
                }

                for(Noeud* n2 : path2) {
                    if (path1.back() == n2) { convergence = n2; break; }
                }
                if(convergence) break;

                for(Noeud* n1 : path1) {
                    if (path2.back() == n1) { convergence = n1; break; }
                }
                if(convergence) break;

                if (path1.back() == last1 && path2.back() == last2) break;
            }

            // --- Simplification (Rien ne change ici non plus) ---
            if (convergence && convergence != branch1 && convergence != branch2) {
                uint32_t cov1 = 0;
                for(auto* n : path1) { if(n==convergence) break; cov1 += n->coverage; }

                uint32_t cov2 = 0;
                for(auto* n : path2) { if(n==convergence) break; cov2 += n->coverage; }

                std::vector<Noeud*>& to_remove = (cov1 >= cov2) ? path2 : path1;
                Noeud* branch_to_keep = (cov1 >= cov2) ? branch1 : branch2;

                for(auto* n : to_remove) {
                    if (n != convergence) {
                        n->removed = true;
                    }
                }

                Noeud* branch_to_remove = (cov1 >= cov2) ? branch2 : branch1;
                disconnectNodes(s, branch_to_remove);

                bool connected = false;
                for(auto* child : s->c) if(child == branch_to_keep) connected = true;
                if(!connected) s->c.push_back(branch_to_keep);

                bubbles_collapsed++;
                changed = true;
            }
        }
    }
    std::cout << "Bulles simplifiees : " << bubbles_collapsed << std::endl;

    return bubbles_collapsed;
}

std::vector<std::string> GraphDBJ::generateContigs() const {
    std::vector<std::string> contigs;
    std::unordered_map<Noeud*, bool> visited;

    // Seuil heuristique : si un chemin est 5x plus couvert que l'autre, on le suit.
    // Je pense que ça on pourrait le mettre en option, faudrait le tester
    const double COVERAGE_RATIO = 1.0;

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

        // Correspondance inverse à celle définie dans addCha()
        if (two_bits == 0) seq += 'A';       // 00
        else if (two_bits == 2) seq += 'C';  // 10
        else if (two_bits == 1) seq += 'G';  // 01
        else seq += 'T';                     // 11
    }
    return seq;
}
