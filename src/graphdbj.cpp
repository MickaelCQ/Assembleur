#include "graphdbj.h"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <algorithm>

/**
 * @brief Helper : Inverse et complémente un BitVector complet.
 * Utile pour gérer le double brin de l'ADN (les lectures peuvent venir du brin Forward ou Reverse).
 */
static BitVector getBitVectorReverseComplement(const BitVector& seq) {
    BitVector rc;
    rc.reserve(seq.size());
    size_t len = seq.size() / 2;

    // Parcours inversé des nucléotides (de la fin vers le début)
    for (size_t i = 0; i < len; ++i) {
        size_t original_idx = len - 1 - i;
        uint8_t nuc = seq.getNucleotideAt(original_idx);

        // Opération bit à bit pour le complément :
        // 00 (A) ^ 11 (3) = 11 (T)
        // 10 (C) ^ 11 (3) = 01 (G)
        uint8_t comp = nuc ^ 3;

        // On repousse les bits dans le nouveau vecteur
        rc.push_back((comp & 2) != 0); // Bit de poids fort
        rc.push_back((comp & 1) != 0); // Bit de poids faible
    }
    return rc;
}

/**
 * @brief Calcule la longueur du chevauchement maximal entre la fin de b1 et le début de b2.
 * @details Cette fonction teste toutes les longueurs de chevauchement possibles de min_len à la taille min des vecteurs.
 * Elle tolère un certain pourcentage d'erreurs (mismatches).
 */
int GraphDBJ::calculateBinaryOverlap(const BitVector& b1, const BitVector& b2, int min_len, double error_percent) {
    int max_ov = 0;
    size_t len1 = b1.size() / 2;
    size_t len2 = b2.size() / 2;
    size_t limit = std::min(len1, len2);

    // On teste chaque longueur de chevauchement potentielle 'len'
    for (size_t len = min_len; len <= limit; ++len) {
        int mismatches = 0;

        // Calcul dynamique du seuil d'erreur autorisé pour cette longueur spécifique
        int allowed_errors = (int)(len * error_percent) + 1; // +1 offre une tolérance minimale de base

        bool match = true;

        // Index de départ dans b1 (on regarde la fin de b1)
        size_t start_idx_b1 = (len1 - len);

        // Comparaison nucléotide par nucléotide
        for (size_t i = 0; i < len; ++i) {
            uint8_t n1 = b1.getNucleotideAt(start_idx_b1 + i);
            uint8_t n2 = b2.getNucleotideAt(i);

            if (n1 != n2) {
                mismatches++;
                // Optimisation : arrêt précoce si trop d'erreurs
                if (mismatches > allowed_errors) {
                    match = false;
                    break;
                }
            }
        }

        // Si le chevauchement est valide, on le mémorise (on cherche le plus long)
        if (match) {
            max_ov = len;
        }
    }
    return max_ov;
}

/**
 * @brief Vérifie si le BitVector 'small' est inclus (contained) quelque part dans 'large'.
 * @details Effectue une recherche de sous-chaîne avec tolérance d'erreurs.
 */
bool GraphDBJ::isContained(const BitVector& large, const BitVector& small, double error_percent) {
    size_t l_len = large.size() / 2;
    size_t s_len = small.size() / 2;

    if (s_len > l_len) return false;

    // Seuil d'erreur basé sur la longueur du fragment inséré
    size_t max_mismatches = (size_t)(s_len * error_percent) + 1;

    // Fenêtre glissante sur 'large'
    for (size_t i = 0; i <= l_len - s_len; ++i) {
        size_t mismatches = 0;
        bool possible = true;

        for (size_t j = 0; j < s_len; ++j) {
            if (large.getNucleotideAt(i + j) != small.getNucleotideAt(j)) {
                mismatches++;
                if (mismatches > max_mismatches) {
                    possible = false;
                    break;
                }
            }
        }

        if (possible) return true; // Trouvé une correspondance valide
    }
    return false;
}

/**
 * @brief Extrait une clé k-mer (uint64) à partir d'un BitVector pour l'indexation.
 * Utilisé par mergeContigs pour créer les hash maps de seeds.
 */
static uint64_t getKmerKey(const BitVector& bv, size_t start_nuc_idx, int k) {
    uint64_t val = 0;
    size_t max_len = bv.size() / 2;

    // Sécurité : ne pas lire hors bornes
    if (start_nuc_idx + k > max_len) return 0;

    for (int i = 0; i < k; ++i) {
        uint8_t nuc = bv.getNucleotideAt(start_nuc_idx + i);
        val = (val << 2) | nuc; // Décalage et ajout des 2 bits
    }
    return val;
}

/**
 * @brief Vérifie l'inclusion à une position précise (optimisation pour mergeContigs).
 */
static bool checkContainmentAt(const BitVector& master, const BitVector& candidate, size_t master_offset, double error_percent) {
    size_t cand_len = candidate.size() / 2;
    size_t master_len = master.size() / 2;

    if (master_offset + cand_len > master_len) return false;

    size_t max_mismatches = (size_t)(cand_len * error_percent) + 1;
    size_t mismatches = 0;

    for (size_t i = 0; i < cand_len; ++i) {
        uint8_t m_nuc = master.getNucleotideAt(master_offset + i);
        uint8_t c_nuc = candidate.getNucleotideAt(i);

        if (m_nuc != c_nuc) {
            mismatches++;
            if (mismatches > max_mismatches) return false;
        }
    }
    return true;
}

/**
 * @brief Fusionne les contigs en utilisant une approche "Seed & Extend".
 * * Algorithme :
 * 1. Phase d'Inclusion : Détecte les petits contigs totalement inclus dans les grands et les marque comme 'absorbed'.
 * 2. Phase d'Extension : Tente d'étendre les bouts des contigs restants en trouvant des chevauchements.
 * Gère les orientations (Forward/Reverse) pour permettre la fusion même si un contig a été assemblé à l'envers.
 */
void GraphDBJ::mergeContigs(std::vector<BitVector>& contigs, int min_overlap, double overlap_error_percent, double contained_error_percent) {
    std::cout << "--- Post-traitement : Fusion 'Deep Seeding' (Bidirectionnelle) ---" << std::endl;

    // Taille du k-mer utilisé pour l'indexation (seed)
    int index_k = (min_overlap > 31) ? 31 : min_overlap;

    // Paramètres pour limiter la recherche (optimisation performance)
    const size_t MAX_SCAN_DEPTH = 5000; // Ne cherche pas de seed au-delà de 5000pb du bout
    const size_t MAX_SEED_DEPTH = 1500;
    const size_t SEED_STRIDE = index_k; // Pas d'échantillonnage des seeds

    std::vector<bool> absorbed(contigs.size(), false);

    // ==============================================================================
    // PHASE 1: CONTAINMENT (Inclusion)
    // ==============================================================================
    {
        // Structure pour stocker où se trouve un k-mer donné
        struct CandInfo { size_t id; bool is_rc; };
        std::unordered_map<uint64_t, std::vector<CandInfo>> start_map;

        // Indexation du DÉBUT de chaque contig
        for (size_t i = 0; i < contigs.size(); ++i) {
            size_t len = contigs[i].size() / 2;
            if (len < (size_t)index_k) continue;

            // Stocke le k-mer de début (Normal et Reverse Complement)
            start_map[getKmerKey(contigs[i], 0, index_k)].push_back({i, false});
            BitVector rc = getBitVectorReverseComplement(contigs[i]);
            start_map[getKmerKey(rc, 0, index_k)].push_back({i, true});
        }

        // Scan de tous les contigs pour voir s'ils contiennent le début d'un autre
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (absorbed[i]) continue;
            size_t m_len = contigs[i].size() / 2;

            for (size_t pos = 0; pos <= m_len - index_k; ++pos) {
                uint64_t key = getKmerKey(contigs[i], pos, index_k);
                if (start_map.find(key) == start_map.end()) continue;

                // Potentiels candidats à l'inclusion
                for (const auto& cand : start_map[key]) {
                    size_t j = cand.id;
                    if (i == j || absorbed[j]) continue;

                    bool is_contained = false;
                    // Vérification rigoureuse bit à bit
                    if (!cand.is_rc) is_contained = checkContainmentAt(contigs[i], contigs[j], pos, contained_error_percent);
                    else {
                        BitVector rc = getBitVectorReverseComplement(contigs[j]);
                        is_contained = checkContainmentAt(contigs[i], rc, pos, contained_error_percent);
                    }

                    if (is_contained) absorbed[j] = true;
                }
            }
        }
    }

    // ==============================================================================
    // PHASE 2: EXTENSION (Overlap)
    // ==============================================================================
    bool global_change = true;
    while (global_change) {
        global_change = false;

        // Indexation dense des débuts de contigs (Seeds)
        struct SeedInfo { size_t id; bool is_rc; size_t offset_in_cand; };
        std::unordered_map<uint64_t, std::vector<SeedInfo>> seed_map;

        // 1. Construction de l'Index
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (absorbed[i]) continue;
            size_t len = contigs[i].size() / 2;
            if (len < (size_t)index_k) continue;

            // On indexe plusieurs positions au début du contig pour permettre un chevauchement partiel
            for (size_t offset = 0; offset < len && offset < MAX_SEED_DEPTH; offset += SEED_STRIDE) {
                if (offset + index_k > len) break;
                seed_map[getKmerKey(contigs[i], offset, index_k)].push_back({i, false, offset});
                BitVector rc = getBitVectorReverseComplement(contigs[i]);
                seed_map[getKmerKey(rc, offset, index_k)].push_back({i, true, offset});
            }
        }

        // 2. Définition de la logique d'extension (Lambda)
        auto try_extend = [&](BitVector& master) -> bool {
            size_t m_len = master.size() / 2;
            if (m_len < (size_t)index_k) return false;

            // On scanne la FIN du master pour voir si elle correspond à un début indexé
            size_t scan_limit = (m_len > MAX_SCAN_DEPTH) ? MAX_SCAN_DEPTH : (m_len - index_k);

            for (size_t offset = 0; offset < scan_limit; ++offset) {
                size_t probe_pos = m_len - index_k - offset; // On recule depuis la fin
                uint64_t key = getKmerKey(master, probe_pos, index_k);

                if (seed_map.find(key) == seed_map.end()) continue;

                for (const auto& seed : seed_map[key]) {
                    size_t j = seed.id;
                    if (absorbed[j]) continue;

                    // Vérifie que l'alignement géométrique est cohérent
                    if (probe_pos < seed.offset_in_cand) continue;

                    size_t align_start = probe_pos - seed.offset_in_cand;

                    // Prépare le candidat (Inverse si nécessaire)
                    const BitVector& candidate_ref = contigs[j];
                    BitVector candRC;
                    const BitVector* to_check = &candidate_ref;
                    if (seed.is_rc) {
                        candRC = getBitVectorReverseComplement(candidate_ref);
                        to_check = &candRC;
                    }

                    // Vérification du chevauchement (Overlap verification)
                    size_t master_rem = m_len - align_start;
                    size_t check_len = std::min(master_rem, to_check->size()/2);
                    size_t max_err = std::max((size_t)1, (size_t)(check_len * overlap_error_percent));
                    size_t mismatches = 0;
                    bool match = true;

                    for(size_t k = 0; k < check_len; ++k) {
                         if (master.getNucleotideAt(align_start + k) != to_check->getNucleotideAt(k)) {
                             mismatches++;
                             if (mismatches > max_err) { match = false; break; }
                         }
                    }

                    // Si match validé et que le candidat apporte de la nouvelle info (plus long)
                    if (match && to_check->size()/2 > master_rem) {
                         // On coupe le master juste avant le début du candidat et on colle le candidat
                         master.resize(align_start * 2);
                         master.append(*to_check);
                         absorbed[j] = true; // Le candidat est absorbé
                         return true; // Succès
                    }
                }
            }
            return false;
        };

        // 3. Boucle principale de tentative d'extension
        for (size_t i = 0; i < contigs.size(); ++i) {
            if (absorbed[i]) continue;

            // Essai 1: Étendre le Master tel quel
            if (try_extend(contigs[i])) {
                global_change = true;
            }
            // Essai 2: Étendre le Master inversé (si le contig a été construit à l'envers)
            else {
                BitVector masterRC = getBitVectorReverseComplement(contigs[i]);
                if (try_extend(masterRC)) {
                    contigs[i] = masterRC; // Met à jour le Master avec la version retournée + étendue
                    global_change = true;
                }
            }
        }
    }

    // Construction du vecteur final    
    for (size_t k = 0; k < contigs.size(); /* rien */) {
      if (absorbed[k]) {
          contigs.erase(contigs.begin() + k);
          absorbed.erase(absorbed.begin() + k);  // ← important : maintenir la synchro
      } else {
          ++k;
      }
    }
}

/**
 * @brief Reconstruit les nucléotides à partir de l'entier uint64 et les ajoute au BitVector.
 */
void GraphDBJ::addKmerToBitVector(BitVector& bv, uint64_t val, int length) const {
    // On itère de haut en bas car les bits de poids fort contiennent les premiers nucléotides
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);

        // Mapping bits -> bools pour BitVector::push_back
        // 0 (00), 1 (01), 2 (10), 3 (11)
        bv.push_back((two_bits & 2) != 0);
        bv.push_back((two_bits & 1) != 0);
    }
}

/**
 * @brief Génère les contigs en parcourant les chemins simples du graphe.
 * Un chemin simple est une suite de noeuds (u,v) où u n'a que v comme enfant et v n'a que u comme parent.
 */
std::vector<BitVector> GraphDBJ::generateContigs() const {
    std::vector<BitVector> contigs;
    std::unordered_map<Noeud*, bool> visited; // Pour ne pas traiter deux fois le même noeud
    const double COVERAGE_RATIO = config.COVERAGE_RATIO;
    const size_t MAX_CONTIG_LEN = config.MAX_CONTIG_LEN;

    // On parcourt tous les noeuds possibles comme point de départ
    for (const auto& pair : nodes_map) {
        Noeud* startNode = pair.second;

        if (startNode->removed || visited[startNode]) continue;

        // Un début de contig est souvent un noeud sans parents ou une bifurcation (plusieurs parents)
        if (startNode->parents.empty() || startNode->parents.size() > 1) {

            // --- INITIALISATION DU NOUVEAU CONTIG ---
            BitVector currentContig;
            // Le premier noeud donne k-1 nucléotides
            addKmerToBitVector(currentContig, startNode->p, k - 1);

            visited[startNode] = true;
            Noeud* curr = startNode;
            size_t sanity_check = 0;

            // --- PARCOURS (TRAVERSAL) ---
            while (sanity_check < MAX_CONTIG_LEN) {
                sanity_check++;
                Noeud* next = nullptr;

                // Cas simple : chemin linéaire (1 enfant)
                if (curr->c.size() == 1) {
                    next = curr->c[0];
                }
                // Cas complexe : bifurcation (plusieurs enfants)
                // On essaie de résoudre par la couverture (heuristic greedy)
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
                    // On ne suit le chemin que s'il est clairement dominant (selon COVERAGE_RATIO)
                    if (best_candidate != nullptr && max_cov >= (second_max_cov * COVERAGE_RATIO)) {
                        next = best_candidate;
                    } else {
                        break; // Ambiguïté trop forte, on coupe le contig ici
                    }
                } else {
                    break; // Cul-de-sac
                }

                if (next == nullptr || next->removed) break;

                // --- EXTENSION DU CONTIG ---
                // On ajoute seulement le DERNIER nucléotide du k-1 mer suivant
                // car les k-2 précédents chevauchent.
                uint64_t val = next->p;
                uint64_t last_nuc_val = val & 3ULL; // Masque binaire (11) pour les 2 derniers bits

                currentContig.push_back((last_nuc_val & 2) != 0);
                currentContig.push_back((last_nuc_val & 1) != 0);

                if (visited[next]) break; // Cycle détecté

                visited[next] = true;
                curr = next;
            }
            contigs.push_back(currentContig);
        }
    }
    return contigs;
}

/**
 * @brief Constructeur principal : Construit le graphe de De Bruijn.
 * @details Parcourt toutes les lectures, extrait tous les k-mers glissants, et crée les noeuds/arêtes.
 * Gère le fait que l'ADN peut être lu dans les deux sens (Canonical k-mers).
 */
GraphDBJ::GraphDBJ(const Convert& converter, int kmer_size, const GraphDBJConfig& conf) : k(kmer_size), config(conf) {
    if (k < 2 || k > 32) throw std::invalid_argument("K invalide (2-32)");

    const BitVector& bv = converter.getBitVector();
    const std::vector<size_t>& read_ends = converter.getEndPos();
    size_t current_read_start = 0;

    // Itération sur chaque lecture (read)
    for (size_t end_pos : read_ends) {
        size_t read_len_nuc = (end_pos - current_read_start) / 2;

        if (read_len_nuc >= (size_t)k) {
            // Fenêtre glissante sur la lecture
            for (size_t i = 0; i <= read_len_nuc - k; ++i) {
                // Extraction des valeurs brutes
                // u_fwd : k-1 mer préfixe
                // v_fwd : k-1 mer suffixe (le prochain noeud)
                uint64_t u_fwd = extractKmerValue(bv, current_read_start + i * 2, k - 1);
                uint64_t v_fwd = extractKmerValue(bv, current_read_start + i * 2 + 2, k - 1); // +2 bits = +1 nucléotide

                // Calcul des compléments inverses pour le graphe bidirectionnel
                uint64_t v_rev = getReverseComplement(v_fwd, k - 1);
                uint64_t u_rev = getReverseComplement(u_fwd, k - 1);

                // Fonction locale pour créer/lier les noeuds
                auto add_edge = [&](uint64_t source_val, uint64_t target_val) {
                    Noeud*& sourceNode = nodes_map[source_val];
                    if (!sourceNode) sourceNode = new Noeud(source_val);

                    Noeud*& targetNode = nodes_map[target_val];
                    if (!targetNode) targetNode = new Noeud(target_val);

                    // Mise à jour stats
                    sourceNode->coverage++;
                    targetNode->coverage++;

                    // Création du lien s'il n'existe pas déjà
                    bool edgeExists = false;
                    for (auto* child : sourceNode->c) { if (child == targetNode) { edgeExists = true; break; } }
                    if (!edgeExists) {
                        sourceNode->c.push_back(targetNode);
                        targetNode->parents.push_back(sourceNode);
                    }
                };

                // Ajout de l'arête Brin Forward : u -> v
                add_edge(u_fwd, v_fwd);
                // Ajout de l'arête Brin Reverse : RC(v) -> RC(u)
                // C'est crucial pour que le graphe soit navigable quel que soit le sens de lecture
                add_edge(v_rev, u_rev);
            }
        }
        current_read_start = end_pos;
    }
}

GraphDBJ::~GraphDBJ() {
    for (auto& pair : nodes_map) delete pair.second;
    nodes_map.clear();
}

std::vector<Noeud*> GraphDBJ::getNodes() const {
    std::vector<Noeud*> result;
    for (const auto& pair : nodes_map) result.push_back(pair.second);
    return result;
}

uint64_t GraphDBJ::extractKmerValue(const BitVector& bv, size_t start_bit_idx, int len_nucleotides) const {
    uint64_t val = 0;
    for (int i = 0; i < len_nucleotides; ++i) {
        size_t pos = start_bit_idx + (i * 2);
        bool b1 = bv.test(pos);
        bool b2 = bv.test(pos + 1);

        val = val << 2; // Décale pour faire place au nouveau nucléotide
        if (b1) val |= 2ULL; // Bit fort
        if (b2) val |= 1ULL; // Bit faible
    }
    return val;
}

uint64_t GraphDBJ::getReverseComplement(uint64_t val, int length) const {
    uint64_t res = 0;
    for (int i = 0; i < length; ++i) {
        uint64_t nuc = val & 3ULL;   // Prend les 2 derniers bits
        uint64_t comp = nuc ^ 3ULL;  // Complément (A<->T, C<->G)
        // Place le complément à la position opposée
        res |= (comp << (2 * (length - 1 - i)));
        val >>= 2; // Passe au nucléotide suivant
    }
    return res;
}

void GraphDBJ::disconnectNodes(Noeud* parent, Noeud* child) {
    if (!parent || !child) return;
    // Idiome "Erase-Remove" pour supprimer un élément d'un vecteur
    parent->c.erase(std::remove(parent->c.begin(), parent->c.end(), child), parent->c.end());
    child->parents.erase(std::remove(child->parents.begin(), child->parents.end(), parent), child->parents.end());
}

/**
 * @brief Simplification : Coupe les "Tips" (Culs-de-sac).
 * @details Un Tip est une chaîne courte qui se termine brusquement.
 * Souvent dû à une erreur de séquençage en fin de lecture.
 */
int GraphDBJ::clipTips() {
    int tips_removed_count = 0;
    bool changed = true;

    const size_t TOPO_MAX_LEN = config.TOPO_MAX_LEN;
    const size_t RCTC_MAX_LEN = config.RCTC_MAX_LEN;
    const double RCTC_RATIO = config.RCTC_RATIO;

    // Boucle tant qu'on trouve des tips à supprimer
    while (changed) {
        changed = false;
        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);

        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue; // Sécurité suppression concurrente
            Noeud* node = nodes_map[key];
            if (node->removed) continue;

            // Condition 1 : C'est une fin de chaîne (pas d'enfants)
            if (node->c.empty()) {

                // Remonter le chemin en arrière pour voir la longueur du tip
                std::vector<Noeud*> tip_path;
                Noeud* curr = node;
                bool is_valid_tip = true;
                uint64_t sum_coverage = 0;

                while (tip_path.size() <= RCTC_MAX_LEN) {
                    tip_path.push_back(curr);
                    sum_coverage += curr->coverage;

                    if (curr->parents.empty()) { is_valid_tip = false; break; } // Tip flottant (isolé), pas attaché au graphe

                    Noeud* parent = curr->parents[0];

                    // Si le parent a plusieurs enfants, c'est le point d'embranchement (l'ancrage)
                    if (parent->c.size() > 1) break;

                    // Si le parent a plusieurs parents, structure complexe, on abandonne
                    if (parent->parents.size() != 1) { is_valid_tip = false; break; }

                    curr = parent;
                }

                if (!is_valid_tip || tip_path.empty()) continue;

                Noeud* last_tip_node = tip_path.back();
                if (last_tip_node->parents.empty()) continue;

                Noeud* anchor = last_tip_node->parents[0]; // Le noeud principal du graphe
                int len = tip_path.size();

                double tip_avg_cov = (double)sum_coverage / len;
                double anchor_cov = (double)anchor->coverage;
                bool should_remove = false;

                // Critère 1 : Très court (Topologique)
                if (len <= TOPO_MAX_LEN) should_remove = true;
                // Critère 2 : Court et couverture faible par rapport à l'ancrage
                else if (len <= RCTC_MAX_LEN && anchor_cov > (tip_avg_cov * RCTC_RATIO)) should_remove = true;

                if (should_remove) {
                    disconnectNodes(anchor, last_tip_node); // Coupe le lien avec le graphe
                    for (Noeud* n : tip_path) n->removed = true; // Marque tout le bras comme supprimé
                    tips_removed_count++;
                    changed = true;
                }
            }
        }
    }
    return tips_removed_count;
}

/**
 * @brief Simplification : Résolution des bulles.
 * @details Une bulle est une divergence (SNP ou erreur indel) qui converge plus loin.
 * A -> B -> D
 * A -> C -> D
 * On garde le chemin le plus couvert (probablement le bon).
 */
int GraphDBJ::resolveBubbles() {
    int search_depth = (int)(k * config.SEARCH_DEPTH_FACTOR);
    int max_passes = config.MAX_PASSES;

    int bubbles_collapsed = 0;
    bool changed = true;
    int pass = 0;

    while (changed && pass < max_passes) {
        changed = false; pass++;
        std::vector<uint64_t> keys;
        keys.reserve(nodes_map.size());
        for(auto& p : nodes_map) keys.push_back(p.first);

        for (uint64_t key : keys) {
            if (nodes_map.find(key) == nodes_map.end()) continue;
            Noeud* s = nodes_map[key];

            // On cherche une bifurcation simple (2 enfants)
            if (s->removed || s->c.size() != 2) continue;

            Noeud* branch1 = s->c[0];
            Noeud* branch2 = s->c[1];

            std::vector<Noeud*> path1, path2;
            path1.push_back(branch1);
            path2.push_back(branch2);

            Noeud* convergence = nullptr;

            // Exploration en largeur limitée (BFS-like) sur les deux branches
            for (int step = 0; step < search_depth; ++step) {
                Noeud* last1 = path1.back();
                // Avance sur path1 si linéaire
                if (last1->c.size() == 1 && !last1->c[0]->removed) path1.push_back(last1->c[0]);

                Noeud* last2 = path2.back();
                // Avance sur path2 si linéaire
                if (last2->c.size() == 1 && !last2->c[0]->removed) path2.push_back(last2->c[0]);

                // Test croisé : est-ce que path1 a touché un noeud de path2 (ou vice versa) ?
                for(Noeud* n2 : path2) { if (path1.back() == n2) { convergence = n2; break; } }
                if(convergence) break;

                for(Noeud* n1 : path1) { if (path2.back() == n1) { convergence = n1; break; } }
                if(convergence) break;

                // Si on ne peut plus avancer, on arrête
                if (path1.back() == last1 && path2.back() == last2) break;
            }

            // Si convergence trouvée (et ce n'est pas une boucle immédiate sur soi-même)
            if (convergence && convergence != branch1 && convergence != branch2) {
                // Calcul de la couverture totale de chaque branche
                uint32_t cov1 = 0; for(auto* n : path1) { if(n==convergence) break; cov1 += n->coverage; }
                uint32_t cov2 = 0; for(auto* n : path2) { if(n==convergence) break; cov2 += n->coverage; }

                // Décision : qui supprimer ?
                std::vector<Noeud*>& to_remove = (cov1 >= cov2) ? path2 : path1;
                Noeud* branch_to_keep = (cov1 >= cov2) ? branch1 : branch2;

                // Suppression logique
                for(auto* n : to_remove) { if (n != convergence) n->removed = true; }

                // Recâblage : Le noeud S pointe maintenant uniquement vers la branche gagnante
                Noeud* branch_to_remove = (cov1 >= cov2) ? branch2 : branch1;
                disconnectNodes(s, branch_to_remove);

                // Vérifie que S est bien connecté à la branche gardée
                bool connected = false;
                for(auto* child : s->c) if(child == branch_to_keep) connected = true;
                if(!connected) s->c.push_back(branch_to_keep);

                bubbles_collapsed++;
                changed = true;
            }
        }
    }
    return bubbles_collapsed;
}

std::string GraphDBJ::kmerToString(uint64_t val, int length) const {
    std::string seq = "";
    for (int i = length - 1; i >= 0; --i) {
        uint64_t mask = 3ULL << (i * 2);
        uint64_t two_bits = (val & mask) >> (i * 2);
        if (two_bits == 0) seq += 'A';
        else if (two_bits == 2) seq += 'C';
        else if (two_bits == 1) seq += 'G';
        else seq += 'T';
    }
    return seq;
}

void GraphDBJ::exportToGFA(const std::string& filename) const {
    std::ofstream out(filename);
    if (!out.is_open()) return;

    // En-tête GFA
    out << "H\tVN:Z:1.0\n";

    // Mapping ID (uint64) -> ID GFA (int lisible 1, 2, 3...)
    std::unordered_map<uint64_t, int> id_map;
    int current_id = 1;

    // Écriture des Segments (S) = Noeuds
    for (const auto& pair : nodes_map) {
        if (!pair.second->removed) id_map[pair.first] = current_id++;
    }
    for (const auto& pair : nodes_map) {
        Noeud* node = pair.second;
        if (node->removed) continue;
        std::string seq = kmerToString(node->p, k - 1);
        // Format: S <id> <seq> DP:i:<coverage>
        out << "S\t" << id_map[node->p] << "\t" << seq << "\tDP:i:" << node->coverage << "\n";
    }

    // Écriture des Liens (L) = Arêtes
    for (const auto& pair : nodes_map) {
        Noeud* src = pair.second;
        if (src->removed) continue;
        for (Noeud* dest : src->c) {
            if (dest->removed) continue;
            if (id_map.find(src->p) != id_map.end() && id_map.find(dest->p) != id_map.end()) {
                // Le CIGAR string est (k-2)M car le chevauchement est de k-2 nucléotides
                out << "L\t" << id_map[src->p] << "\t+\t" << id_map[dest->p] << "\t+\t" << (k - 2) << "M\n";
            }
        }
    }
    out.close();
}
