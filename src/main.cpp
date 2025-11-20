#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <chrono>
#include <cstring>
#include <iomanip>

#include "graphdbj.h"
#include "convert.h"
#include "compare.h"

// Variable globale pour le mode debug
bool DEBUG_MODE = false;

// --- Fonction utilitaire pour mesurer le temps ---
template<typename TimeT = std::chrono::milliseconds>
struct Timer
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
    std::string name;

    Timer(const std::string& taskName) : name(taskName) {
        start = std::chrono::high_resolution_clock::now();
        if (DEBUG_MODE) {
            std::cout << "\n[DEBUG] DEBUT : " << name << std::endl;
        }
    }

    ~Timer() {
        if (DEBUG_MODE) {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<TimeT>(end - start);
            std::cout << "[DEBUG] FIN : " << name << " en " << duration.count() << "ms" << std::endl;
        }
    }
};

// --- Affichage de l'aide ---
void print_usage(const char* prog_name, const GraphDBJConfig& def_conf) {
    std::cout << "Usage: " << prog_name << " <input.fasta> <output_prefix> [OPTIONS]\n\n"
              << "Arguments obligatoires:\n"
              << "  <input.fasta>        Fichier contenant les lectures (FASTA)\n"
              << "  <output_prefix>      Prefixe pour les fichiers de sortie\n\n"
              << "Options Generales:\n"
              << "  -k <int>             Taille des k-mers (defaut: 31)\n"
              << "  --fuse               Activer l'etape de fusion des contigs (defaut: inactif)\n"
              << "  --gfa                Exporter le graphe au format GFA (defaut: non)\n"
              << "  --debug              Afficher les temps d'execution et infos detailles\n"
              << "  --help, -h           Afficher ce message\n\n"
              << "Options de l'Assembleur (GraphDBJ):\n"
              << "  --simplification-passes <int>       Nb max de passes de simplification (defaut: " << def_conf.MAX_PASSES << ")\n\n"
              << "  --popping-passes <int>        Nb max de passes de suppression de tips et bulle (defaut: 10" << ")\n\n"
              << "  --overlap-err <dbl>  % d'erreur autorise pour chevauchement (defaut: " << def_conf.ERROR_PERCENT_OVERLAP << ")\n"
              << "  --contained-err <dbl>% d'erreur autorise pour inclusion (defaut: " << def_conf.ERROR_PERCENT_CONTAINED << ")\n\n"
              << "  --cov-ratio <dbl>    Ratio de couverture pour bifurcations (defaut: " << def_conf.COVERAGE_RATIO << ")\n"
              << "  --tip-ratio <dbl>    Ratio couverture ancrage/bout (defaut: " << def_conf.RCTC_RATIO << ")\n\n"
              << "  --search-depth <dbl> Facteur de profondeur de recherche (defaut: " << def_conf.SEARCH_DEPTH_FACTOR << ")\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    // 1. Gestion des arguments
    if (argc < 3) {
        // On cree une config bidon juste pour afficher les valeurs par defaut dans l'aide
        print_usage(argv[0], GraphDBJConfig(31));
        return 1;
    }

    std::string filename = argv[1];
    std::string output_prefix = argv[2];

    // Paramètres par défaut
    int k_size = 31;
    bool export_gfa = false;
    bool enable_fusion = false; // Par défaut, la fusion est désactivée

    // On stocke les paramètres optionnels pour les appliquer plus tard
    // Valeurs par défaut initiales (seront écrasées par GraphDBJConfig mais utiles pour le parsing)
    int max_passes = -1;
    int max_passes_pop = 10;
    double overlap_err = -1.0;
    double contained_err = -1.0;
    double cov_ratio = -1.0;
    double tip_ratio = -1.0;
    double search_depth = -1.0;

    // Parsing manuel des arguments restants
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-k") {
            if (i + 1 < argc) k_size = std::stoi(argv[++i]);
            else { std::cerr << "Erreur: -k necessite une valeur." << std::endl; return 1; }
        }
        else if (arg == "--fuse") {
            enable_fusion = true;
        }
        else if (arg == "--gfa") {
            export_gfa = true;
        }
        else if (arg == "--debug") {
            DEBUG_MODE = true;
        }
        else if (arg == "--simplification-passes") {
            if (i + 1 < argc) max_passes = std::stoi(argv[++i]);
        }
        else if (arg == "--popping-passes") {
            if (i + 1 < argc) max_passes_pop = std::stoi(argv[++i]);
        }
        else if (arg == "--overlap-err") {
            if (i + 1 < argc) overlap_err = std::stod(argv[++i]);
        }
        else if (arg == "--contained-err") {
            if (i + 1 < argc) contained_err = std::stod(argv[++i]);
        }
        else if (arg == "--cov-ratio") {
            if (i + 1 < argc) cov_ratio = std::stod(argv[++i]);
        }
        else if (arg == "--tip-ratio") {
            if (i + 1 < argc) tip_ratio = std::stod(argv[++i]);
        }
        else if (arg == "--search-depth") {
            if (i + 1 < argc) search_depth = std::stod(argv[++i]);
        }
        else if (arg == "--help" || arg == "-h") {
            print_usage(argv[0], GraphDBJConfig(31));
            return 0;
        }
        else {
            std::cerr << "Argument inconnu: " << arg << std::endl;
            return 1;
        }
    }

    if (k_size < 2 || k_size > 32) {
        std::cerr << "Erreur: K doit etre entre 2 et 32." << std::endl;
        return 1;
    }

    // --- Initialisation ---
    std::cout << "=== Assembleur GraphDBJ ===" << std::endl;
    std::cout << "Entree  : " << filename << std::endl;
    std::cout << "Sortie  : " << output_prefix << ".*" << std::endl;
    std::cout << "K-mer   : " << k_size << std::endl;
    std::cout << "Fusion  : " << (enable_fusion ? "ACTIVEE" : "DESACTIVEE") << std::endl;
    if (DEBUG_MODE) std::cout << "Mode    : DEBUG (Timers actifs)" << std::endl;

    Timer<std::chrono::milliseconds> total_timer("EXECUTION TOTALE");

    // Configuration de l'objet Config
    GraphDBJConfig config(k_size);
    if (max_passes != -1) config.MAX_PASSES = max_passes;
    if (overlap_err != -1.0) config.ERROR_PERCENT_OVERLAP = overlap_err;
    if (contained_err != -1.0) config.ERROR_PERCENT_CONTAINED = contained_err;
    if (cov_ratio != -1.0) config.COVERAGE_RATIO = cov_ratio;
    if (tip_ratio != -1.0) config.RCTC_RATIO = tip_ratio;
    if (search_depth != -1.0) config.SEARCH_DEPTH_FACTOR = search_depth;

    // Affichage de la config si debug
    if (DEBUG_MODE) {
        std::cout << "\n[CONFIG] Passes Max: " << config.MAX_PASSES << "\n"
                  << "[CONFIG] Overlap Error: " << config.ERROR_PERCENT_OVERLAP << "\n"
                  << "[CONFIG] Contained Error: " << config.ERROR_PERCENT_CONTAINED << "\n"
                  << "[CONFIG] Coverage Ratio: " << config.COVERAGE_RATIO << std::endl;
    }

    Convert converter;

    // 2. Lecture FASTA
    try {
        Timer<std::chrono::milliseconds> t("Lecture FASTA");
        converter.processFile(filename);
    } catch (const std::runtime_error& e) {
        std::cerr << "Erreur fatale: " << e.what() << std::endl;
        return 1;
    }

    const auto& read_ends = converter.getEndPos();
    if (read_ends.empty()) {
        std::cerr << "Aucune lecture trouvee dans le fichier." << std::endl;
        return 1;
    }

    // Stats rapides
    {
        CompareKMers comparator(converter.getBitVector(), read_ends, k_size);
        std::cout << "Lectures chargees : " << comparator.get_nReads() << std::endl;
        if (DEBUG_MODE) std::cout << "K-mers theoriques : " << comparator.get_all_nKmers() << std::endl;
    }

    // 3. Construction du Graphe
    GraphDBJ graph(converter, k_size, config);
    std::cout << "Graphe initial construit: " << graph.getNodes().size() << " noeuds." << std::endl;

    // 4. Simplification
    {
        Timer<std::chrono::milliseconds> t("Simplification");
        bool changed = true;
        int iter = 0;
        std::cout << "Simplification en cours..." << std::endl;

        while (changed && iter < max_passes_pop) {
            changed = false; iter++;
            int bubbles = graph.resolveBubbles();
            int tips = graph.clipTips();

            if (bubbles > 0 || tips > 0) changed = true;

            if (DEBUG_MODE && changed) {
                std::cout << "  Passe " << iter << ": " << bubbles << " bulles, " << tips << " pointes." << std::endl;
            }
        }
        std::cout << "Termine en " << iter << " passes." << std::endl;
    }

    // 5. Generation Contigs
    std::vector<BitVector> contigs;
    {
        Timer<std::chrono::milliseconds> t("Generation Contigs");
        contigs = graph.generateContigs();
        std::cout << "Contigs bruts generes : " << contigs.size() << std::endl;
    }

    // 6. Fusion (Overlap/Contained) - Conditionnée par --fuse
    if (enable_fusion) {
        Timer<std::chrono::milliseconds> t("Fusion des Contigs");
        int min_overlap = k_size / 2;
        // On peut rendre min_overlap configurable si besoin, ici hardcodé à k/2

        contigs = GraphDBJ::mergeContigs(
            contigs,
            min_overlap,
            config.ERROR_PERCENT_OVERLAP,
            config.ERROR_PERCENT_CONTAINED
        );
        std::cout << "Contigs finaux apres fusion : " << contigs.size() << std::endl;
    } else {
        std::cout << "Etape de fusion desactivee (utiliser --fuse pour activer)." << std::endl;
    }

    // 7. Export
    {
        Timer<std::chrono::milliseconds> t("Ecriture Fichiers");

        // Export GFA optionnel
        if (export_gfa) {
            std::string gfa_name = output_prefix + ".gfa";
            std::cout << "Export GFA vers " << gfa_name << "..." << std::endl;
            graph.exportToGFA(gfa_name);
        }

        // Export FASTA (Toujours fait)
        std::string fasta_name = output_prefix + ".contigs.fasta";
        std::ofstream out_contigs(fasta_name);
        if (!out_contigs.is_open()) {
             std::cerr << "Erreur: Impossible d'ecrire le fichier " << fasta_name << std::endl;
             return 1;
        }

        int exported_count = 0;
        for (size_t i = 0; i < contigs.size(); ++i) {
            size_t len_bp = contigs[i].size() / 2;
            // Filtre minimal : on ne garde que ce qui est >= 2*k (optionnel, garde-fou)
            if (len_bp >= (size_t)k_size) {
                out_contigs << ">contig_" << i << "_len_" << len_bp << "\n";
                out_contigs << contigs[i].readBitVector() << "\n";
                exported_count++;
            }
        }
        out_contigs.close();
        std::cout << "Export FASTA termine : " << exported_count << " contigs ecrits dans " << fasta_name << std::endl;
    }

    return 0;
}
