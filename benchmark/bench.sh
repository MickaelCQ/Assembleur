#!/bin/bash

# --- Configuration ---
RUNS=10
# Utilise des chemins absolus si possible, ou relatifs safe
DATA_FILE="$(pwd -P)/../data/reads.fasta"
RESULTS_FILE="benchmark_results_raw.tsv"

# Commands to benchmark
# We silence standard output (> /dev/null) to keep the terminal clean
CMD_MINIA="pixi run -e bench minia -in $DATA_FILE -out /tmp/minia_out"
CMD_RAMILASS="pixi run ramilass $DATA_FILE /tmp/ramilass_out --debug --fuse --popping-passes 0 --max-contig-len 11000"
CMD_SPADES="pixi run -e bench spades.py --only-assembler -s $DATA_FILE -o /tmp/spades_out"

# Check for GNU time
TIME_CMD="/usr/bin/time"
if [ ! -f "$TIME_CMD" ]; then
    echo "Error: GNU time not found at $TIME_CMD."
    exit 1
fi

# --- Benchmarking Function ---
run_benchmark() {
    local tool_name="$1"
    local cmd_string="$2"
    local tmp_file="${tool_name}_raw_metrics.txt"

    # >&2 envoie le texte vers le terminal même si on redirige la sortie vers un fichier
    echo "---------------------------------------------------" >&2
    echo "Starting benchmark for: $tool_name ($RUNS runs)" >&2
    echo "${cmd_string}" >&2
    echo "---------------------------------------------------" >&2

    # 1. Exécuter la série de tests
    > "$tmp_file"

    for ((i=1; i<=RUNS; i++)); do
        echo -ne "Run $i/$RUNS... \r" >&2  # La barre de progression va sur l'écran

        # On exécute la commande
        # On redirige stdout et stderr de l'outil vers /dev/null pour ne pas polluer
        $TIME_CMD -f "%e %M" -o "$tmp_file" --append $cmd_string > /dev/null 2>&1
    done

    echo -e "\nDone." >&2

    # 2. Formater les résultats pour le TSV
    # Awk imprime sur stdout (sortie standard), c'est SEULEMENT ça qui ira dans le fichier
    awk -v name="$tool_name" '
    BEGIN {
        i=1;
    }
    {
        # Format: ToolName <tab> RunID <tab> Time <tab> Memory
        printf "%s\t%d\t%s\t%s\n", name, i++, $1, $2;
    }' "$tmp_file"

    # 3. Nettoyer
    rm "$tmp_file"
}

# --- Execution ---

echo "Benchmark initialized. Running $RUNS iterations per tool."

# 1. Initialiser le fichier TSV avec l'en-tête
echo -e "Tool_Name\tRun_ID\tElapsed_Time_s\tMax_RSS_KB" > $RESULTS_FILE

# 2. Run Minia (Les logs vont à l'écran, les données vont dans le fichier)
run_benchmark "Minia" "$CMD_MINIA" >> $RESULTS_FILE

# 3. Run Ramilass
run_benchmark "Ramilass" "$CMD_RAMILASS" >> $RESULTS_FILE

# 4. Run SPAdes
run_benchmark "SPAdes" "$CMD_SPADES" >> $RESULTS_FILE

# Nettoyage final
rm -rf /tmp/spades_out /tmp/minia_out /tmp/ramilass_out > /dev/null 2>&1

echo "---------------------------------------------------"
echo "Benchmark Complete. Raw results saved in $RESULTS_FILE"
echo "---------------------------------------------------"