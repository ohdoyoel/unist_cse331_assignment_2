#!/bin/bash

# Create results directory if it doesn't exist
mkdir -p result

# Define algorithms and datasets
ALGORITHMS=("christofides_edmonds" "christofides_gabow" "mst2approx" "heldkarp")
DATASETS=("ulysses16" "ulysses22" "a280" "xql662" "kz9976")

# Create or clear the results file
echo "dataset,algorithm,run,time,cost,optimal_cost,approximation_ratio,gap" > result/experiment_results.csv

# Function to run an algorithm on a dataset
run_algorithm() {
    local algo=$1
    local dataset=$2
    local run_number=$3

    # Run the algorithm and capture its output
    local output=$(./src/$algo data/${dataset}.tsp data/${dataset}.tour 2>&1)
    local exit_code=$?
    
    if [ $exit_code -eq 0 ] && [ ! -z "$output" ]; then
        # Skip if this is a warning message from Held-Karp
        if echo "$output" | grep -q "Warning: Held-Karp is impractical"; then
            echo "Skipping Held-Karp for large instance" >&2
            return 1
        fi
        
        # Parse the output line by line
        local time=$(echo -n "$output" | grep "Elapsed time:" | head -n1 | tr -d '\r' | sed 's/.*: \(.*\) sec/\1/')
        local cost=$(echo -n "$output" | grep "^Tour cost:" | head -n1 | tr -d '\r' | sed 's/.*: \([0-9.e+-]*\)/\1/')
        local opt_cost=$(echo -n "$output" | grep "Optimal Tour cost:" | head -n1 | tr -d '\r' | sed 's/.*: \([0-9.e+-]*\)/\1/')
        local ratio=$(echo -n "$output" | grep "Approximation Ratio:" | head -n1 | tr -d '\r' | sed 's/.*: \([0-9.e+-]*\)/\1/')
        local gap=$(echo -n "$output" | grep "Gap:" | head -n1 | tr -d '\r' | sed 's/.*: \([0-9]*\)/\1/')
        
        # Debug output
        echo "Debug: Command output:" >&2
        echo "$output" | cat -A >&2
        echo "Debug: Exit code: $exit_code" >&2
        echo "Debug: Parsed values:" >&2
        echo "Time: $time" >&2
        echo "Cost: $cost" >&2
        echo "Optimal cost: $opt_cost" >&2
        echo "Ratio: $ratio" >&2
        echo "Gap: $gap" >&2
        
        # Verify all values are present and numeric (including scientific notation)
        if [[ $time =~ ^[0-9.e+-]+$ ]] && [[ $cost =~ ^[0-9.e+-]+$ ]] && \
           [[ $opt_cost =~ ^[0-9.e+-]+$ ]] && [[ $ratio =~ ^[0-9.e+-]+$ ]] && \
           [[ $gap =~ ^[0-9]+$ ]]; then
            # Combine into a single line
            echo "${dataset},${algo},${run_number},${time},${cost},${opt_cost},${ratio},${gap}"
        else
            echo "Error: Invalid output format from ${algo}" >&2
            return 1
        fi
    else
        echo "Error running ${algo} on ${dataset} (run ${run_number})" >&2
        return 1
    fi
}

# Run experiments
for dataset in "${DATASETS[@]}"; do
    echo "Processing dataset: ${dataset}"
    for algo in "${ALGORITHMS[@]}"; do
        echo "  Running algorithm: ${algo}"
        result=$(run_algorithm $algo $dataset 1)
        if [ $? -eq 0 ]; then
            echo "$result" >> result/experiment_results.csv
        fi
    done
done

echo "Experiments completed. Results are in result/experiment_results.csv" 