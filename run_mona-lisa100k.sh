#!/bin/bash

# Create results directory if it doesn't exist
mkdir -p result

# Define algorithms
ALGORITHMS=("christofides_edmonds" "christofides_gabow" "mst2approx" "heldkarp" "hylos")
DATASET="mona-lisa100k"

# Create or clear the results file
echo "algorithm,run,time,cost" > result/mona-lisa100k_results.csv

# Function to run algorithm on mona-lisa100k
run_algorithm() {
    local algo=$1
    local run_number=$2

    # Run the algorithm and capture its output
    local output=$(./src/$algo data/${DATASET}.tsp data/${DATASET}.tour 2>&1)
    local exit_code=$?
    
    if [ $exit_code -eq 0 ] && [ ! -z "$output" ]; then
        # Skip if this is a warning message from Held-Karp
        if echo -n "$output" | grep -q "Warning: Held-Karp is impractical"; then
            echo "Skipping Held-Karp for large instance" >&2
            return 1
        fi
        
        # Parse the output line by line
        local time=$(echo -n "$output" | grep "Elapsed time:" | head -n1 | tr -d '\r' | sed 's/.*: \(.*\) sec/\1/')
        local cost=$(echo -n "$output" | grep "^Tour cost:" | head -n1 | tr -d '\r' | sed 's/.*Tour cost: \([0-9.e+-]*\).*/\1/')
        
        # Debug output
        echo "Debug: Raw output for ${algo}:" >&2
        echo "$output" >&2
        echo "Debug: Parsed time=$time cost=$cost" >&2
        
        # Verify values are present and numeric (including scientific notation)
        if [[ $time =~ ^[0-9.e+-]+$ ]] && [[ $cost =~ ^[0-9.e+-]+$ ]]; then
            # Combine into a single line
            echo "${algo},${run_number},${time},${cost}"
        else
            echo "Error: Invalid output format from ${algo}" >&2
            echo "Debug: time regex match: ${time}" >&2
            echo "Debug: cost regex match: ${cost}" >&2
            return 1
        fi
    else
        echo "Error running ${algo} on ${DATASET} (run ${run_number})" >&2
        echo "Debug: Algorithm output:" >&2
        echo "$output" >&2
        return 1
    fi
}

# Run experiments
echo "Processing dataset: ${DATASET}"
for algo in "${ALGORITHMS[@]}"; do
    echo "  Running algorithm: ${algo}"
    if [ "$algo" == "hylos" ]; then
        # Run Hylos 10 times
        for run in {1..10}; do
            echo "    Run ${run}/10"
            result=$(run_algorithm $algo $run)
            if [ $? -eq 0 ]; then
                echo "$result" >> result/mona-lisa100k_results.csv
            fi
        done

        # Calculate statistics for Hylos
        echo -e "\nHylos Statistics for ${DATASET}:" > result/mona-lisa100k_hylos_stats.txt
        awk -F',' '
            $1 == "hylos" {
                time[NR] = $3;
                cost[NR] = $4;
                count++;
            }
            END {
                if (count > 0) {
                    # Calculate means
                    time_sum = cost_sum = 0;
                    for (i in time) {
                        time_sum += time[i];
                        cost_sum += cost[i];
                    }
                    time_mean = time_sum/count;
                    cost_mean = cost_sum/count;
                    
                    # Calculate standard deviations
                    time_var = cost_var = 0;
                    for (i in time) {
                        time_var += (time[i] - time_mean)^2;
                        cost_var += (cost[i] - cost_mean)^2;
                    }
                    time_std = sqrt(time_var/(count-1));
                    cost_std = sqrt(cost_var/(count-1));
                    
                    printf("Time_mean=%.6f\n", time_mean);
                    printf("Time_std=%.6f\n", time_std);
                    printf("Cost_mean=%.2f\n", cost_mean);
                    printf("Cost_std=%.2f\n", cost_std);
                } else {
                    print "No_valid_results";
                }
            }
        ' result/mona-lisa100k_results.csv >> result/mona-lisa100k_hylos_stats.txt
    else
        # Run other algorithms once
        result=$(run_algorithm $algo 1)
        if [ $? -eq 0 ]; then
            echo "$result" >> result/mona-lisa100k_results.csv
        fi
    fi
done

echo "Mona Lisa 100K experiments completed. Results are in result/mona-lisa100k_results.csv"
echo "Hylos statistics are in result/mona-lisa100k_hylos_stats.txt" 