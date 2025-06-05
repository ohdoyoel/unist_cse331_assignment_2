#!/bin/bash

# Create results directory if it doesn't exist
mkdir -p result

# Define datasets
DATASETS=("ulysses16" "ulysses22" "a280" "xql662" "kz9976")

# Parse command line arguments
DEBUG_FLAG=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --debug)
            DEBUG_FLAG="--debug"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--debug]"
            exit 1
            ;;
    esac
done

# Create or clear the results file
echo "dataset,algorithm,run,time,cost,optimal_cost,approximation_ratio,gap" > result/hylos_results.csv

# Function to run hylos on a dataset
run_hylos() {
    local dataset=$1
    local run_number=$2

    # Run the algorithm and capture its output
    # 디버그 출력은 stderr로, 결과는 stdout으로 분리
    if [ ! -z "$DEBUG_FLAG" ]; then
        ./src/hylos $DEBUG_FLAG data/${dataset}.tsp data/${dataset}.tour 2>&1 1>/tmp/hylos_output.txt
        cat /tmp/hylos_output.txt > /tmp/hylos_result.txt
    else
        ./src/hylos data/${dataset}.tsp data/${dataset}.tour > /tmp/hylos_result.txt 2>&1
    fi
    
    local exit_code=$?
    local output=$(cat /tmp/hylos_result.txt)
    
    if [ $exit_code -eq 0 ] && [ ! -z "$output" ]; then
        # Parse the output line by line
        local time=$(echo -n "$output" | grep "Elapsed time:" | head -n1 | tr -d '\r' | sed 's/.*: \(.*\) sec/\1/')
        local cost=$(echo -n "$output" | grep "^Tour cost:" | head -n1 | tr -d '\r' | sed 's/.*Tour cost: \([0-9.e+-]*\).*/\1/')
        local optimal_cost=$(echo -n "$output" | grep "Optimal Tour cost:" | head -n1 | tr -d '\r' | sed 's/.*Tour cost: \([0-9.e+-]*\).*/\1/')
        local ratio=$(echo -n "$output" | grep "Approximation Ratio:" | head -n1 | tr -d '\r' | sed 's/.*Ratio: \([0-9.e+-]*\).*/\1/')
        local gap=$(echo -n "$output" | grep "^Gap:" | head -n1 | tr -d '\r' | sed 's/.*Gap: \([0-9]*\).*/\1/')
        
        # Verify all values are present and numeric
        if [[ $time =~ ^[0-9.e+-]+$ ]] && [[ $cost =~ ^[0-9.e+-]+$ ]] && \
           [[ $optimal_cost =~ ^[0-9.e+-]+$ ]] && [[ $ratio =~ ^[0-9.e+-]+$ ]] && \
           [[ $gap =~ ^[0-9]+$ ]]; then
            # Combine into a single line
            echo "${dataset},hylos,${run_number},${time},${cost},${optimal_cost},${ratio},${gap}"
        else
            echo "Error: Invalid output format" >&2
            if [ ! -z "$DEBUG_FLAG" ]; then
                echo "Debug: Raw output:" >&2
                echo "$output" >&2
            fi
            return 1
        fi
    else
        echo "Error running hylos on ${dataset} (run ${run_number})" >&2
        if [ ! -z "$DEBUG_FLAG" ]; then
            echo "Debug: Command output:" >&2
            echo "$output" >&2
        fi
        return 1
    fi
}

# Process each dataset
echo -e "Hylos Statistics:\n" > result/hylos_stats.txt

for dataset in "${DATASETS[@]}"; do
    echo "Processing dataset: ${dataset}"
    
    # Run hylos 10 times
    for run in {1..10}; do
        echo "  Run ${run}/10"
        result=$(run_hylos $dataset $run)
        if [ $? -eq 0 ]; then
            echo "$result" >> result/hylos_results.csv
        fi
    done
    
    # Calculate statistics for this dataset
    echo -e "\n${dataset}:" >> result/hylos_stats.txt
    awk -F',' -v dataset="$dataset" '
        $1 == dataset && $2 == "hylos" {
            times[NR] = $4;
            costs[NR] = $5;
            ratios[NR] = $7;
            gaps[NR] = $8;
            count++;
        }
        END {
            if (count > 0) {
                # Calculate means
                time_sum = cost_sum = ratio_sum = gap_sum = 0;
                for (i in times) {
                    time_sum += times[i];
                    cost_sum += costs[i];
                    ratio_sum += ratios[i];
                    gap_sum += gaps[i];
                }
                time_mean = time_sum/count;
                cost_mean = cost_sum/count;
                ratio_mean = ratio_sum/count;
                gap_mean = gap_sum/count;
                
                # Calculate standard deviations
                time_var = cost_var = ratio_var = gap_var = 0;
                for (i in times) {
                    time_var += (times[i] - time_mean)^2;
                    cost_var += (costs[i] - cost_mean)^2;
                    ratio_var += (ratios[i] - ratio_mean)^2;
                    gap_var += (gaps[i] - gap_mean)^2;
                }
                time_std = sqrt(time_var/(count-1));
                cost_std = sqrt(cost_var/(count-1));
                ratio_std = sqrt(ratio_var/(count-1));
                gap_std = sqrt(gap_var/(count-1));
                
                printf("Time: %.6f ± %.6f sec\n", time_mean, time_std);
                printf("Cost: %.2f ± %.2f\n", cost_mean, cost_std);
                printf("Approximation ratio: %.6f ± %.6f\n", ratio_mean, ratio_std);
                printf("Gap: %.1f ± %.1f\n", gap_mean, gap_std);
            } else {
                print "No valid results";
            }
        }
    ' result/hylos_results.csv >> result/hylos_stats.txt
done

# 임시 파일 정리
rm -f /tmp/hylos_output.txt /tmp/hylos_result.txt

echo "Hylos experiments completed. Results are in result/hylos_results.csv"
echo "Statistics are in result/hylos_stats.txt" 