#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: $0 <dependency_file.txt>"
    exit 1
fi

DEPENDENCY_FILE=$1
PREV_JOB_IDS=()
LAST_JOB_ID=""
START_DIR=$(pwd)

echo "Starting automated Slurm job submission..."
echo "----------------------------------------"

while read -r line; do
    line=$(echo "$line" | xargs)

#    if [ -z "$line" ]; then
#        echo "Found blank line. Resetting dependency chain."
#        PREV_JOB_IDS=()
#        continue
#    fi
    # Skip comment lines (lines starting with #) and blank lines
    if [[ "$line" =~ ^#.* ]] || [ -z "$line" ]; then
        if [[ "$line" =~ ^#.* ]]; then
            echo "Skipping comment line: $line"
        else
            echo "Found blank line. Resetting dependency chain."
            PREV_JOB_IDS=()
        fi
        continue
    fi

    CURRENT_JOB_IDS=()
    DEPENDENCY_ARG=""
    if [ ${#PREV_JOB_IDS[@]} -gt 0 ]; then
        DEPENDENCY_ARG="--dependency=afterok:$(IFS=:; echo "${PREV_JOB_IDS[*]}")"
    fi

    # The corrected logic for handling wildcards
    # It loops through each word (like './mapping/hifi_*.sh') and expands it
    for pattern in $line; do
        # Use a sub-loop to expand the wildcard
        for script_path in $pattern; do
            if [ ! -f "$script_path" ]; then
                echo "Error: Script not found at path: $script_path" >&2
                echo "The wildcard pattern '$pattern' did not match any files." >&2
                exit 1
            fi
            
            SCRIPT_DIR=$(dirname "$script_path")
            SCRIPT_NAME=$(basename "$script_path")
            
            echo "Submitting script: $SCRIPT_NAME from directory: $SCRIPT_DIR"
            
            cd "$SCRIPT_DIR" || { echo "Error: Failed to change directory to $SCRIPT_DIR" >&2; exit 1; }
            
            JOB_ID=$(sbatch --parsable $DEPENDENCY_ARG "$SCRIPT_NAME")
            
            if [ -z "$JOB_ID" ]; then
                echo "Error: sbatch submission failed for $SCRIPT_NAME" >&2
                cd "$START_DIR"
                exit 1
            fi
            
            echo "  -> Submitted with Job ID: $JOB_ID"
            CURRENT_JOB_IDS+=("$JOB_ID")
            LAST_JOB_ID="$JOB_ID"
            
            cd "$START_DIR"
        done
    done

    PREV_JOB_IDS=("${CURRENT_JOB_IDS[@]}")
    echo "Current stage submitted. Next stage will depend on: ${PREV_JOB_IDS[*]}"
    echo "----------------------------------------"

done < "$DEPENDENCY_FILE"

echo "All jobs submitted. The last submitted job ID is: $LAST_JOB_ID"
