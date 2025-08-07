#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <dependency_file.txt>"
    exit 1
fi

DEPENDENCY_FILE=$1
PREV_JOB_IDS=()
LAST_JOB_ID=""
START_DIR=$(pwd)
WRAPPER_SCRIPT="/path/to/your/job_wrapper.sh" # EDIT THIS: Path to your wrapper script

echo "Starting automated Slurm job submission..."
echo "----------------------------------------"

while read -r line; do
    # Remove leading/trailing whitespace
    line=$(echo "$line" | xargs)

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
        # Construct the dependency argument with all previous job IDs
        DEPENDENCY_ARG="--dependency=afterok:$(IFS=:; echo "${PREV_JOB_IDS[*]}")"
    fi

    # Loop through each word (e.g., './mapping/hifi_*.sh') and expand it
    for pattern in $line; do
        # Use a sub-loop to expand the wildcard and get all matching files
        for script_path in $pattern; do
            if [ ! -f "$script_path" ]; then
                echo "Error: Script not found at path: $script_path" >&2
                echo "The wildcard pattern '$pattern' did not match any files." >&2
                exit 1
            fi

            # --- Key change: Submitting the wrapper script ---
            # We now submit the wrapper script and pass the original script path as an argument.
            echo "Submitting wrapper script for: $script_path"
            
            # The sbatch command now calls the wrapper script, passing the original script path.
            # The wrapper is responsible for cd'ing to the correct directory and executing the script.
            JOB_ID=$(sbatch --parsable $DEPENDENCY_ARG "$WRAPPER_SCRIPT" "$script_path")

            if [ -z "$JOB_ID" ]; then
                echo "Error: sbatch submission failed for wrapper script" >&2
                exit 1
            fi
            
            echo "  -> Submitted with Job ID: $JOB_ID"
            CURRENT_JOB_IDS+=("$JOB_ID")
            LAST_JOB_ID="$JOB_ID"
        done
    done

    # Update the dependencies for the next stage
    PREV_JOB_IDS=("${CURRENT_JOB_IDS[@]}")
    echo "Current stage submitted. Next stage will depend on: ${PREV_JOB_IDS[*]}"
    echo "----------------------------------------"

done < "$DEPENDENCY_FILE"

echo "All jobs submitted. The last submitted job ID is: $LAST_JOB_ID"
