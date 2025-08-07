#!/usr/bin/env python3

import argparse
import sys

def cluster_reads_by_similarity(similarity_file):
    """
    Reads a similarity file and clusters reads into sets.

    Args:
        similarity_file (str): Path to the input file with two columns of read names.

    Returns:
        list: A list of sets, where each set represents a cluster of similar reads.
    """
    read_sets = []
    read_to_set_id = {}

    try:
        with open(similarity_file, 'r') as f:
            for line_number, line in enumerate(f, 1):
                # 跳过空行和以 # 开头的注释行
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                try:
                    read1, read2 = line.split()
                except ValueError:
                    print(f"Skipping malformed line {line_number}: '{line}'", file=sys.stderr)
                    continue

                set_id1 = read_to_set_id.get(read1)
                set_id2 = read_to_set_id.get(read2)

                if set_id1 is None and set_id2 is None:
                    # 如果两个 read 都是新的，创建一个新的集合
                    new_set = {read1, read2}
                    read_sets.append(new_set)
                    new_id = len(read_sets) - 1
                    read_to_set_id[read1] = new_id
                    read_to_set_id[read2] = new_id
                elif set_id1 is not None and set_id2 is None:
                    # 如果 read1 已存在，将 read2 加入其集合
                    read_sets[set_id1].add(read2)
                    read_to_set_id[read2] = set_id1
                elif set_id1 is None and set_id2 is not None:
                    # 如果 read2 已存在，将 read1 加入其集合
                    read_sets[set_id2].add(read1)
                    read_to_set_id[read1] = set_id2
                elif set_id1 != set_id2:
                    # 如果两个 read 都已存在，但属于不同的集合，则将这两个集合合并
                    # 将 set_id2 的所有元素移动到 set_id1，并更新它们的映射
                    set2 = read_sets[set_id2]
                    for read in set2:
                        read_sets[set_id1].add(read)
                        read_to_set_id[read] = set_id1
                    
                    # 清空 set_id2 的集合，避免重复处理
                    read_sets[set_id2] = set()

    except FileNotFoundError:
        print(f"Error: The file '{similarity_file}' was not found.", file=sys.stderr)
        return []

    # 清除合并后为空的集合
    final_read_sets = [s for s in read_sets if s]
    return final_read_sets


def print_results(read_sets):
    """
    Prints the clustered read sets to standard output.
    """
    print(f"Found {len(read_sets)} clusters.")
    for i, read_set in enumerate(read_sets):
        print(f"Cluster {i+1} ({len(read_set)} reads):")
        print(f"  {', '.join(sorted(list(read_set)))}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster reads based on a similarity file.")
    parser.add_argument("input_file", help="Path to the input similarity file.")
    args = parser.parse_args()

    clusters = cluster_reads_by_similarity(args.input_file)
    print_results(clusters)
