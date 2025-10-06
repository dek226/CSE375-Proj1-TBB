#time ./kmeans_serial < poker_train.data > serial_output.txt
#time ./kmeans_parallel < poker_train.data > parallel_output.txt
#You should see the same clusters (within floating-point tolerance) but a lower execution time for the parallel version, especially on the 23 MB dataset.

# 6. Choosing K and tuning performance
#Start with K = 10 to cluster by the poker hand class count.

#For stress tests, try K = 50, K = 100 — more clusters means more work, so you’ll better see speedup.

#Adjust max_iterations (e.g. 100 or 200) for convergence control.

#python3 convert_poker_to_kmeans.py
#You’ll get:
#poker_train.data


# convert_poker_to_kmeans.py
import csv

input_file = "poker-hand-training-true.data"
output_file = "poker_train.data"

with open(input_file, "r") as fin, open(output_file, "w") as fout:
    reader = csv.reader(fin)
    data = [row for row in reader if row]

    total_points = len(data)
    total_values = 10  # first 10 columns only
    K = 10             # can choose (e.g., 10 clusters)
    max_iterations = 100
    has_name = 0

    # header line
    fout.write(f"{total_points} {total_values} {K} {max_iterations} {has_name}\n")

    for row in data:
        # only write first 10 numeric attributes
        features = row[:10]
        fout.write(" ".join(features) + "\n")
