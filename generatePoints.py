import csv
import random

points = [
    [0, 0, 0],
    [1, 2, 0],
    [1, 9, 0],
    [2, 5, 0],
    [2, 6, 0],
    [3, 5, 0],
    [3, 6, 0],
    [6, 5, 0],
    [6, 6, 0],
    [7, 1, 0],
    [7, 10, 0]
]

num_points = 1000

output_file = "points.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    writer.writerow([num_points])
    
    for _ in range(num_points):
        random_point = random.choice(points)
        writer.writerow(random_point)

print("done")