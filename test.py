import csv

# Data to be written to the CSV file
data = [
    ["real", "imag"],
    [0.5, 0.3],
    [-0.7, 0.8],
    [1.0, -0.5],
    [0.2, 0.6],
    [-0.9, -0.3],
    [0.1, -0.8],
    [-0.4, 0.4],
    [0.3, -0.2],
    [-0.6, 0.7],
    [0.8, -0.1],
]

# Specify the filename
filename = "Filters.csv"

# Write data to the CSV file
with open(filename, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(data)

print(f"CSV file '{filename}' created successfully.")
