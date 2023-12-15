import matplotlib.pyplot as plt

# Read best fitness values from the text file
with open("best_fitness_values.txt", "r") as file:
    best_fitness_values = [float(line) for line in file]

# Plot the best fitness values
plt.figure(figsize=(10, 5))
plt.plot(best_fitness_values, marker='X', linestyle='-', color = 'g')
plt.title("Best Fitness Values Over Generations")
plt.xlabel("Generation")
plt.ylabel("Best Fitness Value")
plt.grid(True)

# Save the plot as an image or display it
plt.savefig("best_fitness_plot.png")
plt.show()  # Uncomment if you want to display the plot in real-time

print("Best fitness plot saved as 'best_fitness_plot.png'")
