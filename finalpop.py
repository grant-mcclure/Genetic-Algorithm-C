import numpy as np
import matplotlib.pyplot as plt

# Define the Ackley function
def ackley(x, y):
    a = -20 * np.exp(-0.2 * np.sqrt(0.5 * (x**2 + y**2)))
    b = -np.exp(0.5 * (np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y)))
    return a + b + 20 + np.exp(1)

# Create a mesh grid for the x and y coordinates
x = np.linspace(-5, 5, 400)
y = np.linspace(-5, 5, 400)
X, Y = np.meshgrid(x, y)

# Calculate the Ackley function values on the mesh grid
Z = ackley(X, Y)

# Create a contour plot of the Ackley function
plt.contourf(X, Y, Z, levels=100, cmap='viridis')
plt.colorbar()
plt.title("Contour Plot of Ackley Function")

# Plot the final population as red points
final_population = np.loadtxt("final_population.txt")
variable1_values = final_population[:, 0]
variable2_values = final_population[:, 1]
plt.scatter(variable1_values, variable2_values, c='red', marker='o', label="Final Population")

plt.xlabel('Variable 1')
plt.ylabel('Variable 2')
plt.legend()

plt.show()
