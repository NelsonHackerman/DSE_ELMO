import matplotlib.pyplot as plt
import numpy as np

phi = (1+np.sqrt(5))/2
ratio = phi
norbits = 36

angle = 2*np.pi/ratio

plt.figure()  # Create a new figure
ax = plt.subplot(111, polar=True)  # Add a polar subplot

for i in range(36):
    r = np.linspace(0,2,10)
    theta1 = np.ones(10)*i*angle
    #theta2 = np.ones(10)*i*(angle)+np.pi
    ax.plot(theta1, r, color='blue')
    #ax.plot(theta2, r, color='blue')


# Customize the plot
ax.set_title("Custom Polar Plot", va='bottom')  # Set title
ax.set_rticks([0.5, 1.0, 1.5, 2.0])  # Radial ticks
ax.grid(True)  # Enable grid lines

# Add a legend
ax.legend(loc='upper right')
plt.savefig("Mission operations/Polar_ground_track")
