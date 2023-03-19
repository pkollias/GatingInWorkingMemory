import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def main():
    ax = plt.axes(projection='3d')
    ax.plot3D([1, 2, 3], [1, 2, 3], [1, 3, 2])
    plt.show()

main()
