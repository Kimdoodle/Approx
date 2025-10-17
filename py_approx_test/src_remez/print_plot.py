import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from .basic import evalP, error_func

def plot_coeff(coeff: list, a: float, b: float):
    x = np.linspace(a, b, 10000)
    y = evalP(coeff, x)
    plt.plot(x, y, label='p(x)', color='b')

def plot_F(evalF, intervals: list):
    err = 1e-10
    for i, interval in enumerate(intervals):
        a, b = interval
        length = b - a
        n_points = max(200, min(5000, int(1000 * (1/length))))
        x = np.linspace(a + err, b - err, n_points)
        y = np.array([evalF(xi) for xi in x])
        lw = 1.0 if length > 1e-3 else 2.5
        
        plt.plot(x, y, 
                 label='f(x)' if i == 0 else None, 
                 color='r', linewidth=lw)
        if length < 1e-4:
            plt.scatter(x, y, color='red', s=8)

        # x = np.linspace(interval[0]+err, interval[1]-err, 1000)
        # y = np.array([evalF(xi) for xi in x])
        # plt.plot(x, y, label='f(x)' if i == 0 else None, color='r')

def plot_error(evalF, coeff, intervals: list):
    err = 1e-10
    for i, interval in enumerate(intervals):
        x = np.linspace(interval[0]+err, interval[1]-err, 1000)
        y = np.array([error_func(coeff, xi, evalF(xi)) for xi in x])
        plt.plot(x, y, label='error x' if i == 0 else None, color='green')

def plot_max_points(x: list, y: list):
    plt.scatter(x, y, label='local max err')

def remove_images():
    image_dir = Path("./image")
    if not image_dir.exists():
        image_dir.mkdir(parents=True, exist_ok=True)
        return
    for png_file in image_dir.glob("*.png"):
        png_file.unlink()

def draw_graph(n:int, iter: int, err: float, image_mode: str):
    plt.grid(True)
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    if image_mode == "print":
        plt.show()
    else:
        title = f"n={n}, iteration {iter}, err={err}"
        plt.title(title)
        plt.savefig(f"./image/{n}_{iter}.png", dpi=300)

def draw(coeff, evalF, a, intervals, iter, pre, max_err):
    plt.clf()
    plot_coeff(coeff, -a, a)
    plot_F(evalF, intervals)
    for start, end in intervals:
        y1 = evalP(coeff, start)
        y2 = evalP(coeff, end)
        plt.text(start, y1+0.2, f"({start}, {y1})", ha='center', color="black", fontsize=5, fontweight='bold')
        plt.text(end, y2-0.1, f"({end}, {y2})", ha='center', color="black", fontsize=5, fontweight='bold')
    plt.grid(True)
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("y")
    title = f"iteration={iter}, precision={pre}, max_err={max_err}"
    plt.title(title)
    plt.savefig(f"./image/{iter}_{pre}.png", dpi=300)