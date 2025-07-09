#!/usr/bin/env python3
"""
Load the wavefunction from output/wavefunction.npz
and animate the probability density |ψ(x,t)|²
in the 1D harmonic oscillator.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def animate_wavefunction(npz_path='output/wavefunction.npz', interval=50):
    """
    Animate the time-evolving probability density |ψ(x,t)|².

    Parameters
    ----------
    npz_path : str
        Path to the .npz file containing 'x', 't', and 'psi'.
    interval : int
        Delay between frames in milliseconds.

    Returns
    -------
    anim : matplotlib.animation.FuncAnimation
        The animation object.
    """
    # Load data
    data = np.load(npz_path)
    x   = data['x']          # shape (N,)
    t   = data['t']          # shape (NT,)
    psi = data['psi']        # shape (NT, N)

    # Compute probability density
    density = np.abs(psi)**2

    # Set up figure
    fig, ax = plt.subplots()
    line, = ax.plot(x, density[0], lw=2)
    # create a text object for time display
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes,
                        ha='left', va='top')

    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(0, density.max() * 1.1)
    ax.set_xlabel('x')
    ax.set_ylabel(r'$|\psi(x,t)|^2$')

    def init():
        line.set_ydata(density[0])
        time_text.set_text(f't = {t[0]:.2f}')
        return line, time_text

    def update(frame):
        line.set_ydata(density[frame])
        time_text.set_text(f't = {t[frame]:.2f}')
        return line, time_text

    # Use blit=False to ensure text updates correctly
    anim = FuncAnimation(
        fig,
        update,
        frames=len(t),
        init_func=init,
        interval=interval,
        blit=False  # disable blitting for title/text
    )

    plt.show()
    return anim


if __name__ == '__main__':
    animate_wavefunction()
