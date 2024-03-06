from truncatedhernquist.potential import TruncatedHernquistPotential
import numpy as np


def test_potentialdemo():
    pot = TruncatedHernquistPotential(m=1.,c=1.,rmax=10.)
    E = pot.energy([1., 0.1, -0.41])
    assert np.all(np.isfinite(E))

if __name__ == '__main__':
    test_potentialdemo()
    print("Done.")