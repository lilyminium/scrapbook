PCM_KEYWORDS = {
    "pcm": "true",
    "pcm_scf_type": "total",
    "pcm__input": r"""
        Units = Angstrom
        Medium {
            SolverType = CPCM
            Solvent = water
        }

        Cavity {
            RadiiSet = Bondi # Bondi | UFF | Allinger
            Type = GePol
            Scaling = True # radii for spheres scaled by 1.2
            Area = 0.3
            Mode = Implicit
        }
        """
}

PROTOCOLS = {"wavefunction": "orbitals_and_eigenvalues"}