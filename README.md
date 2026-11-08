# Soufflerie thermique 2D

Simulateur CFD 2D didactique dans le navigateur.
- Écoulement : Lattice Boltzmann D2Q9, collisionneur MRT (Lallemand & Luo) + LES Smagorinsky, forçage de Guo.
- Thermique : D2Q5, Boussinesq ; diffusivité numérique relevée seulement là où le Péclet de maille dépasse 6.
- Échangeurs : zone poreuse (Δp ∝ v²) + source de chaleur (T de batterie ou puissance imposée).
- Ventilateurs : débit imposé par force asservie (régulateur intégral).
- Moteur GPU (WebGPU, `gpu.js`) avec repli automatique sur le processeur (`solver.js`). Les deux donnent les mêmes résultats.

## Utiliser
Ouvrir `index.html` dans un navigateur récent (Chrome/Edge 113+, Safari 26+ pour WebGPU).

## Fichiers
- `solver.js` : noyau processeur (aussi utilisable sous Node).
- `gpu.js` : noyau WebGPU (shaders WGSL), même interface que `solver.js`.
- `index.src.html` : interface, avec les marqueurs `/*SOLVER*/` et `/*GPU*/`.
- `build.py` : régénère `index.html` (`python3 build.py`).
- `tests/` : essais Node (`node tests/test.js`) et comparaison GPU/processeur
  (`python3 -m http.server 8765` à la racine, puis `python3 tests/cmp.py`, nécessite Playwright).

## Limites
Viscosité de base relevée au minimum stable du maillage : affiner la grille rapproche de la physique.
Débits et puissances pour la profondeur réglée dans l'onglet Domaine. 2D : turbulence qualitative.
