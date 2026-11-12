# Soufflerie thermique 2D

Simulateur CFD 2D didactique dans le navigateur.
- Écoulement : Lattice Boltzmann D2Q9, collisionneur MRT (Lallemand & Luo) + LES Smagorinsky, forçage de Guo.
- Thermique : D2Q5, Boussinesq ; diffusivité numérique relevée seulement là où le Péclet de maille dépasse 6.
- Échangeurs : zone poreuse (Δp ∝ v²) + source de chaleur (T de batterie ou puissance imposée).
- Ventilateurs : débit imposé par force asservie (régulateur intégral).
- Moteur GPU (WebGPU, `gpu.js`) avec repli automatique sur le processeur (`solver.js`). Les deux donnent les mêmes résultats.

## Utiliser
Ouvrir `index.html` dans un navigateur récent (Chrome/Edge 113+, Safari 26+ pour WebGPU).

## Cas exemples
Menu « Cas » en haut de l'interface :
- **Pompe à chaleur, unité extérieure** : évaporateur à −3 °C et ventilateur à débit imposé.
- **Radiateur dans une pièce** : boucle de convection naturelle entre radiateur et fenêtre froide.
- **Conduit avec batterie chaude** : convection forcée, puissance imposée, obstacle en aval.
- **Panache thermique dans l'eau** : résistance chauffante au fond d'une cuve.
- **Cellules de Rayleigh-Bénard** : couche d'air entre plaque chaude (35 °C) et plaque froide (5 °C) ;
  un petit défaut sur la plaque chaude amorce l'instabilité et les rouleaux de convection.
- **Allée de Kármán derrière un tube chauffé** : sillage alterné rendu visible par la température.
- **Baie de serveurs et climatiseur** : 3 kW dissipés dans la baie, soufflage froid vers le bas,
  recirculation d'air chaud dans la salle.
- **Cheminée solaire** : mur absorbeur à 55 °C derrière un vitrage, tirage thermique sans ventilateur.
- **Cavité différentiellement chauffée** : paroi verticale chaude face à une paroi froide,
  cas de référence de la convection naturelle (de Vahl Davis).
- **Convecteur sous une fenêtre** : le panache chaud remonte le long du vitrage froid et contre la descente d'air froid.
- **Dissipateur à ailettes** : semelle et ailettes à 70 °C en convection naturelle.
- **Plancher chauffant contre radiateur** : deux pièces identiques côte à côte, chauffées différemment.
- **Rideau d'air de porte** : jet descendant entre l'extérieur (zone maintenue à −5 °C par un échangeur)
  et une pièce chauffée ; mettre le débit du rideau à zéro pour voir l'air froid entrer par le bas.
- **Faisceau de tubes en quinconce** : tubes chauds dans un conduit ventilé, sillages et mélange.
- **Lac en hiver** : eau à 10 °C sous la glace, avec le fluide « Eau froide » dont la densité est maximale
  à 4 °C. L'eau refroidie plonge d'abord, puis l'eau sous 4 °C, plus légère, reste en surface.
- **Rue canyon et façade ensoleillée** : vent de 3 m/s au-dessus de deux immeubles, tourbillon piégé dans la rue.
- **Ballon d'eau chaude stratifié** : la résistance chauffe le bas, l'eau chaude s'accumule en haut.
- **Bassin chauffé, à moitié couvert** : l'air ne reçoit de chaleur qu'au-dessus de la partie découverte.
  (Seul l'air est simulé : la surface de l'eau est une paroi à 28 °C.)

## Fichiers
- `solver.js` : noyau processeur (aussi utilisable sous Node).
- `gpu.js` : noyau WebGPU (shaders WGSL), même interface que `solver.js`.
- `index.src.html` : interface, avec les marqueurs `/*SOLVER*/` et `/*GPU*/`.
- `build.py` : régénère `index.html` (`python3 build.py`).
- `tests/` : essais Node (`npm test`), vérification de stabilité de tous les cas exemples
  (`npm run test:presets`, ou `node tests/presets.js rb` pour un seul cas) et comparaison GPU/processeur
  (`npm run serve`, puis `python3 tests/cmp.py`).

## Dépendances
Aucune pour utiliser le simulateur. Pour les tests : Node.js 24 ou plus récent ; pour la comparaison
GPU/processeur, Playwright (`pip install -r requirements.txt && playwright install chromium`).

## Limites
Fluide « Eau froide » : ρ = ρmax (1 − k (T − 4 °C)²), k = 8·10⁻⁶ K⁻², pour les phénomènes autour de 4 °C.
Viscosité de base relevée au minimum stable du maillage : affiner la grille rapproche de la physique.
Débits et puissances pour la profondeur réglée dans l'onglet Domaine. 2D : turbulence qualitative.
