# Cosserat Shell — Revue mathématique et physique de la formulation

> Auteur : Y. Adagolodjo (DEFROST / INRIA)  
> Date : 2026-05-26  
> Objet : Audit de l'implémentation `src/cosserat_shell/` par rapport au papier  
> *A Geometrically Exact FEM Framework for Cosserat Shells* (`cosserat_shell_v0.pdf`)

---

## Table des matières

1. [Ce qui est mathématiquement correct](#1-ce-qui-est-mathématiquement-correct)
2. [Erreur 1 — Loi constitutive : couplage de Poisson absent](#2-erreur-1--loi-constitutive--couplage-de-poisson-absent)
3. [Erreur 2 — Cisaillement transverse appliqué aux mauvais DOF](#3-erreur-2--cisaillement-transverse-appliqué-aux-mauvais-dof)
4. [Approximation 1 — Strain angulaire sans la carte tangente](#4-approximation-1--strain-angulaire-sans-la-carte-tangente)
5. [Approximation 2 — Raideur géométrique $K_G$ incomplète](#5-approximation-2--raideur-géométrique-kg-incomplète)
6. [Approximation 3 — Type de locking traité](#6-approximation-3--type-de-locking-traité)
7. [Résumé et priorités de correction](#7-résumé-et-priorités-de-correction)

---

## 1. Ce qui est mathématiquement correct

Les éléments suivants sont fidèles au papier et correctement implémentés :

- **Strain twist** $\xi_{t\alpha} = (g_t^{-1}\,\partial g_t/\partial X_\alpha)^\vee$ et mesure de
  déformation $E_{t\alpha} = \xi_{t\alpha} - \xi_{0\alpha}$ (Éqs. 17–19).

- **Forme faible** $G_{\mathrm{int}}$ avec l'opérateur tangent
  $K_\alpha\eta = \partial_\alpha\eta + \mathrm{ad}_{\xi_{t\alpha}}\eta$
  (Lemme 1, Éq. 75), correctement discrétisé par la B-matrice.

- **Raideur matérielle**
  $K_M = \int B_\alpha^T D^{\alpha\beta} B_\beta\,j_0\,dA$ (Éq. 55).

- **Pseudo-inverse** $X_0^* = (X_0^T X_0)^{-1}X_0^T$ et relation
  $\varepsilon = E \cdot X_0^*$ (Éq. 21).

- **Mise à jour Newton** $g_{\mathrm{new}} = g_{\mathrm{old}}\cdot\exp(\hat\eta)$
  avec le Jacobien gauche $J_L$ dans la partie translation (Algorithm 1).

- **Anti-locking** : évaluation des strains au centroïde de l'élément,
  constants par élément (Éq. 63).

---

## 2. Erreur 1 — Loi constitutive : couplage de Poisson absent

### Diagnostic

C'est l'**erreur physique la plus importante** de l'implémentation actuelle.

Pour un matériau isotrope, les stress resultants membranaires $N_{11}$ et $N_{22}$
sont couplés par le coefficient de Poisson $\nu$ :

$$
N_{11} = \frac{Eh}{1-\nu^2}\,\varepsilon_{11} + \frac{\nu Eh}{1-\nu^2}\,\varepsilon_{22}
$$

Dans la formulation 6D, $\varepsilon_{11}$ provient de $E_{t1}$ (direction $X_1$)
et $\varepsilon_{22}$ provient de $E_{t2}$ (direction $X_2$).
Ce couplage passe donc par $D^{12}$, que l'implémentation actuelle met à zéro.

### Ce qui est implémenté (incorrect)

```
D^{12} = D^{21} = 0   →   ν_effectif = 0
```

La coque se comporte comme si $\nu = 0$, tant en membrane qu'en flexion.

### Ce qui devrait être implémenté

En identifiant les composantes physiques dans la convention $[\phi;\rho]$
(angulaire en tête, linéaire en queue) pour une plaque plate avec $R_0 = I$ :

| Terme physique | Source dans $E_{t\beta}$ | Cible dans $S^\alpha$ | Matrice |
|---|---|---|---|
| $\nu\frac{Eh}{1-\nu^2}\,\varepsilon_{22} \to N_{11}$ | $E_{t2}[\rho_y]$ — indice 4 | $S^1[\rho_x]$ — indice 3 | $D^{12}[3,4]$ |
| $\nu\frac{Eh}{1-\nu^2}\,\varepsilon_{11} \to N_{22}$ | $E_{t1}[\rho_x]$ — indice 3 | $S^2[\rho_y]$ — indice 4 | $D^{21}[4,3]$ |
| $\nu D\,\kappa_{22} \to M_{11}$ | $E_{t2}[\phi_y]$ — indice 1 | $S^1[\phi_x]$ — indice 0 | $D^{12}[0,1]$ |
| $\nu D\,\kappa_{11} \to M_{22}$ | $E_{t1}[\phi_x]$ — indice 0 | $S^2[\phi_y]$ — indice 1 | $D^{21}[1,0]$ |

avec $D = Eh^3/[12(1-\nu^2)]$ (raideur de flexion).

### Structure correcte des $D^{\alpha\beta}$

$$
D^{11} = \begin{pmatrix} \mathcal{D}_{11} & 0 \\ 0 & \mathcal{A}_{11} \end{pmatrix},
\qquad
D^{12} = D^{21} = \begin{pmatrix} \mathcal{D}_{12} & 0 \\ 0 & \mathcal{A}_{12} \end{pmatrix}
\neq 0
$$

où les blocs $3\times3$ sont (convention locale : 0=normal à la direction, 1=dans le plan, 2=transverse) :

$$
\mathcal{D}_{12} = \nu D \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix},
\qquad
\mathcal{A}_{12} = \frac{\nu Eh}{1-\nu^2} \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{pmatrix}
$$

> **Note :** Le cisaillement plan symétrisé $N_{12} = N_{21} = Gh\,\varepsilon_{12}$
> nécessite aussi un terme dans $D^{12}$ pour relier $E_{t1}[\rho_y]$
> et $E_{t2}[\rho_x]$ (cf. symétrie de Cauchy $N_{12} = N_{21}$).

### Impact

Sans cette correction, toute simulation avec $\nu \neq 0$ donnera des résultats
**quantitativement faux** — en particulier la déflexion transverse d'une plaque
sera surestimée d'un facteur $\approx 1/(1-\nu^2)$ par rapport à la solution analytique.

---

## 3. Erreur 2 — Cisaillement transverse appliqué aux mauvais DOF

### Diagnostic

Dans `ConstitutiveLaw.h`, la correction de cisaillement transverse `ksGh`
est ajoutée aux indices globaux **[4]** ($\rho_y$) et **[5]** ($\rho_z$)
pour $D^{11}$ et $D^{22}$ de façon identique.

Cependant, la signification physique de ces indices **dépend de la direction $\alpha$**.

### Décomposition physique des $\rho_{t\alpha}$

Pour une plaque plate avec $R_0 = I$ :

$$
\rho_{t1} = R_t^T \frac{\partial\varphi_t}{\partial X_1}
\approx \begin{bmatrix} 1+\varepsilon_{11} \\ \varepsilon_{12} \\ \gamma_{13} \end{bmatrix}
\quad\Rightarrow\quad
E_{t1,\mathrm{lin}} = \begin{bmatrix} \varepsilon_{11} \\ \varepsilon_{12} \\ \gamma_{13} \end{bmatrix}
$$

$$
\rho_{t2} = R_t^T \frac{\partial\varphi_t}{\partial X_2}
\approx \begin{bmatrix} \varepsilon_{21} \\ 1+\varepsilon_{22} \\ \gamma_{23} \end{bmatrix}
\quad\Rightarrow\quad
E_{t2,\mathrm{lin}} = \begin{bmatrix} \varepsilon_{21} \\ \varepsilon_{22} \\ \gamma_{23} \end{bmatrix}
$$

### Effet de l'implémentation actuelle

| Direction | Indice global [4] = $\rho_y$ | Indice global [5] = $\rho_z$ |
|---|---|---|
| $\alpha=1$ : $D^{11}[4,4]$ += `ksGh` | cisaillement **plan** $\varepsilon_{12}$ ❌ | cisaillement **transverse** $\gamma_{13}$ ✅ |
| $\alpha=2$ : $D^{22}[4,4]$ += `ksGh` | extension $\varepsilon_{22}$ ❌❌ | cisaillement **transverse** $\gamma_{23}$ ✅ |

- Pour $\alpha=1$, l'indice [4] correspond au cisaillement plan $\varepsilon_{12}$,
  qui devrait être traité par $G\cdot h$ (pas $\kappa G h$).
- Pour $\alpha=2$, l'indice [4] correspond à l'extension $\varepsilon_{22}$,
  et lui ajouter `ksGh` **pollue directement la raideur d'extension**.

### Correction

La correction transverse doit être **sélective** : appliquer `ksGh` uniquement
sur le DOF correspondant à $\gamma_{1\alpha}$ (hors-plan), soit l'indice **[5]**
($\rho_z$) pour les **deux** directions $\alpha=1$ et $\alpha=2$.
Il faut supprimer l'ajout sur l'indice [4].

```cpp
// Correct : seulement ρ_z = cisaillement transverse (indice 5) pour les deux directions
D11(5, 5) += ksGh;   // D^{11} : γ_13
D22(5, 5) += ksGh;   // D^{22} : γ_23
// NE PAS ajouter sur [4] pour D^{11} ni D^{22}
```

---

## 4. Approximation 1 — Strain angulaire sans la carte tangente

### Formule implémentée

```cpp
phi_talpha ≈ Σ_i  (∂N^i/∂X_α)(0,0) · log(R_base^{-1} · R_i)
```

### Formule exacte

La dérivée de l'interpolation SO(3) nécessite la **carte tangente** (right Jacobian
inverse) de l'application exponentielle :

$$
\phi_{t\alpha} = J_R^{-1}\!\left(\sum_i N^i \log(R_{\mathrm{base}}^{-1} R_i)\right)
\cdot \sum_i \frac{\partial N^i}{\partial X_\alpha}\,\log(R_{\mathrm{base}}^{-1} R_i)
$$

où le Jacobien droit inverse de SO(3) est :

$$
J_R^{-1}(\omega) = I + \frac{1}{2}[\omega]_\times
+ \left(\frac{1}{\theta^2} - \frac{1+\cos\theta}{2\theta\sin\theta}\right)[\omega]_\times^2,
\qquad \theta = \|\omega\|
$$

avec $J_R^{-1}(\omega) \xrightarrow{\theta\to 0} I + \frac{1}{2}[\omega]_\times + \frac{1}{12}[\omega]_\times^2$.

### Quand est-ce que ça pose problème ?

L'approximation $J_R^{-1} \approx I$ est valide pour de **petites rotations inter-nœuds**.
Elle induit une erreur d'ordre $O(\theta^2)$ sur le strain angulaire, où $\theta$ est
l'angle de rotation entre nœuds voisins.

| Rotation inter-nœuds | Erreur sur $\phi_{t\alpha}$ |
|---|---|
| $< 10°$ | $< 1\%$ — négligeable |
| $30°$ | $\approx 4\%$ — acceptable |
| $60°$ | $\approx 15\%$ — significatif |
| $> 90°$ | erreur majeure, risque de non-convergence Newton |

### Impact pratique

Pour les benchmarks du papier (roll-up complet à $2\pi$, grande torsion),
la carte tangente est **indispensable** pour obtenir la convergence quadratique
de Newton et des résultats précis.

---

## 5. Approximation 2 — Raideur géométrique $K_G$ incomplète

### Formule implémentée

$$
K_G^e = \int\sum_\alpha B_\alpha^T\,\bigl(-\mathrm{ad}_{S^\alpha}^T\bigr)\,B_\alpha\,j_0\,dA
$$

### Formule exacte (Éq. 79–80 du papier)

La contribution géométrique au bilinéaire $DG_{\mathrm{int}}\cdot(\eta,\kappa)$ est :

$$
a_G(\eta,\kappa) = \int\sum_\alpha \langle S^\alpha,\,\mathrm{ad}_{K_\alpha\eta}\,\kappa\rangle\,j_0\,dA
$$

En développant avec les constantes de structure de $\mathfrak{se}(3)$,
la matrice élémentaire correcte est :

$$
K_G^e = \int\sum_\alpha B_\alpha^T\,\widetilde{M}_{S^\alpha}\,B_\alpha\,j_0\,dA
$$

où $\widetilde{M}_{S^\alpha}$ est la **matrice de représentation adjointe** de $S^\alpha$,
définie par :

$$
\bigl(\widetilde{M}_{S^\alpha}\bigr)_{ik} = \sum_j f_{jki}\,S^\alpha_j
$$

avec $f_{jki}$ les constantes de structure de $\mathfrak{se}(3)$.
Pour la convention $[\phi;\rho]$ (angulaire en tête) :

$$
\widetilde{M}_{S^\alpha} = \begin{pmatrix} [\phi_{S^\alpha}]_\times & 0 \\ [\rho_{S^\alpha}]_\times & [\phi_{S^\alpha}]_\times \end{pmatrix}
= \mathrm{ad}_{S^\alpha}
$$

### Différence avec l'implémentation

L'implémentation utilise $-\mathrm{ad}_{S^\alpha}^T$ alors que la formule exacte
donne $\mathrm{ad}_{S^\alpha}$ (sans transposition, sans signe moins).

Pour $\mathfrak{se}(3)$ avec métrique plate $\langle a,b\rangle = a^Tb$ :

$$
-\mathrm{ad}_{S^\alpha}^T \neq \mathrm{ad}_{S^\alpha}
\quad \text{(car SE(3) n'est pas semi-simple)}
$$

Les deux coïncident seulement si $\mathrm{ad}$ est skew-adjoint, ce qui est vrai
pour les groupes compacts (SO(3)) mais **pas pour SE(3)**.

### Impact pratique

L'erreur sur $K_G$ est proportionnelle aux stress $S^\alpha$ (pré-contrainte).
Elle affecte :
- La **convergence quadratique** de Newton (dégradée en convergence linéaire)
- La **précision de la raideur** pour de grandes déformations
- La **symétrie de $K$** au-delà des points d'équilibre (Appendice C)

Pour des simulations quasi-statiques avec chargement incrémental et petits pas,
l'impact reste modéré car $K_G \ll K_M$.

---

## 6. Approximation 3 — Type de locking traité

### Ce que le papier revendique

Le papier affirme que l'évaluation au centroïde élimine le **shear-locking**
sans coût supplémentaire.

### Analyse physique

Il existe en réalité plusieurs types de locking pour les éléments de coque :

| Type de locking | Mécanisme | Évaluation centroïde | MITC | EAS / $\bar{B}$ |
|---|---|---|---|---|
| **Membrane locking** | extension parasite en flexion pure | ✅ éliminé | ✅ | ✅ |
| **Shear locking** | cisaillement parasite en flexion mince ($h\to0$) | ⚠️ partiellement | ✅ | ✅ |
| **Volumetric locking** | incompressibilité ($\nu\to0.5$) | ❌ non traité | ❌ | ✅ |
| **Curvature locking** | courbure parasite pour éléments distordus | ⚠️ partiellement | ✅ | — |

L'évaluation au centroïde traite surtout le **membrane locking** (réduction d'intégration
effective des termes membranaires), pas le shear locking au sens strict de Reissner-Mindlin.

Pour le shear locking classique, les méthodes de référence sont :
- **MITC4** (*Mixed Interpolation of Tensorial Components*) : interpolation mixte
  des composantes de cisaillement transverse
- **Reduced integration 1-point** : simple mais génère des modes sablier (*hourglass*)

### Conséquences pratiques

- Pour des coques d'épaisseur modérée ($h/L > 0.05$) : l'approximation centroïde suffit
- Pour des coques minces ($h/L < 0.01$) : possible rigidité parasite résiduelle
- Pour $\nu \to 0.5$ (quasi-incompressible) : locking volumétrique non traité,
  résultats potentiellement bloqués

---

## 7. Résumé et priorités de correction

### Tableau récapitulatif

| Priorité | Problème | Fichier concerné | Impact physique | Difficulté |
|---|---|---|---|---|
| 🔴 **1** | $D^{12} = 0$ — effet Poisson absent | `ConstitutiveLaw.h` | Résultats faux pour $\nu\neq0$ | Moyenne |
| 🔴 **2** | Cisaillement transverse sur mauvais indices | `ConstitutiveLaw.h` | Raideur de flexion/extension polluée | Faible |
| 🟡 **3** | Carte tangente $J_R^{-1}$ absente | `SE3ShellKinematics.h` | Erreur en grandes rotations ($>30°$) | Élevée |
| 🟡 **4** | $K_G$ : $-\mathrm{ad}^T$ vs $\mathrm{ad}$ | `ShellElementComputer.h` | Convergence Newton dégradée | Faible |
| 🟡 **5** | Locking volumétrique non traité | `ConstitutiveLaw.h` + FEM | Problème pour $\nu\to0.5$ | Élevée |

### Ordre de correction recommandé

**Étape 1 (avant tout test numérique) :**
Corriger la loi constitutive — erreurs 1 et 2 dans `ConstitutiveLaw.h`.
Sans cette correction, même un cas simple de plaque cantilever donnera des
résultats quantitativement faux dès que $\nu \neq 0$.

**Étape 2 (validation sur benchmarks du papier) :**
- Corriger la matrice $K_G$ (erreur 4, changement mineur)
- Valider sur cantilever, roll-up et torsion (Sections 6.1–6.3 du papier)

**Étape 3 (grandes déformations) :**
- Ajouter la carte tangente $J_R^{-1}$ (approximation 3)
- Nécessaire pour reproduire les benchmarks à $2\pi$ (roll-up complet)

**Étape 4 (coques minces et quasi-incompressibles) :**
- Implémenter MITC4 ou une formulation EAS pour le locking résiduel
- Utile pour $h/L < 0.01$ ou $\nu > 0.45$

---

*Référence : `cosserat_shell_v0.pdf` (39 pages). Document créé le 2026-05-26.*
