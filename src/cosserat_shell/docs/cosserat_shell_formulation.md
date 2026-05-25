# Cosserat Shell — Formulation mathématique complète

> **Référence :** *A Geometrically Exact Finite Element Framework for Cosserat Shells*  
> Auteurs : M.J. & R.C. — NSERC grant DGECR-2019-00085  
> Notes rédigées par Y. Adagolodjo pour l'implémentation dans `plugin.Cosserat`  
> Dernière révision : 2026-05-25

---

## Convention liegroups (à lire en premier)

Toute la formulation utilise la convention définie dans
`src/liegroups/STRAIN_CONVENTION.md` :

```
ζ = [φ ; ρ]ᵀ  ∈ ℝ⁶
     └─┬─┘  └─┬─┘
  angulaire  linéaire
  head<3>()  tail<3>()
```

| Indice | Symbole | Signification physique |
|--------|---------|------------------------|
| 0 | φ_x | torsion (rotation autour de $X_1$) |
| 1 | φ_y | flexion autour Y |
| 2 | φ_z | flexion autour Z |
| 3 | ρ_x | extension le long de $X_1$ (déviation par rapport à 1) |
| 4 | ρ_y | cisaillement Y |
| 5 | ρ_z | cisaillement Z |

> ⚠️ **Note :** le commentaire Doxygen dans `SE3.h` (`(vx, vy, vz, ωx, ωy, ωz)`) est
> **obsolète**. La convention effective du code est **toujours** [angulaire ; linéaire].
> Voir `STRAIN_CONVENTION.md` §« Correspondance avec se(3) ».

---

## Table des matières

1. [Configuration géométrique](#1-configuration-géométrique)
2. [Strains locaux — twists de déformation](#2-strains-locaux)
3. [Gradient de déformation local $F_e$](#3-gradient-de-déformation-local)
4. [Énergie interne et forme forte](#4-énergie-interne-et-forme-forte)
5. [Stress resultants et loi constitutive](#5-stress-resultants-et-loi-constitutive)
6. [Matrice constitutive $D^{\alpha\beta}$](#6-matrice-constitutive)
7. [Forme faible — principe des travaux virtuels](#7-forme-faible)
8. [Raideur tangente](#8-raideur-tangente)
9. [Chargement magnétique](#9-chargement-magnétique)
10. [Discrétisation FEM](#10-discrétisation-fem)
11. [Algorithme de Newton](#11-algorithme-de-newton)
12. [Système linéaire global](#12-système-linéaire-global)
13. [Appendices — preuves](#13-appendices-preuves)
14. [API liegroups disponible](#14-api-liegroups-disponible)
15. [État naturel $X_0$ — cas de la plaque plate](#15-état-naturel-x0-cas-de-la-plaque-plate)
16. [Domaine de validité](#16-domaine-de-validité)
17. [Paramètres des benchmarks](#17-paramètres-des-benchmarks)

---

## 1. Configuration géométrique

### Domaine de référence

La coque est un **2-manifold** $\mathcal{A} \subset \mathbb{R}^2$, paramétré par les
coordonnées curvilignes $(X_1, X_2)$.

### Configuration actuelle

$$
g_t : \mathcal{A} \to SE(3), \qquad
g_t(X_1, X_2) = \begin{pmatrix} R_t & \varphi_t \\ 0 & 1 \end{pmatrix} \in SE(3) \subset \mathbb{R}^{4\times4}
$$

où :
- $R_t(X_1,X_2) \in SO(3)$ — orientation du repère matériel local
- $\varphi_t(X_1,X_2) \in \mathbb{R}^3$ — position du point matériel

### Configuration de référence (non déformée)

$$
g_0 : \mathcal{A} \to SE(3), \qquad g_0(X_1,X_2) = \begin{pmatrix} R_0 & \varphi_0 \\ 0 & 1 \end{pmatrix}
$$

### Mesure d'aire de référence

$$
j_0(X_1,X_2) = \left\| \frac{\partial \varphi_0}{\partial X_1} \times \frac{\partial \varphi_0}{\partial X_2} \right\| \in \mathbb{R}^+
\qquad \text{(Jacobien scalaire de l'aire)}
$$

---

## 2. Strains locaux

### Strain twist actuel (body-frame)

Pour $\alpha \in \{1, 2\}$ :

$$
\zeta_{t\alpha} = \left(g_t^{-1} \frac{\partial g_t}{\partial X_\alpha}\right)^\vee \in \mathbb{R}^6 \cong \mathfrak{se}(3)
\qquad \text{(Éq. 17)}
$$

Forme explicite en blocs, en appliquant l'opérateur vee (${}^\vee$) à la matrice $4\times4$ :

$$
g_t^{-1} \frac{\partial g_t}{\partial X_\alpha} = \begin{pmatrix} R_t^T \dfrac{\partial R_t}{\partial X_\alpha} & R_t^T \dfrac{\partial \varphi_t}{\partial X_\alpha} \\ 0 & 0 \end{pmatrix}
$$

Le bloc supérieur gauche $R_t^T \partial R_t/\partial X_\alpha$ est une matrice **antisymétrique 3×3**.
L'opérateur vee (${}^\vee$) extrait le vecteur axial :

$$
\zeta_{t\alpha} = \begin{bmatrix}
\underbrace{\left(R_t^T \dfrac{\partial R_t}{\partial X_\alpha}\right)^\vee}_{\phi_{t\alpha} \,\in\, \mathbb{R}^3 \; \text{(angulaire)}} \\[10pt]
\underbrace{R_t^T \dfrac{\partial \varphi_t}{\partial X_\alpha}}_{\rho_{t\alpha} \,\in\, \mathbb{R}^3 \; \text{(linéaire)}}
\end{bmatrix}
\in \mathbb{R}^6
$$

Convention lib liegroups : **indices 0-2 = angulaire** ($\phi$), **indices 3-5 = linéaire** ($\rho$).

### Strain twist de référence

$$
\zeta_{0\alpha} = \left(g_0^{-1} \frac{\partial g_0}{\partial X_\alpha}\right)^\vee \in \mathbb{R}^6
\qquad \text{(Éq. 18)}
$$

### Strain différentiel (mesure de déformation)

$$
E_{t\alpha} = \zeta_{t\alpha} - \zeta_{0\alpha} \in \mathbb{R}^6
\qquad \text{(Éq. 19)}
$$

> **Remarque :** Pour une plaque plate avec $R_0 = I$ et $\varphi_0$ linéaire, on a
> $\zeta_{01} = [0,0,0,1,0,0]^T$ et $\zeta_{02} = [0,0,0,0,1,0]^T$ (voir §15).

---

## 3. Gradient de déformation local

### Matrices de configuration

$$
X_t = \bigl[\zeta_{t1} \,\big|\, \zeta_{t2}\bigr] \in \mathbb{R}^{6\times2},
\qquad
X_0 = \bigl[\zeta_{01} \,\big|\, \zeta_{02}\bigr] \in \mathbb{R}^{6\times2}
$$

### Pseudo-inverse gauche de $X_0$

$X_0 \in \mathbb{R}^{6\times2}$ est supposée de rang 2 (colonnes linéairement indépendantes).
Son **pseudo-inverse gauche** (left pseudo-inverse) est :

$$
X_0^* = \bigl(X_0^T X_0\bigr)^{-1} X_0^T \in \mathbb{R}^{2\times6}
\qquad \text{tel que } X_0^* X_0 = I_2
$$

> ⚠️ **Formule corrigée :** c'est bien $(X_0^T X_0)^{-1} X_0^T$ (produit $2\times2$
> inversé), **pas** $X_0^T(X_0 X_0^T)^{-1}$ qui serait de dimension incompatible.
> L'implémentation dans `ShellElementComputer.h` utilise correctement
> `XtX.ldlt().solve(X0.transpose())`.

### Gradient de déformation local

$$
F_e = X_t \cdot X_0^* \in \mathbb{R}^{6\times6}
\qquad \text{(Éq. 20)}
$$

### Tenseur de déformation (strain tensor)

$$
\varepsilon = F_e - I_e = E \cdot X_0^*
\qquad \text{(Éq. 21)}
$$

où $E = [E_{t1} \,|\, E_{t2}] \in \mathbb{R}^{6\times2}$ et
$I_e = X_0 \cdot X_0^* \in \mathbb{R}^{6\times6}$ est le **projecteur** sur l'espace tangent de référence
(pas l'identité $I_6$ sauf si $X_0$ est de rang plein en lignes, ce qui est faux pour une surface 2D).

> **Note d'implémentation :** $X_0^*$ ne dépend que de $g_0$ — calculé une fois en `init()`.

---

## 4. Énergie interne et forme forte

### Énergie de déformation (fonctionnelle scalaire)

$$
W_{\mathrm{int}} = \frac{1}{2} \int_{\mathcal{A}} \sum_{\alpha,\beta=1}^{2}
\langle E_{t\alpha},\, D^{\alpha\beta} E_{t\beta} \rangle \, j_0 \, dA
\qquad \text{(Éq. 35)}
$$

où $D^{\alpha\beta} \in \mathbb{R}^{6\times6}$ est la matrice constitutive (voir §6).

### Équilibre — forme forte (PDE)

$$
-\sum_{\alpha=1}^{2} \frac{1}{j_0} \frac{\partial(j_0 S^\alpha)}{\partial X_\alpha}
+ \sum_{\alpha=1}^{2} \mathrm{ad}^*_{\zeta_{t\alpha}} S^\alpha + f_{\mathrm{ext}} = 0
\qquad \text{(Éq. 30)}
$$

où $\mathrm{ad}^*_\xi$ est le **co-adjoint** défini par $\langle \mathrm{ad}^*_\xi \eta, \zeta \rangle = \langle \eta, \mathrm{ad}_\xi \zeta \rangle$.

---

## 5. Stress resultants et loi constitutive

### Stress resultants $S^\alpha$ (intégrés sur l'épaisseur)

$$
S^\alpha = \sum_{\beta=1}^{2} D^{\alpha\beta} \cdot E_{t\beta} \in \mathbb{R}^6
\qquad \text{(Éq. 39)}
$$

Interprétation physique de $S^\alpha = [\phi_{S^\alpha};\, \rho_{S^\alpha}]$ :

| Composante | Indices | Signification |
|------------|---------|---------------|
| $\phi_{S^\alpha}$ | 0-2 | Moments (flexion + torsion) par unité de longueur |
| $\rho_{S^\alpha}$ | 3-5 | Forces (extension + cisaillement) par unité de longueur |

---

## 6. Matrice constitutive

### Structure par blocs

En accord avec la convention liegroups, $D^{\alpha\beta} \in \mathbb{R}^{6\times6}$ est **bloc-diagonal** :

$$
D^{\alpha\beta} = \begin{pmatrix} C_{\mathrm{bend}}^{\alpha\beta} & 0 \\ 0 & C_{\mathrm{memb}}^{\alpha\beta} \end{pmatrix}
\qquad \text{(Éqs. 40–41)}
$$

- Bloc **0:3×0:3** ($\phi$-$\phi$) : raideur de **flexion/torsion**
- Bloc **3:6×3:6** ($\rho$-$\rho$) : raideur **membranaire + cisaillement transverse**

### Matériau isotrope linéaire

Pour $E$ (Young), $\nu$ (Poisson), $h$ (épaisseur), $G = E/(2(1+\nu))$,
$\kappa = 5/6$ (facteur de Timoshenko) :

**Raideur membranaire** $C_{\mathrm{memb}} \in \mathbb{R}^{3\times3}$ :

$$
C_{\mathrm{memb}} = \frac{Eh}{1-\nu^2}
\begin{pmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{pmatrix}
$$

**Raideur de flexion** $C_{\mathrm{bend}} \in \mathbb{R}^{3\times3}$ :

$$
C_{\mathrm{bend}} = \frac{Eh^3}{12(1-\nu^2)}
\begin{pmatrix} 1 & \nu & 0 \\ \nu & 1 & 0 \\ 0 & 0 & \frac{1-\nu}{2} \end{pmatrix}
$$

**Cisaillement transverse** (ajouté à la diagonale de $C_{\mathrm{memb}}$, composantes
$\rho_y$, $\rho_z$ = indices 4 et 5 globalement) :

$$
k_s = \kappa \, G \, h \qquad \text{(par unité d'aire)}
$$

**Couplage $\alpha$-$\beta$** pour un matériau **isotrope** sur surface plane :

$$
D^{12} = D^{21} = 0
\qquad \text{(pas de couplage entre directions } X_1 \text{ et } X_2\text{)}
$$

> ⚠️ $D^{12} \neq 0$ pour des matériaux anisotropes ou des surfaces initialement courbées.

---

## 7. Forme faible

### Principe des travaux virtuels

Pour tout champ de variation admissible $\kappa : \mathcal{A} \to \mathfrak{se}(3)$
avec $\hat\kappa = g_t^{-1} \delta g_t$ :

$$
G_{\mathrm{int}}(g_t,\kappa) = \int_{\mathcal{A}} \sum_{\alpha=1}^{2}
\left\langle S^\alpha,\; \frac{\partial \kappa}{\partial X_\alpha} + \mathrm{ad}_{\zeta_{t\alpha}} \kappa \right\rangle j_0 \, dA
\qquad \text{(Éq. 32)}
$$

### Condition d'équilibre

$$
G_{\mathrm{int}}(g_t,\kappa) + G_{\mathrm{ext}}(g_t,\kappa) = 0 \qquad \forall\,\kappa
\qquad \text{(Éq. 33)}
$$

### Travaux virtuels extérieurs

Pour une force distribuée $f \in \mathbb{R}^3$ et un couple distribué $m \in \mathbb{R}^3$ :

$$
G_{\mathrm{ext}} = -\int_{\mathcal{A}} \bigl(\langle m, \kappa_\phi \rangle + \langle f, \kappa_\rho \rangle\bigr) \, j_0 \, dA
$$

où $\kappa = [\kappa_\phi;\, \kappa_\rho]$ avec $\kappa_\phi$ = partie angulaire (indices 0-2),
$\kappa_\rho$ = partie linéaire (indices 3-5).

---

## 8. Raideur tangente

### Dérivation (Proposition 1, Éq. 55)

La dérivée directionnelle de $G_{\mathrm{int}}$ dans la direction $\eta$ s'écrit :

$$
DG_{\mathrm{int}} \cdot \eta =
\int_{\mathcal{A}} \sum_{\alpha=1}^{2} \left[
\left\langle \sum_\beta D^{\alpha\beta} K_\beta\eta,\; K_\alpha \kappa \right\rangle
+ \left\langle S^\alpha,\; \mathrm{ad}_{K_\alpha\eta} \kappa \right\rangle
\right] j_0 \, dA
\qquad \text{(Éq. 79)}
$$

avec l'opérateur tangent $K_\alpha\eta = \dfrac{\partial\eta}{\partial X_\alpha} + \mathrm{ad}_{\zeta_{t\alpha}} \eta$
(Lemme 1, Éq. 75).

### Raideur matérielle $K_M$

$$
K_M = \int_{\mathcal{A}} \sum_{\alpha,\beta} B_\alpha^T D^{\alpha\beta} B_\beta \, j_0 \, dA
$$

où $B_\alpha \in \mathbb{R}^{6\times N_e}$ est la discrétisation de $K_\alpha$ (voir §10).

### Raideur géométrique $K_G$

Le terme $\langle S^\alpha, \mathrm{ad}_{K_\alpha\eta} \kappa \rangle$ se réécrit via l'Éq. 80 :

$$
\langle S^\alpha, \mathrm{ad}_{K_\alpha\eta} \kappa \rangle
= -\langle \mathrm{ad}^*_\kappa S^\alpha, K_\alpha\eta \rangle
= -\kappa^T \bigl(-\mathrm{ad}_{S^\alpha}\bigr)^T \, B_\alpha \eta
$$

La contribution à la raideur géométrique est donc :

$$
K_G = \int_{\mathcal{A}} \sum_{\alpha=1}^{2} B_\alpha^T \bigl(-\mathrm{ad}_{S^\alpha}^T\bigr) B_\alpha \, j_0 \, dA
$$

où $\mathrm{ad}_{S^\alpha} \in \mathbb{R}^{6\times6}$ est la **matrice du crochet de Lie**
$[\,S^\alpha,\, \cdot\,]$, définie en §14.

### Symétrie à l'équilibre (Remarque 2, Appendice C)

La partie anti-symétrique de $K_G$ est nulle aux points d'équilibre statique :
$K = K_M + K_G$ est donc **symétrique** à l'équilibre, ce qui permet l'utilisation de
solveurs symétriques (LDLT, Cholesky).

---

## 9. Chargement magnétique

### Magnétisation rémanente (frame matériel)

$B_0^r \in \mathbb{R}^3$ : vecteur de flux magnétique rémanent dans le **repère matériel**
de la configuration de référence (propriété matérielle, pas une position).

### Champ extérieur appliqué

$B^a \in \mathbb{R}^3$ : champ magnétique externe uniforme (repère global).

### Magnétisation dans la configuration courante

$$
B_t^r = R_t R_0^T B_0^r \in \mathbb{R}^3
\qquad \text{(Éq. 45)}
$$

### Couple magnétique (body couple)

$$
m_t = \frac{1}{\mu_0} B_t^r \times B^a \in \mathbb{R}^3
\qquad \text{(Éq. 46)}
$$

où $\mu_0 = 4\pi \times 10^{-7}$ T·m/A.

### Travail virtuel magnétique

$$
G_{\mathrm{mag}} = -\int_{\mathcal{A}} \langle m_t, \kappa_\phi \rangle \, j_0 \, dA
\qquad \text{(Éq. 47)}
$$

où $\kappa_\phi$ = partie angulaire (indices 0-2) de $\kappa$.

### Raideur magnétique $K_{MG}$

La linéarisation de $G_{\mathrm{mag}}$ contribue à la raideur tangente via :

$$
(K_{MG} - K_M)\,\eta = F_U - F_M
\qquad \text{(Éq. 70)}
$$

- $F_U$ : résidu des forces élastiques internes
- $F_M$ : résidu des forces magnétiques
- $K_{MG}$ : raideur géométrique magnétique (contribution de $\partial m_t/\partial R_t$)

---

## 10. Discrétisation FEM

### Élément de référence

**Quadrilatère isoparamétrique à 4 nœuds** (Q4) sur le domaine $[-1,1]^2$ :

```
Nœuds locaux :  (X_1, X_2) = (-1,-1), (+1,-1), (+1,+1), (-1,+1)
```

### Fonctions de forme bilinéaires

$$
N^i(x,y) = \tfrac{1}{4}(1 + x\,x_i)(1 + y\,y_i), \quad i = 1,\ldots,4
\qquad \text{(Éq. 61)}
$$

### Interpolation de la configuration

Position (linéaire en $\mathbb{R}^3$) :

$$
\varphi_t(x,y) = \sum_{i=1}^{4} N^i(x,y)\, \varphi_t^i
$$

Orientation (interpolation dans $\mathfrak{so}(3)$ — geodésique sur SO(3)) :

$$
R_t(x,y) = R_{\mathrm{base}} \cdot \exp\!\left(\sum_{i=1}^{4} N^i(x,y)\,\log\!\bigl(R_{\mathrm{base}}^{-1} R_t^i\bigr)\right)
$$

où $R_{\mathrm{base}} = R_t^1$ (nœud 0 comme frame de référence locale).

> ⚠️ **Singularité :** $\log(R)$ n'est pas défini de façon unique pour $\|R - I\| \geq \pi$
> (rotations d'angle $\geq\pi$ entre nœuds voisins). Cette interpolation est valide
> en pratique tant que les rotations inter-éléments restent $< \pi$.

### Intégration de Gauss

Quadrature $2\times2$ (intégration complète) :

$$
\text{Points : } \bigl(\pm 1/\sqrt{3},\, \pm 1/\sqrt{3}\bigr), \qquad \text{Poids : } w_g = 1
$$

### Anti-locking — évaluation au centroïde

Pour éviter le **shear-locking**, les strains $\zeta_{t\alpha}$ sont évalués au
**centroïde de l'élément** $(x,y) = (0,0)$ et supposés **constants par élément** :

$$
\zeta_{t\alpha}\big|_e = \zeta_{t\alpha}(0,0)
\qquad \text{(Éq. 63 — centroid strain evaluation)}
$$

Ces strains constants sont utilisés pour **tous les points de Gauss** de l'élément,
ce qui élimine le locking sans coût de calcul supplémentaire.
Analogue à l'hypothèse de déformation constante dans les éléments CST.

### DOF et B-matrices

Pour un maillage de $n_{\mathrm{nœuds}}$ nœuds ($6$ DOF par nœud) :

$$
N_{\mathrm{DOF}} = 6\, n_{\mathrm{nœuds}}
$$

La **B-matrice** discrétisant $K_\alpha\eta$ sur un élément Q4 ($N_e = 24$ DOF) est :

$$
B_\alpha \in \mathbb{R}^{6 \times 24}, \qquad
B_\alpha\big|_{:,\,6i:6i+6} = \frac{\partial N^i}{\partial X_\alpha}\, I_6 + N^i\, \mathrm{ad}_{\zeta_{t\alpha}}
\qquad \text{(discrétisation de } \partial_\alpha + \mathrm{ad}_{\zeta_{t\alpha}})
$$

---

## 11. Algorithme de Newton

```
Initialisation : g_t^0 = g_0  (ou état quasi-statique précédent)

Pour k = 0, 1, 2, … jusqu'à ||R^k|| < ε_tol :

  Boucle sur les éléments e :
    1. Évaluer les strains au centroïde : ζ_{tα}^e  ←  ζ_{tα}(0,0; g_t^k)
    2. Strains différentiels : E^e = X_t^e − X_0^e
    3. Stress resultants : S^{α,e} = Σ_β D^{αβ} E_β^e
    4. B-matrices : B_α^e(x_g, y_g) pour chaque point de Gauss g
    5. Résidu élémentaire : f_e = Σ_α Σ_g B_α^{e,T} S^{α,e} j_0^g w_g
    6. Raideur élémentaire : K_e = K_M^e + K_G^e (Gauss integration)
  Fin boucle éléments

  7. Assemblage global : R^k, K^k
  8. Résoudre : K^k · η = −R^k
  9. Mise à jour Lie (nœud par nœud) :
       R^{k+1} = R^k · exp_SO3(δω),    δω = η[0:3]
       φ^{k+1} = φ^k + δv,             δv = η[3:6]
```

### Mise à jour SE(3) explicite (Éq. Algorithm 1)

$$
g_{\mathrm{new}} = g_{\mathrm{old}} \cdot \exp_{SE(3)}(\hat\eta),
\qquad
\exp_{SE(3)}\!\begin{pmatrix}\delta\omega \\ \delta v\end{pmatrix}
= \begin{pmatrix} \exp_{SO(3)}(\delta\omega) & J_L(\delta\omega)\,\delta v \\ 0 & 1 \end{pmatrix}
$$

où $J_L(\omega)$ est le **Jacobien gauche de SO(3)** :

$$
J_L(\omega) = I + \frac{1-\cos\theta}{\theta^2}[\omega]_\times + \frac{\theta - \sin\theta}{\theta^3}[\omega]_\times^2,
\quad \theta = \|\omega\|
$$

> ⚠️ C'est $J_L$ (Jacobien **gauche**), pas $J_r$ (Jacobien droit).
> Dans la lib liegroups, `SE3::computeExp` utilise `expCosseratGeneral` qui
> implémente exactement cette formule.

---

## 12. Système linéaire global

### Assemblage standard FEM

$$
K_{\mathrm{global}}\, \eta_{\mathrm{global}} = F_{\mathrm{global}}
$$

Les matrices élémentaires $K_e$ (24×24) et vecteurs $f_e$ (24×1) sont assemblés
via la table de connectivité.

### Conditions aux limites de Dirichlet

Imposer $\eta = 0$ sur les DOF bloqués (nœuds encastrés, 6 DOF par nœud, ou
sous-ensemble selon les conditions physiques).

### Avec actuation magnétique (Éq. 70)

$$
(K_{MG} - K_M)\,\eta = F_U - F_M
$$

---

## 13. Appendices — preuves

### Appendice A — Variation de $\zeta_{t\alpha}$ (Lemme 1)

**Résultat** (Éq. 75) :

$$
\delta\hat\zeta_{t\alpha} = \frac{\partial\hat\kappa}{\partial X_\alpha} + [\hat\zeta_{t\alpha},\hat\kappa]
= \left(\frac{\partial\kappa}{\partial X_\alpha} + \mathrm{ad}_{\zeta_{t\alpha}}\kappa\right)^\wedge
$$

**Dérivation** : depuis la définition $\hat\kappa = g_t^{-1}\delta g_t$ (Éq. 28) et la règle de chaîne (Éq. 73) :

$$
\delta\hat\zeta_{t\alpha} = -g_t^{-1}\delta g_t\, g_t^{-1}\frac{\partial g_t}{\partial X_\alpha}
+ g_t^{-1}\frac{\partial(\delta g_t)}{\partial X_\alpha}
$$

Puis en substituant $g_t^{-1}\partial(\delta g_t)/\partial X_\alpha = \partial\hat\kappa/\partial X_\alpha + g_t^{-1}(\partial g_t/\partial X_\alpha)\hat\kappa$ (Éq. 74) :

$$
\delta\hat\zeta_{t\alpha} = \frac{\partial\hat\kappa}{\partial X_\alpha} + [\hat\zeta_{t\alpha},\hat\kappa] \qquad \square
$$

### Appendice B — Raideur tangente (Proposition 1)

La dérivée $\partial/\partial\varepsilon|_{\varepsilon=0}\,G_{\mathrm{int}}(\mathbf{g}_\varepsilon, \kappa)$
avec $\mathbf{g}_\varepsilon = g_t \cdot \exp(\varepsilon\hat\eta)$ donne (Éq. 79) :

$$
= \int_{\mathcal{A}} \sum_\alpha
\Bigl[\langle\textstyle\sum_\beta D^{\alpha\beta} K_\beta\eta, K_\alpha\kappa\rangle
+ \langle S^\alpha, \mathrm{ad}_{K_\alpha\eta}\kappa\rangle\Bigr] j_0\,dA
$$

Puis via (Éq. 80) $\langle S^\alpha, \mathrm{ad}_{K_\alpha\eta}\kappa\rangle
= -\langle \widetilde{\mathrm{ad}}^*_{S^\alpha}\kappa, K_\alpha\eta\rangle$ : $\square$

### Appendice C — Symétrie de $K_G$ à l'équilibre (Remarque 2)

La partie antisymétrique de $DG_{\mathrm{int}}$ vaut après intégration par parties (Éq. 83) :

$$
\mathrm{Skew}\!\left[\frac{\partial G_{\mathrm{int}}}{\partial\varepsilon}\right]
= \int_{\mathcal{A}} \left\langle\left(-\frac{1}{j_0}\frac{\partial(j_0 S^\alpha)}{\partial X_\alpha}
+ \mathrm{ad}^*_{\zeta_{t\alpha}} S^\alpha\right),\, \mathrm{ad}_\kappa\eta\right\rangle j_0\,dA
$$

À l'équilibre, le terme entre parenthèses = $f_{\mathrm{ext}}$ (force résiduelle),
qui est nul à l'équilibre statique sans forces extérieures. Donc $K_G$ est symétrique. $\square$

---

## 14. API liegroups disponible

Voir `src/liegroups/SE3.h`, `SO3.h` et `STRAIN_CONVENTION.md`.

```cpp
using SO3d = sofa::component::cosserat::liegroups::SO3<double>;
using SE3d = sofa::component::cosserat::liegroups::SE3<double>;

// ── SO(3) ──────────────────────────────────────────────────────────────────
SO3d R    = SO3d::identity();
SO3d R2   = SO3d::exp(phi);         // phi ∈ ℝ³  (axe×angle)
Eigen::Vector3d phi = R.log();      // → ℝ³
SO3d Rinv = R.inverse();            // R^T
SO3d Rc   = R * R2;                 // composition
Eigen::Vector3d v2 = R.act(v);      // R·v
Eigen::Matrix3d M  = R.toRotationMatrix();

// ── SE(3) ──────────────────────────────────────────────────────────────────
// Convention strain = [φ; ρ] = [angulaire (head<3>); linéaire (tail<3>)]
SE3d g    = SE3d::identity();
SE3d g2   = SE3d::exp(xi);          // xi ∈ ℝ⁶ = [φ; ρ]
SE3d gCos = SE3d::expCosserat(xi, s); // avec élongation nominale (pour poutre)
Eigen::Matrix<double,6,1> xi = g.log();
SE3d ginv = g.computeInverse();
SE3d gc   = g.compose(g2);          // g·g2
Eigen::Matrix<double,6,6> Ad = g.adjoint(); // Ad_g ∈ ℝ^{6×6}
```

### Matrice $\mathrm{ad}_\xi$ — opérateur Lie bracket

Pour $\xi = [\phi;\,\rho] \in \mathbb{R}^6$ (convention lib : angulaire en tête) :

$$
\mathrm{ad}_\xi = \begin{pmatrix} [\phi]_\times & 0 \\ [\rho]_\times & [\phi]_\times \end{pmatrix}
\in \mathbb{R}^{6\times6}
$$

$$
\mathrm{ad}_\xi\cdot\eta = \begin{bmatrix} \phi\times\phi_\eta \\ \rho\times\phi_\eta + \phi\times\rho_\eta \end{bmatrix}
$$

### Co-adjoint $\mathrm{ad}^*_\xi$

$$
\mathrm{ad}^*_\xi = -\mathrm{ad}_\xi^T =
\begin{pmatrix} -[\phi]_\times & -[\rho]_\times \\ 0 & -[\phi]_\times \end{pmatrix}
$$

---

## 15. État naturel $X_0$ — cas de la plaque plate

Pour une plaque plate de référence avec :
- $R_0 = I_3$ (identité) en tout point
- $\varphi_0(X_1,X_2) = (X_1, X_2, 0)^T$

les strains de référence au centroïde d'un élément sont :

$$
\zeta_{01} = \left(g_0^{-1}\frac{\partial g_0}{\partial X_1}\right)^\vee
= \begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \\ 0 \\ 0 \end{bmatrix},
\qquad
\zeta_{02} = \left(g_0^{-1}\frac{\partial g_0}{\partial X_2}\right)^\vee
= \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 1 \\ 0 \end{bmatrix}
$$

Donc :

$$
X_0 = \begin{pmatrix} 0 & 0 \\ 0 & 0 \\ 0 & 0 \\ 1 & 0 \\ 0 & 1 \\ 0 & 0 \end{pmatrix},
\qquad
X_0^* = (X_0^T X_0)^{-1}X_0^T = \begin{pmatrix} 0&0&0&1&0&0 \\ 0&0&0&0&1&0 \end{pmatrix}
$$

**Vérification :** $X_0^* X_0 = I_2$ ✓

**Interprétation :** Dans l'état naturel, les DOF angulaires (indices 0-2) sont nuls
(pas de courbure), les extensions sont unitaires dans $X_1$ et $X_2$ (indices 3 et 4).
Le strain différentiel $E_{t\alpha} = \zeta_{t\alpha} - \zeta_{0\alpha}$ mesure les
**déviations** par rapport à cet état.

---

## 16. Domaine de validité

Le modèle combine **grandes rotations** (via $SE(3)$) avec une **loi constitutive linéaire**
($D^{\alpha\beta}$ constant).

| Hypothèse | Domaine valide | Hors domaine |
|-----------|---------------|--------------|
| Loi constitutive linéaire | petites déformations matérielles ($\|E\| \ll 1$) | hyper-élasticité, plasticité |
| Grandes rotations | rotations arbitraires via $SO(3)$ | aucune limite |
| Isotropie | matériaux homogènes isotropes | fibres, composites (→ $D^{12} \neq 0$) |
| Interpolation $SO(3)$ | rotations inter-nœuds $< \pi$ | maillages trop grossiers, grandes déformations locales |
| Anti-locking centroïde | courbure constante par élément | gradients de courbure élevés (→ raffiner) |

Typiquement valide pour : **silicone, caoutchouc, coques élastiques** sous grands déplacements
avec petites déformations (type Reissner-Mindlin géométriquement exact).

---

## 17. Paramètres des benchmarks

| Test | Géométrie | $E$ (kPa) | $\nu$ | $h$ (mm) | Éléments | Chargement |
|------|-----------|-----------|-------|----------|----------|------------|
| Cantilever | $100\times10$ mm | 1000 | 0.3 | 0.1 | $10\times1$ | Force au bout |
| Bending roll-up | rectangle | 1000 | 0.0 | — | $10\times1$ | Moment distribué ($2\pi$) |
| Torsion | rectangle | 1000 | 0.3 | — | $10\times4$ | Angle au bout |
| Arc circulaire | arc 90° | — | — | — | 20 | Force latérale |
| Plaque magnétique $L/h=41$ | $30\times5$ mm | — | 0.3 | — | 30 | $B^a=0.05$ T, $B^r=0.143$ T, $\mu=303$ kPa, $\lambda=7300$ kPa |
| Soft gripper magnétique | plaques flexibles | — | — | — | — | champ antiparallèle |

> Les cases `—` indiquent des paramètres non explicitement listés dans le papier
> (dépendent des expériences Zhao et al. 2019 ou Dadgar-Rad et al. 2022).

---

*Référence : `cosserat_shell_v0.pdf` (39 pages). Révisé le 2026-05-25.*
