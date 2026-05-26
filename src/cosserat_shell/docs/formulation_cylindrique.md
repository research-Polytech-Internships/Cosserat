# Cosserat Shell — Formulation en coordonnées cylindriques

> Auteur : Y. Adagolodjo (DEFROST / INRIA)  
> Date : 2026-05-26  
> Prérequis : `cosserat_shell_formulation.md` + `formulation_review.md`  
> Motivation : coques tubulaires (cathéters, robots souples, tubes endoscopiques)

---

## Table des matières

1. [Motivation et géométrie cible](#1-motivation-et-géométrie-cible)
2. [Paramétrisation et repère de référence](#2-paramétrisation-et-repère-de-référence)
3. [Configuration de référence $g_0(\theta, z)$](#3-configuration-de-référence)
4. [Strains twists de référence $\xi_{0\alpha}$](#4-strains-twists-de-référence)
5. [Métrique et jacobien d'aire](#5-métrique-et-jacobien-daire)
6. [Strains différentiels et énergie](#6-strains-différentiels-et-énergie)
7. [Loi constitutive en coordonnées cylindriques](#7-loi-constitutive)
8. [Opérateur tangent $K_\alpha$ et connexion géométrique](#8-opérateur-tangent-et-connexion-géométrique)
9. [B-matrices en coordonnées cylindriques](#9-b-matrices)
10. [Corrections issues de `formulation_review.md`](#10-corrections-issues-de-formulation_reviewmd)
11. [Cas particuliers](#11-cas-particuliers)
12. [Résumé des modifications par rapport au cas cartésien](#12-résumé-des-modifications)
13. [Notes FEM et implémentation](#13-notes-fem-et-implémentation)

---

## 1. Motivation et géométrie cible

La formulation cartésienne $(X_1, X_2)$ est adaptée aux plaques planes ou aux coques
faiblement courbées. Pour les **coques tubulaires cylindriques** — cathéters,
robots continus tubulaires, manchons pneumatiques — le choix naturel est
les **coordonnées cylindriques** $(\theta, z)$.

Les avantages sont :
- Maillage naturel sur la surface cylindrique (quadrilatères $\Delta\theta \times \Delta z$)
- Conditions aux limites axisymétriques simplifiées
- Couplage membrane-flexion dû à la courbure initiale capturé **automatiquement**
  par $\xi_{0\theta} \neq 0$ dans le cadre SE(3)

> **Avantage clé du cadre SE(3) :** les symboles de Christoffel de la géométrie
> cylindrique (termes de connexion de la dérivée covariante) sont automatiquement
> encodés dans $\mathrm{ad}_{\xi_{0\alpha}}$ — aucun terme de courbure additionnel
> n'est à ajouter manuellement à la forme faible.

---

## 2. Paramétrisation et repère de référence

### Domaine de référence

$$
\mathcal{A} = [0,\,2\pi] \times [0,\,L]
\qquad \text{avec} \quad (X_1, X_2) = (\theta, z)
$$

### Repère matériel local

À chaque point $(\theta, z)$ du cylindre de rayon $R$, on définit
les vecteurs de base physiques :

$$
\hat{e}_\theta = \begin{pmatrix}-\sin\theta\\\cos\theta\\0\end{pmatrix},
\qquad
\hat{e}_z = \begin{pmatrix}0\\0\\1\end{pmatrix},
\qquad
\hat{e}_r = \begin{pmatrix}\cos\theta\\\sin\theta\\0\end{pmatrix}
$$

- $\hat{e}_\theta$ : direction **circonférentielle** (tangente au cercle)
- $\hat{e}_z$ : direction **axiale**
- $\hat{e}_r$ : **normale** sortante à la surface (direction radiale)

Le repère matériel de référence est donc :

$$
R_0(\theta) = \bigl[\hat{e}_\theta \;\big|\; \hat{e}_z \;\big|\; \hat{e}_r\bigr]
= \begin{pmatrix}
-\sin\theta & 0 & \cos\theta \\
 \cos\theta & 0 & \sin\theta \\
 0          & 1 & 0
\end{pmatrix}
$$

avec la convention lib liegroups $[\phi;\rho]$ : la **première colonne** est alignée
avec $X_1 = \theta$, la **deuxième** avec $X_2 = z$.

---

## 3. Configuration de référence

La configuration de référence $g_0 : \mathcal{A} \to SE(3)$ est :

$$
g_0(\theta, z) = \begin{pmatrix} R_0(\theta) & \varphi_0(\theta,z) \\ 0 & 1 \end{pmatrix},
\qquad
\varphi_0(\theta,z) = \begin{pmatrix} R\cos\theta \\ R\sin\theta \\ z \end{pmatrix}
$$

**Dérivées partielles** (nécessaires pour $\xi_{0\alpha}$) :

$$
\frac{\partial\varphi_0}{\partial\theta} = R\begin{pmatrix}-\sin\theta\\\cos\theta\\0\end{pmatrix} = R\,\hat{e}_\theta,
\qquad
\frac{\partial\varphi_0}{\partial z} = \begin{pmatrix}0\\0\\1\end{pmatrix} = \hat{e}_z
$$

$$
\frac{\partial R_0}{\partial\theta} = \begin{pmatrix}
-\cos\theta & 0 & -\sin\theta \\
-\sin\theta & 0 &  \cos\theta \\
 0          & 0 & 0
\end{pmatrix},
\qquad
\frac{\partial R_0}{\partial z} = 0
$$

---

## 4. Strains twists de référence

### Calcul explicite de $\xi_{0\theta}$

**Partie angulaire** (vee de $R_0^T\,\partial R_0/\partial\theta$) :

$$
R_0^T\,\frac{\partial R_0}{\partial\theta}
= \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ -1 & 0 & 0 \end{pmatrix}
= [\hat{y}_{\mathrm{mat}}]_\times
\qquad\Rightarrow\qquad \phi_{0\theta} = (0,\,1,\,0)^T
$$

Interprétation : courbure autour de la **deuxième direction matérielle** ($\hat{e}_z$),
soit la courbure axiale du cylindre. Valeur dimensionnelle : $1/R$ par unité d'arc.

**Partie linéaire** ($R_0^T\,\partial\varphi_0/\partial\theta$) :

$$
R_0^T\,\frac{\partial\varphi_0}{\partial\theta} = R_0^T\,(R\hat{e}_\theta) = R\,(1,0,0)^T
\qquad\Rightarrow\qquad \rho_{0\theta} = (R,\,0,\,0)^T
$$

Interprétation : extension de $R$ mètres par radian dans la première direction matérielle.

### Calcul explicite de $\xi_{0z}$

**Partie angulaire** : $\partial R_0/\partial z = 0$ $\Rightarrow$ $\phi_{0z} = (0,0,0)^T$

**Partie linéaire** :

$$
R_0^T\,\frac{\partial\varphi_0}{\partial z} = R_0^T\,\hat{e}_z = (0,\,1,\,0)^T
\qquad\Rightarrow\qquad \rho_{0z} = (0,\,1,\,0)^T
$$

Interprétation : extension unité dans la **deuxième direction matérielle** ($\hat{e}_z$).

### Résumé — strains de référence

$$
\boxed{
\xi_{0\theta} = \begin{pmatrix}0\\1\\0\\R\\0\\0\end{pmatrix},
\qquad
\xi_{0z} = \begin{pmatrix}0\\0\\0\\0\\1\\0\end{pmatrix}
}
$$

Convention : $\xi = [\phi;\rho]$ (angulaire en tête, indices 0-2 ;
linéaire en queue, indices 3-5).

**Comparaison avec la plaque plate** (§15 de `cosserat_shell_formulation.md`) :

| | Plaque plate | Cylindre |
|---|---|---|
| $\xi_{01}$ | $[0,0,0,1,0,0]^T$ | $[0,1,0,R,0,0]^T$ |
| $\xi_{02}$ | $[0,0,0,0,1,0]^T$ | $[0,0,0,0,1,0]^T$ |

La différence clé : $\phi_{0\theta} = (0,1,0)^T \neq 0$ encode la **courbure initiale du cylindre**.
La direction $z$ est identique à la plaque plate (pas de courbure axiale).

---

## 5. Métrique et jacobien d'aire

### Tenseur métrique

En coordonnées cylindriques $(\theta, z)$, le tenseur métrique de la surface est :

$$
g_{\alpha\beta} = \begin{pmatrix} \|\partial\varphi_0/\partial\theta\|^2 & 0 \\ 0 & 1 \end{pmatrix}
= \begin{pmatrix} R^2 & 0 \\ 0 & 1 \end{pmatrix}
$$

### Jacobien d'aire

$$
j_0 = \sqrt{\det g_{\alpha\beta}} = R
$$

L'élément d'aire est $dA = j_0\,d\theta\,dz = R\,d\theta\,dz$, ce qui donne bien
l'aire d'un cylindre : $\int_0^{2\pi}\int_0^L R\,d\theta\,dz = 2\pi R L$.

### Facteurs d'échelle des strains

Les strains twists $\xi_{t\theta}$ sont **par radian** (pas par mètre).
Les strains physiques (par unité d'arc-length) sont :

$$
\tilde{\xi}_{t\theta} = \frac{1}{R}\,\xi_{t\theta}
\quad\text{[par mètre]},
\qquad
\tilde{\xi}_{tz} = \xi_{tz}
\quad\text{[par mètre]}
$$

Ces facteurs d'échelle doivent être absorbés dans les matrices constitutives.

---

## 6. Strains différentiels et énergie

### Strains différentiels

$$
E_{t\theta} = \xi_{t\theta} - \xi_{0\theta} \in \mathbb{R}^6,
\qquad
E_{tz} = \xi_{tz} - \xi_{0z} \in \mathbb{R}^6
$$

### Énergie de déformation en coordonnées cylindriques

$$
W_{\mathrm{int}} = \frac{1}{2}\int_0^{2\pi}\int_0^L
\sum_{\alpha,\beta \in \{\theta,z\}}
\langle E_{t\alpha},\,D_{\mathrm{cyl}}^{\alpha\beta}\,E_{t\beta}\rangle\,R\,d\theta\,dz
$$

La présence de $j_0 = R$ dans l'intégrale **et** dans la définition des
$D_{\mathrm{cyl}}^{\alpha\beta}$ doit être traitée conjointement pour que
l'énergie soit dimensionnellement correcte.

---

## 7. Loi constitutive

### Principe de la mise à l'échelle

Pour que l'énergie coïncide avec la formulation cartésienne en arc-length $s_\theta = R\theta$ :

$$
\int_0^{2\pi}\int_0^L \langle E_{t\theta}, D_{\mathrm{cyl}}^{\theta\theta} E_{t\theta}\rangle\,R\,d\theta\,dz
= \int_0^{2\pi R}\int_0^L \langle \tilde{E}_{ts}, \tilde{D} \tilde{E}_{ts}\rangle\,ds\,dz
$$

Avec $\tilde{E}_{ts} = E_{t\theta}/R$ et $ds = R\,d\theta$ :

$$
\Rightarrow \quad D_{\mathrm{cyl}}^{\theta\theta} = \frac{1}{R^2}\,D_{\mathrm{Cart}}
$$

De même :

$$
D_{\mathrm{cyl}}^{zz} = D_{\mathrm{Cart}},
\qquad
D_{\mathrm{cyl}}^{\theta z} = D_{\mathrm{cyl}}^{z\theta} = \frac{1}{R}\,D_{\mathrm{Cart}}
$$

### Matrices $D^{\alpha\beta}_{\mathrm{cyl}}$ explicites

Pour un matériau isotrope d'épaisseur $h$, Young $E$, Poisson $\nu$,
facteur de Timoshenko $\kappa = 5/6$ :

Définissons les raideurs de base (en arc-length) :

$$
A_m = \frac{Eh}{1-\nu^2},\quad
B_m = \frac{\nu Eh}{1-\nu^2},\quad
C_m = \frac{Gh}{2},\quad
D_b = \frac{Eh^3}{12(1-\nu^2)},\quad
C_b = \frac{Gh^3}{24},\quad
k_s = \kappa Gh
$$

**Matrice $D^{\theta\theta}_{\mathrm{cyl}}$ — direction circonférentielle :**

$$
D^{\theta\theta}_{\mathrm{cyl}} = \frac{1}{R^2}
\begin{pmatrix}
D_b & 0 & 0 & 0 & 0 & 0 \\
0 & C_b & 0 & 0 & 0 & 0 \\
0 & 0 & C_b & 0 & 0 & 0 \\
0 & 0 & 0 & A_m & 0 & 0 \\
0 & 0 & 0 & 0 & C_m & 0 \\
0 & 0 & 0 & 0 & 0 & k_s
\end{pmatrix}
$$

**Matrice $D^{zz}_{\mathrm{cyl}}$ — direction axiale :**

$$
D^{zz}_{\mathrm{cyl}} =
\begin{pmatrix}
C_b & 0 & 0 & 0 & 0 & 0 \\
0 & D_b & 0 & 0 & 0 & 0 \\
0 & 0 & C_b & 0 & 0 & 0 \\
0 & 0 & 0 & C_m & 0 & 0 \\
0 & 0 & 0 & 0 & A_m & 0 \\
0 & 0 & 0 & 0 & 0 & k_s
\end{pmatrix}
$$

**Matrice $D^{\theta z}_{\mathrm{cyl}}$ — couplage de Poisson et cisaillement :**

$$
D^{\theta z}_{\mathrm{cyl}} = D^{z\theta}_{\mathrm{cyl}} = \frac{1}{R}
\begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & C_b & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & B_m & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

> **Note :** $D^{\theta z}_{\mathrm{cyl}}[4,3] = B_m/R$ encode le **couplage de Poisson**
> entre extension axiale $\varepsilon_{zz}$ (depuis $E_{tz}[\rho_y]$, indice 4)
> et force circonférentielle $N_{\theta\theta}$ (dans $S^\theta[\rho_x]$, indice 3).
> Ce terme était absent de la formulation cartésienne incorrecte (erreur 1 de
> `formulation_review.md`) et est ici explicitement inclus.

### Interprétation des indices (convention $[\phi;\rho]$)

La correspondance entre indices 6D et grandeurs physiques dépend de la direction $\alpha$ :

| Indice | Symbole | Direction $\theta$ | Direction $z$ |
|--------|---------|-------------------|---------------|
| 0 | $\phi_x$ | torsion circonférentielle $\tau_\theta$ | torsion axiale $\tau_z$ |
| 1 | $\phi_y$ | **courbure $\kappa_{\theta z}$** (autour de $\hat{e}_z$) | courbure $\kappa_{zz}$ |
| 2 | $\phi_z$ | courbure $\kappa_{\theta r}$ (hors-plan) | courbure $\kappa_{zr}$ |
| 3 | $\rho_x$ | extension $\varepsilon_{\theta\theta}$ / $R$ | cisaillement plan $\varepsilon_{z\theta}$ |
| 4 | $\rho_y$ | cisaillement plan $\varepsilon_{\theta z}$ | extension $\varepsilon_{zz}$ |
| 5 | $\rho_z$ | **cisaillement transverse** $\gamma_{\theta r}$ | **cisaillement transverse** $\gamma_{zr}$ |

> ⚠️ Le **cisaillement transverse physique** (hors-plan) est **toujours** l'indice **5** ($\rho_z$)
> pour les deux directions — cohérent avec la correction de `formulation_review.md` §3.

---

## 8. Opérateur tangent et connexion géométrique

### L'opérateur $K_\alpha\eta$ en cylindrique

$$
K_\theta\eta = \frac{\partial\eta}{\partial\theta} + \mathrm{ad}_{\xi_{t\theta}}\,\eta,
\qquad
K_z\eta = \frac{\partial\eta}{\partial z} + \mathrm{ad}_{\xi_{tz}}\,\eta
$$

### La matrice $\mathrm{ad}_{\xi_{0\theta}}$ — connexion cylindrique

Pour $\xi_{0\theta} = [0,1,0;\,R,0,0]^T$ :

$$
\mathrm{ad}_{\xi_{0\theta}} = \begin{pmatrix} [\phi_{0\theta}]_\times & 0 \\ [\rho_{0\theta}]_\times & [\phi_{0\theta}]_\times \end{pmatrix}
$$

avec :

$$
[\phi_{0\theta}]_\times = \begin{pmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ -1 & 0 & 0 \end{pmatrix},
\qquad
[\rho_{0\theta}]_\times = \begin{pmatrix} 0 & 0 & 0 \\ 0 & 0 & -R \\ 0 & R & 0 \end{pmatrix}
$$

soit :

$$
\mathrm{ad}_{\xi_{0\theta}} =
\begin{pmatrix}
 0 & 0 & 1 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 0 \\
-1 & 0 & 0 & 0 & 0 & 0 \\
 0 & 0 & 0 & 0 & 0 & 1 \\
 0 & 0 & -R & 0 & 0 & 0 \\
 0 & R & 0 & 0 & 0 & 0
\end{pmatrix}
$$

**Interprétation physique de $\mathrm{ad}_{\xi_{0\theta}}$ :**

Ce bloc encode la **dérivée covariante** en $\theta$ sur le cylindre.
En particulier, le terme $-R$ en position $(5, 3)$ (ligne $\rho_z$, colonne $\rho_x$)
représente le couplage entre extension circonférentielle et cisaillement transverse —
c'est l'effet de la courbure cylindrique sur le calcul des forces.

### La matrice $\mathrm{ad}_{\xi_{0z}}$ — direction axiale

Pour $\xi_{0z} = [0,0,0;\,0,1,0]^T$ :

$$
\mathrm{ad}_{\xi_{0z}} =
\begin{pmatrix}
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0
\end{pmatrix}
$$

Seul le terme $(3,5)$ est non nul : $[\rho_{0z}]_\times = [0,1,0]_\times$.
La direction axiale est essentiellement plate (comme la plaque cartésienne).

---

## 9. B-matrices en coordonnées cylindriques

La B-matrice discrétisant $K_\alpha\eta$ sur un élément Q4 reste formellement identique
à la formulation cartésienne (§10 de `cosserat_shell_formulation.md`) :

$$
B_\alpha\big|_{:,\,6i:6i+6} = \frac{\partial N^i}{\partial X_\alpha}\,I_6
+ N^i\,\mathrm{ad}_{\xi_{t\alpha}}
\qquad \in \mathbb{R}^{6\times24}
$$

**Ce qui change en cylindrique :**

Le terme $\mathrm{ad}_{\xi_{t\alpha}}$ contient maintenant la **courbure initiale $1/R$**
via $\xi_{0\theta}$. Au repos ($\xi_{t\theta} = \xi_{0\theta}$) :

$$
B_\theta\big|_{\mathrm{repos}} =
\sum_i \left[\frac{\partial N^i}{\partial\theta}\,I_6 + N^i\,\mathrm{ad}_{\xi_{0\theta}}\right]
$$

Les termes en $\mathrm{ad}_{\xi_{0\theta}}$ jouent le rôle des **termes de courbure**
que la formulation classique de coque ajoute manuellement à l'opérateur différentiel.

**Anti-locking en cylindrique :** L'évaluation au centroïde reste valide.
Pour un élément $[\theta_0 - \Delta\theta/2,\,\theta_0 + \Delta\theta/2]\times[z_0-\Delta z/2,\,z_0+\Delta z/2]$,
le centroïde est $(\theta_0, z_0)$ et $\xi_{t\alpha}$ y est constant par élément.

---

## 10. Corrections issues de `formulation_review.md`

### Erreur 1 — Couplage de Poisson (maintenant corrigé)

Le couplage de Poisson entre $\varepsilon_{\theta\theta}$ et $\varepsilon_{zz}$
est capturé par $D^{\theta z}_{\mathrm{cyl}}[4,3] = B_m/R$ (voir §7).

Sans ce terme, la coque cylindrique sous pression interne ne présente pas
d'allongement axial, ce qui est physiquement faux.

### Erreur 2 — Cisaillement transverse (confirmé)

L'indice physique du cisaillement transverse est **[5]** ($\rho_z$)
pour **les deux** directions $\theta$ et $z$ — confirmé par le tableau §7.

La correction `ksGh` s'applique **uniquement** aux éléments diagonaux $[5,5]$
de $D^{\theta\theta}_{\mathrm{cyl}}$ et $D^{zz}_{\mathrm{cyl}}$ (inclus dans les matrices ci-dessus).

### Approximation 3 — Carte tangente $J_R^{-1}$ (critique en cylindrique)

En cylindrique, la rotation inter-nœuds en $\theta$ est typiquement
$\Delta\phi \approx \Delta\theta$ (pour un cylindre). Pour un maillage de $N_\theta$
éléments : $\Delta\phi \approx 2\pi/N_\theta$.

| $N_\theta$ (éléments) | $\Delta\phi$ | Erreur sur $\phi_{t\theta}$ |
|---|---|---|
| 4 | $90°$ | $\sim 25\%$ — **inacceptable** |
| 8 | $45°$ | $\sim 8\%$ — significatif |
| 16 | $22.5°$ | $\sim 2\%$ — acceptable |
| 32 | $11.25°$ | $\sim 0.5\%$ — bon |

> ⚠️ **La carte tangente $J_R^{-1}$ est obligatoire en coordonnées cylindriques**
> pour des maillages raisonnables ($N_\theta < 32$). En cartésien elle est optionnelle ;
> en cylindrique elle est critique.

---

## 11. Cas particuliers

### 11.1 Chargement axisymétrique ($\partial/\partial\theta = 0$)

Toutes les quantités sont indépendantes de $\theta$. L'opérateur tangent devient :

$$
K_\theta\eta = \mathrm{ad}_{\xi_{t\theta}}\,\eta
\qquad (\partial_\theta = 0)
$$

Le problème se réduit à une dimension ($z$ uniquement), avec un terme de "couplage hoop"
via $\mathrm{ad}_{\xi_{0\theta}}$. L'assemblage est 1D en $z$ avec des B-matrices simplifiées.

**Application :** pressurisation d'un tube, rétraction axiale, torsion uniforme.

### 11.2 Tube mince ($h/R \ll 1$) — limite membrane

Quand $h/R \to 0$, la raideur de flexion $D_b = Eh^3/[12(1-\nu^2)] \to 0$.
La formulation se réduit à une **membrane cylindrique** avec uniquement les termes $A_m$, $B_m$, $C_m$.

Les équations d'équilibre de la membrane cylindrique de Laplace :
$$
\frac{N_{\theta\theta}}{R} = p \qquad \text{(pression interne } p\text{)}
$$
doivent être récupérées. Vérification : $N_{\theta\theta} = A_m\,\varepsilon_{\theta\theta} + B_m\,\varepsilon_{zz}$.

### 11.3 Tube de Cosserat (limite 1D, largeur $W \ll L$)

Pour un tube avec un seul anneau d'éléments en $\theta$ et $N$ éléments en $z$,
on retrouve la formulation d'une **poutre de Cosserat** :
- DOF : $6 \times (N+1)$ (nœuds sur la ligne médiane)
- Strains : $\xi_{tz}$ uniquement (strains le long de $z$)
- $D^{zz}_{\mathrm{cyl}}$ fournit la raideur EA, EI_y, EI_z, GJ

### 11.4 Coque hémisphérique / non-cylindrique

Pour une coque de révolution générale (sphère, paraboloïde, etc.),
les coordonnées $(\theta, s)$ avec $s$ = abscisse curviligne méridionale donnent :

$$
\xi_{0\theta} = \begin{pmatrix}0\\\sin\psi\\0\\r(\psi)\\ 0\\0\end{pmatrix},
\qquad
\xi_{0s} = \begin{pmatrix}0\\1\\0\\0\\1\\0\end{pmatrix}
$$

où $r(\psi)$ est le rayon en fonction de l'angle méridional $\psi$
et $\sin\psi$ remplace le terme de courbure dans $\phi_{0\theta}$.
La formulation SE(3) s'étend naturellement à ces géométries sans modification
de la forme faible — seule $g_0$ change.

---

## 12. Résumé des modifications par rapport au cas cartésien

| Quantité | Cartésien (plaque) | Cylindrique ($R$) |
|---|---|---|
| Domaine $\mathcal{A}$ | $[0,L_1]\times[0,L_2]$ | $[0,2\pi]\times[0,L]$ |
| $\varphi_0$ | plan XY | cylindre rayon $R$ |
| $R_0$ | $I_3$ | $[\hat{e}_\theta|\hat{e}_z|\hat{e}_r](\theta)$ |
| $\xi_{01}$ | $[0,0,0,1,0,0]^T$ | $[0,1,0,R,0,0]^T$ |
| $\xi_{02}$ | $[0,0,0,0,1,0]^T$ | $[0,0,0,0,1,0]^T$ |
| $j_0$ | $1$ | $R$ |
| $D^{11}$ | $D_{\mathrm{Cart}}$ | $D_{\mathrm{Cart}}/R^2$ |
| $D^{22}$ | $D_{\mathrm{Cart}}$ | $D_{\mathrm{Cart}}$ |
| $D^{12}$ | $0$ (incorrect) | $D_{\mathrm{Cart}}^{\mathrm{Poisson}}/R$ |
| $\mathrm{ad}_{\xi_{01}}$ | blocs nuls | blocs avec $1/R$ encodé |
| Carte tangente $J_R^{-1}$ | optionnelle | **obligatoire** |

---

## 13. Notes FEM et implémentation

### Maillage

Le maillage naturel en $(\theta, z)$ est une grille de quadrilatères Q4 :

```
z
↑  ┌──┬──┬──┬──┐
   │  │  │  │  │
   ├──┼──┼──┼──┤
   │  │  │  │  │
   └──┴──┴──┴──┘ → θ (0 à 2π, périodique)
```

La **périodicité** en $\theta$ impose des conditions de raccordement :
nœuds à $\theta = 0$ et $\theta = 2\pi$ sont **identiques** — à traiter
par une contrainte de periodicité ou en dupliquant les nœuds avec contrainte.

### Paramètre `cylinderRadius` à ajouter

Dans `CosseratShellForceField`, ajouter :

```cpp
sofa::Data<SReal> d_cylinderRadius;   // R [m], 0 = mode plaque plate
```

Logique dans `init()` :
- Si `d_cylinderRadius == 0` : formulation cartésienne (actuelle)
- Si `d_cylinderRadius > 0` : formulation cylindrique — calcul de $g_0(\theta,z)$,
  $j_0 = R$, et mise à l'échelle des $D^{\alpha\beta}_{\mathrm{cyl}}$

### Ordre des nœuds par élément

En cylindrique, l'ordre standard (counter-clockwise en vue $(\theta, z)$) :
```
Nœuds : (θ-Δθ/2, z-Δz/2) → (θ+Δθ/2, z-Δz/2) → (θ+Δθ/2, z+Δz/2) → (θ-Δθ/2, z+Δz/2)
```

Le jacobien d'élément est $j_e = R_e\,\Delta\theta\,\Delta z/4$ avec $R_e$ le rayon local.

### Validation recommandée

1. **Test de pression interne uniforme** : vérifier $N_{\theta\theta} = pR$, $N_{zz} = pR/2$
   (formules de membrane de Laplace)
2. **Test de flexion axiale** : tige creuse soumise à un moment aux extrémités,
   comparer avec la solution analytique (rigidité de flexion $EI = E\pi R^3 h$)
3. **Test de torsion** : angle de torsion $\phi = TL/(GJ)$ avec $J = 2\pi R^3 h$
4. **Comparaison avec les benchmarks du papier** (Section 6) en géométrie cylindrique

---

*Référence : `cosserat_shell_formulation.md` + `formulation_review.md` + `cosserat_shell_v0.pdf`.*  
*Document créé le 2026-05-26.*
