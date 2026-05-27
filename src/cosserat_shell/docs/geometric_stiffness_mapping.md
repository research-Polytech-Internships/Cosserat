# Raideur géométrique du mapping Cosserat

> Auteur : Y. Adagolodjo (DEFROST / INRIA)
> Date : 2026-05-26
> Contexte : analyse de `DiscreteCosseratMapping.h/.inl/.cpp` et lien avec `ShellElementComputerFull.h`

---

## 1. État des lieux dans le code

### Architecture du composant

Le mapping est implémenté en trois fichiers :

| Fichier | Contenu |
|---|---|
| `DiscreteCosseratMapping.h` | Déclarations, dont `applyDJT` **déclarée vide** |
| `DiscreteCosseratMapping.inl` | Implémentation générique (template `TIn1, TIn2, TOut`) |
| `DiscreteCosseratMapping.cpp` | Spécialisations explicites `Vec6Types` + enregistrement factory |

Deux instanciations sont enregistrées dans SOFA :
- `Vec3Types, Rigid3Types, Rigid3Types` — strains angulaires uniquement (3 DOF)
- `Vec6Types, Rigid3Types, Rigid3Types` — strains complets angulaires + linéaires (6 DOF)

### Méthodes et état par spécialisation

| Méthode SOFA | Rôle | Vec3 (`.inl`) | Vec6 (`.cpp`) |
|---|---|---|---|
| `apply()` | $x = \phi(q)$ | ✅ | ✅ (hérité du `.inl`) |
| `applyJ()` | $\dot{x} = J\dot{q}$ | ✅ | ✅ spécialisé |
| `applyJT()` (forces) | $f_q = J^T f_x$ | ✅ | ✅ spécialisé |
| `applyJT()` (contraintes) | pullback contraintes | ✅ | ✅ spécialisé |
| `applyDJT()` | **raideur géométrique** | ❌ **corps vide `{}`** | ❌ **corps vide `{}`** |

`applyDJT` est **déclarée** dans le header mais avec un corps vide — SOFA accepte
le composant sans erreur de compilation mais $K_G^{\mathrm{mapping}} = 0$ de facto.

### Comparaison Vec3Types vs Vec6Types

Les deux spécialisations partagent la même structure de Jacobienne mais diffèrent
dans le sous-espace de strain activé :

| Aspect | Vec3 (`.inl`) | Vec6 (`.cpp`) |
|---|---|---|
| DOFs strains `TIn1` | `Vec3` — strains **angulaires** uniquement (3 composantes) | `Vec6` — strains **complets** (6 composantes : angulaires + linéaires) |
| `matB_trans` dans `applyJT` | `Mat3x6` — sélecteur 3×6 (3 premières lignes actives) | `Mat6x6 = I_6` — identité 6×6 (toutes composantes actives) |
| `node_Xi_dot` dans `applyJ` | `Vec3` (`in1_vel[i-1]`, 3 composantes) | `Vec6` (6 composantes — voir bug §1.3) |
| Output forces `TOut` | `Rigid3Types` (6 DOF) | `Rigid3Types` (6 DOF) |
| Référence de base `TIn2` | `Rigid3Types` | `Rigid3Types` |
| `apply()` | hérité du `.inl` | hérité du `.inl` |
| `applyJ()` | générique `.inl` | **spécialisé** `.cpp` |
| `applyJT()` forces | générique `.inl` | **spécialisé** `.cpp` |
| `applyJT()` contraintes | générique `.inl` | **spécialisé** `.cpp` |
| `applyDJT()` | ❌ corps vide | ❌ corps vide |

### Bug identifié dans `applyJ` Vec6Types

**Fichier** : `DiscreteCosseratMapping.cpp`, dans
`applyJ<Vec6Types, Rigid3Types, Rigid3Types>`, boucle de calcul des vitesses nodales.

**Code actuel (bugué)** :
```cpp
Vec6 node_Xi_dot;
for (unsigned int u = 0; u < 6; u++)
    node_Xi_dot(i) = in1_vel[i-1][u];   // ← BUG : i est la variable du
                                          //          boucle externe (nœuds) !
```

**Code corrigé** :
```cpp
Vec6 node_Xi_dot;
for (unsigned int u = 0; u < 6; u++)
    node_Xi_dot(u) = in1_vel[i-1][u];   // ← CORRECT : u est la composante
```

**Conséquence** : dans la boucle `u = 0..5`, `node_Xi_dot(i)` est écrasé
5 fois (avec `i` fixé au rang du nœud courant). Toutes les autres composantes
de `node_Xi_dot` restent à zéro. La vitesse nodale `eta_node_i` est donc
calculée avec un vecteur de strain-velocity incorrect — seul `node_Xi_dot(i)`
prend la valeur de `in1_vel[i-1][5]`, toutes les autres composantes sont nulles.

Ce bug est présent dans la spécialisation Vec6Types uniquement (la version Vec3
dans `.inl` utilise directement `in1_vel[i-1]` sans boucle manuelle).

**Commit suggéré** :
```
[cosserat-shell] fix applyJ Vec6 : node_Xi_dot(i) -> node_Xi_dot(u)
```

---

## 2. Définition exacte de `applyDJT`

### Contexte : décomposition de la raideur totale

Pour un mapping $\phi : q \to x = \phi(q)$ de Jacobienne $J(q) = \partial\phi/\partial q$,
le principe des travaux virtuels donne le pullback des forces :

$$f_q = J(q)^T\,f_x$$

La raideur tangente totale dans l'espace parent $q$ est :

$$K_{\mathrm{total}} = \frac{\partial f_q}{\partial q}
= \underbrace{J^T K_x J}_{\text{terme matériel}} + \underbrace{\frac{\partial (J^T)}{\partial q} f_x}_{\text{terme géométrique}}$$

**`applyDJT` calcule uniquement le terme géométrique** :

$$\mathrm{applyDJT}(\delta q) = \left.\frac{d}{d\varepsilon}\right|_{\varepsilon=0} J(q + \varepsilon\,\delta q)^T \cdot f_x$$

### Ce que calcule exactement `applyDJT`

Ce n'est **pas** :
- $\dfrac{\partial J^T}{\partial q}$ seul — ce serait un tenseur d'ordre 3, non applicable directement
- $\dfrac{\partial (J^T f_x)}{\partial q}$ au sens total — cela ferait aussi varier $f_x$ avec $q$

C'est **la dérivée directionnelle du produit $J^T f_x$ en maintenant $f_x$ constant** :

$$\left[\mathrm{applyDJT}(\delta q)\right]_k
= \sum_{i,j} \frac{\partial J_{ij}(q)}{\partial q_k}\,\delta q_k \cdot (f_x)_j$$

En d'autres termes : on fait varier **$J^T$ seul** (la géométrie du mapping change avec $q$),
avec le vecteur de forces enfant $f_x$ fixé à sa valeur courante.

### Interprétation physique

Ce terme encode le fait que **la direction du pullback des forces change** quand la
configuration $q$ évolue. Pour un mapping SE(3) (Cosserat), cela correspond aux
termes centrifuges / gyroscopiques qui apparaissent lorsque la poutre tourne
significativement entre deux itérations Newton.

---

## 3. Conséquences de l'absence de `applyDJT`

| Régime | Impact |
|---|---|
| Petites déformations | Négligeable — $\partial J^T/\partial q \approx 0$ |
| Grandes rotations (> ~30°) | Convergence Newton dégradée (plus d'itérations) |
| Chargements importants / post-flambage | Possible divergence de Newton |
| Dynamique implicite avec grands pas de temps | Instabilité potentielle |

---

## 4. Cohérence avec la raideur géométrique du ForceField

La raideur géométrique totale du système se décompose en **deux contributions** :

```
K_total = K_ff + K_mapping
```

| Contribution | Source | État actuel | État Full |
|---|---|---|---|
| $K_G^{\mathrm{ff}}$ (ForceField) | `ShellElementComputer.h` | ⚠️ $-\mathrm{ad}^T$ (approximatif) | ✅ `ShellElementComputerFull.h` |
| $K_G^{\mathrm{mapping}}$ (Mapping) | `DiscreteCosseratMapping.inl` | ❌ absent (`applyDJT` manquant) | ⏳ à implémenter |

La priorité d'implémentation est :
1. **`ShellElementComputerFull.h`** — corrige $K_G^{\mathrm{ff}}$ ✅ (fait)
2. **`applyDJT` dans `DiscreteCosseratMapping`** — ajoute $K_G^{\mathrm{mapping}}$ ⏳

---

## 5. Forme de `applyDJT` pour un mapping SE(3)

Pour le mapping Cosserat $\phi : (\xi, g_0) \to g_{\mathrm{frame}}$
avec $g_{\mathrm{frame}} = g_0 \cdot \exp(\xi_1) \cdots \exp(\xi_k) \cdot \exp(\xi_{\mathrm{frame}})$,
la dérivée $\frac{\partial J^T}{\partial q}$ fait intervenir les **adjoints de SE(3)**.

Schématiquement, pour chaque frame $i$ :

$$\frac{\partial}{\partial \xi_k}\bigl[J(\xi)^T f_i\bigr]
= \mathrm{Ad}_{g_{k \to i}}^T\,\left[\mathrm{ad}_{\tilde\xi_k}^T\,f_i\right]$$

où :
- $\mathrm{Ad}_{g_{k \to i}}$ est l'adjoint de SE(3) du frame $k$ au frame $i$
- $\mathrm{ad}_{\tilde\xi_k}$ est la matrice adjoint (petite) de la strain twist courante $\xi_k$
- $f_i$ est la force appliquée au frame $i$ (tenue constante)

Ce terme est structurellement analogue à la raideur géométrique du ForceField :
$$K_G^{\mathrm{ff}} = B^T\,\mathrm{ad}_{S^\alpha}\,B$$
mais s'applique au niveau du mapping entre espaces de configuration.

---

## 6. Signification de $T$ dans la formulation

### Rôle : tangente de l'exponentielle SE(3)

$T_k$ est la **dérivée de l'application exponentielle** $\exp : \mathfrak{se}(3) \to SE(3)$
évaluée à la configuration courante $\xi_k L_k$. Elle répond à la question :

> *Si je perturbe le strain $\xi_k$ d'un incrément $\delta\xi_k$, quel twist infinitésimal $v$
> cela produit-il dans l'algèbre de Lie ?*

$$v_k = T_k \cdot \delta\xi_k \quad \in \mathfrak{se}(3) \cong \mathbb{R}^6$$

### Définition formelle

Pour une section de longueur $L_k$ et de strain courant $\xi_k$ :

$$T_k = \left.\frac{d}{d\varepsilon}\right|_{\varepsilon=0}
\exp\!\bigl((\xi_k + \varepsilon\,\delta\xi_k)\,L_k\bigr)
\cdot \exp(-\xi_k\,L_k)$$

C'est la **Jacobienne droite de l'exponentielle** (right Jacobian), notée $J_R$ en
littérature, mise à l'échelle par $L_k$ :

$$T_k = J_R(\xi_k L_k) \cdot L_k \in \mathbb{R}^{6\times 6}$$

Elle a une forme close (série de Lie) :

$$T_k = L_k\,\mathbf{I}_6
+ \frac{L_k^2}{2}\,\mathrm{ad}_{\xi_k}
+ \frac{L_k^3}{6}\,\mathrm{ad}_{\xi_k}^2 + \ldots$$

Pour $\|\xi_k\| \to 0$ (poutre droite non déformée) : $T_k \to L_k\,\mathbf{I}_6$.

### Dans le code

| Symbole | Variable SOFA | Signification |
|---|---|---|
| $T_k^{\text{node}}$ | `m_nodesTangExpVectors[k+1]` | Tangente section $k$ à longueur complète $L_k$ |
| $T_s^{\text{frame}}$ | `m_framesTangExpVectors[s]` | Tangente section de la frame $s$ à longueur partielle $s_{\text{rel}}$ |

Calculé par `computeTangExpImplementation` dans `BaseCosseratMapping.inl`.

### Pourquoi $T$ apparaît dans `applyDJT`

Dans `applyJT`, la force sur le strain $k$ est :

$$f_k = B^T \underbrace{T_k^T}_{\text{pullback tangent}} \underbrace{F_{\text{tot},k}}_{\text{force transportée}}$$

Quand on perturbe $\xi_k \to \xi_k + \varepsilon\,\delta\xi_k$, la map exponentielle change
et donc l'adjoint $\mathrm{coAd}(E_k)$ change aussi. La variation est :

$$\delta\bigl(\mathrm{coAd}(E_k^{-1})\bigr) \cdot F
= \mathrm{ad}_{T_k\,\delta\xi_k}^T \cdot \mathrm{coAd}(E_k^{-1})\,F$$

C'est pourquoi $T_k \delta\xi_k$ apparaît comme argument de $\mathrm{ad}(\cdot)^T$ :
c'est le **twist infinitésimal** engendré par la perturbation du strain, et c'est lui qui
"tord" la direction du pullback des forces.

### En résumé

$$\boxed{v_k = T_k\,\delta\xi_k \in \mathfrak{se}(3)}$$

$v_k$ est le **twist de vitesse local** à la section $k$ correspondant à l'incrément de
strain $\delta\xi_k$. Il mesure "dans quelle direction la géométrie du mapping bouge"
quand on perturbe la configuration — c'est précisément ce qui pilote la correction
géométrique $\mathrm{ad}(v_k)^T \cdot F$.

---

## 8. Résumé

| Question | Réponse |
|---|---|
| `applyDJT` calcule la dérivée de quoi ? | De $J^T f_x$ par rapport à $q$, avec $f_x$ **fixé** |
| C'est $\partial J^T / \partial q$ seul ? | Non — c'est la dérivée directionnelle du **produit** $J^T f_x$ |
| C'est $\partial (J^T f) / \partial q$ total ? | Non — $f$ ne varie pas (tenu constant à sa valeur courante) |
| Implémenté dans `DiscreteCosseratMapping` ? | ❌ Non — à implémenter |
| Impact pour petites déformations ? | Négligeable |
| Impact pour grandes rotations ? | Convergence Newton dégradée / divergence possible |

---

*Références : `formulation_review.md` §5, `ShellElementComputerFull.h` §K_G*
*Document créé le 2026-05-26.*
