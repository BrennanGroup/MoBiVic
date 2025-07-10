# MoBiVic
MoBiVic (Mono- and Bi-functionalised VEHICLe) is a virtual database of aromatic heterocycles built by functionalising the VEHICLe dataset of [Pitt *et al.* (2009)](https://pubs.acs.org/doi/full/10.1021/jm801513z) with small, low molecular weight substituents.

## Repository Structure
The datasets are provided as `.json` files, with keys as RegIDs (each heterocycle's unique identifier in the database), and the values as the canonicalised SMILES strings. The following datasets are included in the Data directory (for details of their construction please see below):

| Name                          | Description                                                                                                                                                                                               |
|-------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `VEHICLe.json`                | The original VEHICLe database of Pitt *et al.* RegIDs correspond to the original RegIDs of the dataset                                                                                                    |
| `all_monofunctionalised.json` | The full dataset of monofunctionalised VEHICLe heterocycles. This is unfiltered, and so includes heterocycles with no remaining exit-vectors, and those that are likely to be highly reactive or unstable. |
| `all_bifunctionalised.json`   | The full dataset of bifunctionalised VEHICLe heterocycles. This is unfiltered, and so includes heterocycles with no remaining exit-vectors, and those that are likely to be highly reactive or unstable.  |
| `MoBiVic.json`    | The combined library of 546 271 unfunctionalised, monofunctionalised, and bifunctionalised heterocycles. This has been filtered according to the rules outlined below to remove potentially unstable molecules, and those without exit-vectors. |

Also included for reference is `functionalise_library.py`, the Python script used to generate these libraries.

## Construction
The VEHICLe database of [Pitt *et al.* (2009)](https://pubs.acs.org/doi/full/10.1021/jm801513z) was used as a starting point for the construction of the MoBiVic database. VEHICLe is a virtual database of 24 867 heterocycles conforming to the following criteria:

- Mono- or bi-cyclic systems of 5 or 6 membered rings (thus 6,6-,5,5-, and 5,6-bicycles are included);
- only containing the elements C, N, S, O, H;
- obeying Hückel's *4n+2* rule of aromaticity; and
- exocyclic C-O bonds exist as carbonyls rather than hydroxyls in all cases (keto or lactam tautomers only).

To functionalise the VEHICLe heterocycles, separate sets of substituents were chosen for bonding to carbon and nitrogen, based on a [2017 analysis of ChEMBL substituents by Astex](https://pubs.acs.org/doi/10.1021/acs.jmedchem.7b00809) and amended based on chemical experience. These are shown in panels A and B of the below figure.

For each heterocycle in VEHICLe, all the available exit-vectors (C-H or N-H bonds) were identified (for example, in the figure below indole has 7 exit-vectors: 6 C-H bonds and 1 N-H bond). Each substituent was then considered sequentially, and bonded to each available exit vector on the heterocycle in turn, resulting in a set of singly functionalised heterocycles where every available substituent is bonded to every available exit-vector exactly once. For the example of indole below this generates 44 monofunctionalised heterocycles (42 with substituents bonded to carbona and 2 with substituents bonded to nitrogen). 

<div style="text-align: center">
<img src="assets/substituents.png" alt="The substituents used to functionalise the VEHICLe database" width="400"/>
</div>

This process was then applied again to each of the monofunctionalised heterocycles to generate a library of bifunctionalised heterocycles. The number of heterocycles in each library is given below.

| Library            | Number of heterocycles |
|--------------------|------------------------|
| VEHICLe            | 24 867                 |
| Monofunctionalised | 336 816                |
| Bifunctionalised   | 486 220                |

### Filtering
Included within these libraries are heterocycles that are unlikely to ever stably exist, are too reactive to be useful in most applications, or that have no remaining exit-vectors and thus cannot be further functionalised (and are therefore not practically useful).

The combined library of VEHICLe, monofunctionalised, and bifunctionalised heterocycles was therefore curated to remove the molecules described above. Molecules were removed from the library if they met one or more of the following rules:

- The proportion of heavy atoms that are carbon is less than 20% (tetrazole was manually excluded from this filter);
- molecules with >3 aromatic nitrogens bonded together;
- molecules with 3 or more carbonyls;
- cyclic anhydrides; or
- cyclic thioesters.

Those with no remaining exit-vectors were also removed from the library.

This left a library of 546 271 heterocycles (containing both unfunctionalised, monofunctionalised, and bifunctionalised heterocycles), which we have termed the MoBiVic database.
