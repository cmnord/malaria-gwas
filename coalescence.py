import dendropy
import pdb
from ete3 import PhyloTree, TreeStyle

SPECIES_TREE_FILE = 'data/receptor_genes/{}/species_tree.{}'
GENE_TREE_FILE = 'data/receptor_genes/{}/{}_all_primates.{}'


def reconcile_dendropy(protein):
    species_tree = dendropy.Tree.get_from_path(
        SPECIES_TREE_FILE.format(protein, 'nex'), schema="nexus")
    gene_tree = dendropy.Tree.get_from_path(
        GENE_TREE_FILE.format(protein, protein, 'nex'), schema="nexus")

    mapping = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        species_tree.taxon_namespace, num_contained=1,
        contained_taxon_label_fn=lambda taxon, index: taxon)

    dendropy.model.reconcile.ContainingTree(
        species_tree, mapping.domain_taxon_namespace, mapping, contained_trees=gene_tree)

    # tree = dendropy.Tree.get_from_path("./data/receptor_genes/sample_nexus", schema="nexus")
    # tree.write(path="./data/receptor_genes/sample_nexus", schema="nexus")

    species_tree.print_plot()


def reconcile_all_proteins():
    proteins = ['band3', 'basigin', 'duffy', 'epha2', 'gypc', 'scarb1', 'tfr1']
    for protein in proteins:
        print(protein)
        reconcile_etetoolkit(protein)


def reconcile_etetoolkit(protein):
    species_tree = PhyloTree(SPECIES_TREE_FILE.format(
        protein, 'nh'), format=1, sp_naming_function=lambda name: name)
    gene_tree = PhyloTree(GENE_TREE_FILE.format(
        protein, protein, 'nh'), format=1, sp_naming_function=lambda name: name)
    recon_tree, events = gene_tree.reconcile(species_tree)
    recon_tree.render("phylotree.png")


def main():
    reconcile_dendropy('band3')


if __name__ == "__main__":
    main()
