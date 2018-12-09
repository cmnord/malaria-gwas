import dendropy
from ete3 import PhyloTree

def dendro():
    species_tree = dendropy.Tree.get_from_path("./data/receptor_genes/species_tree.nex", schema="nexus")
    
    proteins = ['band3', 'basigin', 'duffy', 'epha2', 'gypc', 'scarb1', 'tfr1']
    gene_trees = []
    for prot in proteins:
        tree = dendropy.Tree.get_from_path("./data/receptor_genes/" + prot + "/" + prot + "_all_primates.nex", schema="nexus")
        gene_trees.append(tree)
    
    contained = [0]*len(species_tree.taxon_namespace)
    for i in  range(len(species_tree.taxon_namespace)):
        species = species_tree.taxon_namespace[i]
        for tree in gene_trees:
            if str(species) in [str(x) for x in tree.taxon_namespace]:
                contained[i] += 1
    
    mapping = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(containing_taxon_namespace=species_tree.taxon_namespace, num_contained=contained)
    
    dendropy.model.reconcile.ContainingTree(species_tree, mapping.domain_taxon_namespace, mapping, contained_trees=species_tree)

    #tree = dendropy.Tree.get_from_path("./data/receptor_genes/sample_nexus", schema="nexus")
    #tree.write(path="./data/receptor_genes/sample_nexus", schema="nexus")

    species_tree.print_plot()

def phylo():
    proteins = ['band3', 'basigin', 'duffy', 'epha2', 'gypc', 'scarb1', 'tfr1']
    gene_trees = []
    for prot in proteins:
        tree = PhyloTree(alignment="./data/receptor_genes/" + prot + "/" + prot + "_all_primates.fa", alg_format="fasta")
        gene_trees.append(tree)
        print(gene_trees)

    species_tree = PhyloTree("(((galeopterus variegatus, nomascus leucogenys, (homo sapiens, (pongo abelii, pongo pygmaeus),(pan troglodytes, pan paniscus), gorilla), (piliocolobus tephrosceles, colobus angolensis palliatus, (rhinopithecus roxellana, rhinopithecus bieti), mandrillus leucophaeus, theropithecus gelada, papio anubis, (macaca nemestrina, macaca mulatta, macaca fascicularis), (chlorocebus sabaeus, chlorocebus aethiops), cercocebus atys), plecturocebus moloch, (aotus nancymaae, aotus trivirgatus)), saguinus labiatus, saimiri boliviensis, cebus capucinus), callithrix jacchus);")
    print(species_tree)
    recon_trees = []
    for tree in gene_trees:
        recon_trees.append(tree.reconcile(species_tree)[0])

    for tree in recon_trees:
        print(tree)

def main():
    dendro()

if __name__ == "__main__":
    main()