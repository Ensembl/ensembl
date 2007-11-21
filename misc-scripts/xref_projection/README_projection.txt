Gene name and GO term projection in Ensembl
===========================================


Gene names (also known as display xrefs) are projected between species in Ensembl using homology information obtained from Ensembl Compara. Details of the way the homologies are obtained may be found at the following URL:

http://www.ensembl.org/info/data/compara/homology_method.html


Gene name projection
--------------------

In Ensembl, most genes have a name (e.g. BRCA2) which is assigned from one of the external identifiers (xrefs) associated with that gene, based on a priority system. The display name of the gene is also known as the display xref. Because of the comparatively large amount of external information available for human and mouse, a much larger proportion of human and mouse genes have names than do genes in other species. For this reason we "project" the names between species.
 
The current projections which are done are:

Human to all mammals (1:1 orthologues only)
Mouse to rat (1:1 orthologues only)
Human to all fish (using 1-many ortholgues as well)

The projection works as follows:

  Fetch all homologies from Compara
  For each homology, extract "from" (e.g. human) and "to" (e.g. dog) genes
  If the "to" gene has no display xref, or a low-priority one (see below): 
    Assign "from" gene's name and description to "to" gene
    Set status of "to" gene to KNOWN_BY_PROJECTION

If the "to" gene has an existing display xref, this is not overwritten, except in the case where the "to" gene has a name derived from a RefSeq_dna_predicted or RefSeq_peptide_predicted xref, and the "from" gene has a name derived from an HGNC (in the case of human) or MGI (in the case of mouse) xref.

Note that only homologies of type "ortholog_one2one" and "apparent_ortholog_one2one" are used from Compara. Other types of homology (e.g. one to many) are not used for the projection. This is the case for both the gene name and GO term projection.

As of release 49 (February 2008), only HGNC symbols will be projected from human to the other species. Projections from mouse will continue to use gene names from other sources, e.g. UniProtKB/Swissprot.


GO term projection
------------------

As with gene names, there are many more GO (Gene Ontology) annotations for human, and to a lesser extent, mouse, than other species. For this reason we also project GO terms between species. Not all evidence types are projected; currently only the following evidence types are projected:  IDA, IEP, IGI, IMP, IPI.

The following projections are done currently:

Human to all mammals
Mouse to all mammals
Drosophila to anopheles

All GO term projections use 1:1 orthologues only.

The algorithm for the projection is similar to that for gene names:

  Fetch all homologies from Compara
  For each homology, extract "from" (e.g. human) and "to" (e.g. dog) genes
    Assign "from" gene's GO terms to "to" gene
