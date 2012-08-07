package XrefMapper::neurospora_crassa;

use  XrefMapper::BasicMapper;

use vars qw(@ISA);

@ISA = qw(XrefMapper::BasicMapper);


sub gene_display_xref_sources {
    my $self     = shift;

    my @list = qw(
                 Uniprot_genename
               );
    
    my %ignore;
    
    return [\@list,\%ignore];
}

sub transcript_display_xref_sources {
    my $self     = shift;

    my @list = qw(
                 Uniprot_genename
               );
    
    my %ignore;
    
    return [\@list,\%ignore];
}


1;
