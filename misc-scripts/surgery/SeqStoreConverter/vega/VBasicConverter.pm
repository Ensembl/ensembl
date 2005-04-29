package SeqStoreConverter::vega::VBasicConverter;

use strict;
use warnings;
use SeqStoreConverter::BasicConverter;
use vars qw(@ISA);
@ISA = qw(SeqStoreConverter::BasicConverter);

sub remove_supercontigs {
    my $self = shift;
    
    my $target = $self->target();
    my $dbh    = $self->dbh();
    $self->debug("Vega specific - removing supercontigs from $target");

    $dbh->do("DELETE FROM $target.meta ". 
             "WHERE meta_value like '%supercontig%'");

    $dbh->do("DELETE FROM $target.coord_system ".
             "WHERE name like 'supercontig'");
    
    $dbh->do("DELETE $target.a ".
             "FROM $target.assembly a, $target.seq_region sr ". 
             "WHERE sr.coord_system_id = 2 ".
             "and a.asm_seq_region_id = sr.seq_region_id");

    $dbh->do("DELETE FROM $target.seq_region ".
             "WHERE coord_system_id = 2");
}

sub copy_other_vega_tables {
    my $self = shift;
    $self->copy_tables(
        # vega tables
        "gene_synonym",
        "transcript_info",
        "current_gene_info",
        "current_transcript_info",
        "author",
        "gene_name",
        "transcript_class",
        "gene_remark",
        "gene_info",
        "evidence",
        "transcript_remark",
        "clone_remark",
        "clone_info",
        "clone_info_keyword",
        "clone_lock",
        "assembly_tag",
    );
    $self->copy_current_clone_info;
}

sub copy_current_clone_info {
    my $self=shift;
    my $source = $self->source();
    my $target = $self->target();
    my $sth = $self->dbh()->prepare
        ("INSERT INTO $target.current_clone_info(clone_id,clone_info_id) SELECT * FROM $source.current_clone_info");
    $sth->execute();
    $sth->finish();    
}

sub update_genscan {
    my $self = shift;
    $self->debug("Vega specific - updating analysis name for Genscans");
    my $target = $self->target();
    my $sth = $self->dbh()->prepare
        ("UPDATE $target.analysis set logic_name = 'Vega_Genscan' where logic_name = 'Genscan'");
    $sth->execute();
    $sth->finish();    
}       

sub update_clone_info {
    my $self = shift;
    return;
}

sub copy_internal_clone_names {
    my $self = shift;
    return;
}

sub copy_assembly_exception {
    my $self = shift;

    # copy assembly_exception table
    $self->debug('Vega specific - copying assembly_exception table');
    $self->copy_tables('assembly_exception');

    my $source = $self->source();
    my $target = $self->target();
    my $dbh    = $self->dbh();

    # fix seq_region_id in assembly_exception
    $self->debug('Vega specific - Updating seq_region_id in assembly_exception table');
    $dbh->do(qq(
        UPDATE $target.assembly_exception, $target.tmp_chr_map
        SET assembly_exception.seq_region_id = tmp_chr_map.new_id
        WHERE assembly_exception.seq_region_id = tmp_chr_map.old_id
    ));
    $dbh->do(qq(
        UPDATE $target.assembly_exception, $target.tmp_chr_map
        SET assembly_exception.exc_seq_region_id = tmp_chr_map.new_id
        WHERE assembly_exception.exc_seq_region_id = tmp_chr_map.old_id
    ));

    # fix seq_region.length if necessary (this is the case if you have an
    # assembly_exception at the end of a chromosome)
    my $sth1 = $dbh->prepare(qq(
        UPDATE $target.seq_region SET length = ? WHERE seq_region_id = ?
    ));
    my $sth2 = $dbh->prepare(qq(
        SELECT
                sr.seq_region_id,
                sr.length,
                max(ae.seq_region_end)
        FROM
                $target.seq_region sr,
                $target.assembly_exception ae
        WHERE   sr.seq_region_id = ae.seq_region_id
        GROUP BY ae.seq_region_id
    ));
    $sth2->execute;
    while (my ($sr_id, $sr_length, $max_ae_length) = $sth2->fetchrow_array) {
        if ($max_ae_length > $sr_length) {
            $self->debug("  Updating seq_region.length for $sr_id (old $sr_length, new $max_ae_length)");
            $sth1->execute($max_ae_length, $sr_id);
        }
    }
}



