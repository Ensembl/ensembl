=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package SeqStoreConverter::vega::DanioRerio;

use strict;
use warnings;

use SeqStoreConverter::DanioRerio;
use SeqStoreConverter::vega::VBasicConverter;
use vars qw(@ISA);

@ISA = qw( SeqStoreConverter::DanioRerio SeqStoreConverter::vega::VBasicConverter);

sub copy_internal_clone_names {
    my $self = shift;

    my $target = $self->target();
    my $source = $self->source();
    my $dbh    = $self->dbh();
    $self->debug("Vega danio specific - copying internal clone names to seq_region_attrib");

#get id for 'fpc_clone_id' attribute

    $dbh->do("INSERT INTO $target.attrib_type (code,name,description)".
            "values ('fpc_clone_id','fpc clone','clone id used for linking to Zebrafish webFPC')");

    my ($attrib_id) = $dbh->selectrow_array("Select attrib_type_id from $target.attrib_type where code = 'fpc_clone_id'");
    warn "No attrib id found\n" unless defined($attrib_id);

#get clone details
    my $select1_sth = $dbh->prepare
        ("SELECT seq_region_id, name from $target.seq_region where coord_system_id = 3;");
    $select1_sth->execute();
    my ($seq_region_id, $embl_name);
    $select1_sth->bind_columns(\$seq_region_id, \$embl_name);

    my $clone_name;
    my $select2_sth = $dbh->prepare("select name from $source.clone where embl_acc= ?");

    my $insert_sth = $dbh->prepare("insert into $target.seq_region_attrib values (?,$attrib_id,?)");

    while ($select1_sth->fetch()) {
        $embl_name =~ s/([\d\w]+).*/$1/;
        $select2_sth->bind_param(1,$embl_name);
        $select2_sth->execute;
        $insert_sth->bind_param(1,$seq_region_id);
        while (my ($clone_name) = $select2_sth->fetchrow_array()) {
            $insert_sth->bind_param(2,$clone_name);
            $insert_sth->execute();
        }
    }
}

sub update_clone_info {
    my $self = shift;
    my $target_cs_name = shift;

    my $target = $self->target();
    my $source = $self->source();
    my $dbh    = $self->dbh();

    # clone_info, current_clone_info
    $self->debug("Vega Danio_specific - Transforming clone_id into seq_region_id for clone_info and current_clone_info");

    foreach my $table_name ('clone_info','current_clone_info') {
        my $select_st1 = 
            "SELECT ctg.name, ctg.clone_id " .
            "FROM   $source.contig ctg, $source.$table_name ci " .
            "WHERE  ctg.clone_id = ci.clone_id " .
            "AND ctg.name not like 'ctg%' " .
            "AND ctg.name not like 'NA%'";

        my $query_results1 = $dbh->selectall_arrayref($select_st1);

        my $i = 0;
        foreach my $contig_name (@$query_results1) {
            my $embl_acc = $contig_name->[0];
            my $select_st2 = 
                "SELECT sr.seq_region_id " .
                "FROM $target.seq_region sr " . 
                "WHERE sr.name = '$embl_acc'";
            my @query_results2 = $dbh->selectrow_array($select_st2);
            push @{$query_results1->[$i]},@query_results2;
            $i++;
        }

        foreach my $clone (@$query_results1) {
            my $seq_reg_id = $clone->[2];
            my $clone_id = $clone->[1];

            my $update_query = 
                "UPDATE $target.$table_name " .
                "SET clone_id = '$seq_reg_id' " .
                "WHERE clone_id = '$clone_id'";
            $dbh->do($update_query);
        }
        my $alter_struct_1 = 
            "ALTER table $target.$table_name " .
            "CHANGE clone_id seq_region_id int(10) not null";
        my $alter_struct_2 = 
            "ALTER table $target.$table_name " .
            "add unique index (seq_region_id)";
        $dbh->do($alter_struct_1);
        $dbh->do($alter_struct_2);
    }

    # assembly_tag
    $self->debug("Vega Danio_specific - Transforming contig_id into seq_region_id for assembly_tag");

    # first remove orphans from assembly_tag table (i.e. entries pointing to 
    # non-existing contigs)
    my $numrows = $dbh->do(qq(
        DELETE at
        FROM $source.assembly_tag at
        LEFT JOIN $source.contig c ON c.contig_id = at.contig_id
        WHERE c.contig_id IS NULL
    ));
    $self->debug("  Deleted $numrows orphans from assembly_tag");

    my $select_st3 = 
        "SELECT ctg.name, ctg.contig_id " .
        "FROM   $source.contig ctg, $source.assembly_tag at " .
        "WHERE  ctg.contig_id = at.contig_id " .
        "AND ctg.name not like 'ctg%' " .
        "AND ctg.name not like 'NA%'";

    my $query_results3 = $dbh->selectall_arrayref($select_st3);

    my $j = 0;
    foreach my $contig_name (@$query_results3) {
        my $embl_acc = $contig_name->[0];
        my $select_st4 = 
            "SELECT sr.seq_region_id " .
            "FROM $target.seq_region sr " . 
            "WHERE sr.name = '$embl_acc'";
        my @query_results4 = $dbh->selectrow_array($select_st4);
        push @{$query_results3->[$j]}, @query_results4;
        $j++;
    }

    foreach my $contig (@$query_results3) {
        my $seq_reg_id = $contig->[2];
        my $contig_id = $contig->[1];

        my $update_query = 
            "UPDATE $target.assembly_tag " .
            "SET contig_id = '$seq_reg_id' " .
            "WHERE contig_id = '$contig_id'";
        $dbh->do($update_query);
    }

    $dbh->do("  ALTER TABLE $target.assembly_tag
                CHANGE contig_id seq_region_id int(10) UNSIGNED NOT NULL");

    $dbh->do("  ALTER TABLE $target.assembly_tag
                CHANGE contig_start seq_region_start int(10)");

    $dbh->do("  ALTER TABLE $target.assembly_tag
                CHANGE contig_end seq_region_end int(10)");

    $dbh->do("  ALTER TABLE $target.assembly_tag
                CHANGE contig_strand seq_region_strand tinyint(1)");

}

1;
