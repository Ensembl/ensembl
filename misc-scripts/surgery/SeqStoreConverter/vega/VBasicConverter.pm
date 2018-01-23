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
        "assembly_tag",
   );
    eval { $self->copy_current_clone_info; };
	warn $@ if $@;
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

# reset gene, transcript, gene_description and external_db tables back to 30
sub back_patch_schema {
	my $self = shift;
	$self->debug ("Patching gene, transcript, gene_description and external_db tables back to sch-30");
	my $target = $self->target;
	my $dbh = $self->dbh;
	$dbh->do("DROP TABLE $target.gene");
	$dbh->do( qq(CREATE TABLE $target.gene (
                 `gene_id` int(10) unsigned NOT NULL auto_increment,
                 `type` varchar(40) NOT NULL default '',
                 `analysis_id` int(11) default NULL,
                 `seq_region_id` int(10) unsigned NOT NULL default '0',
                 `seq_region_start` int(10) unsigned NOT NULL default '0',
                 `seq_region_end` int(10) unsigned NOT NULL default '0',
                 `seq_region_strand` tinyint(2) NOT NULL default '0',
                 `display_xref_id` int(10) unsigned default NULL,
                  PRIMARY KEY  (`gene_id`),
                  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
                  KEY `xref_id_index` (`display_xref_id`),
                  KEY `analysis_idx` (`analysis_id`)
                  ) ENGINE=MyISAM DEFAULT CHARSET=latin1
                ));
	$dbh->do("DROP TABLE $target.transcript");
	$dbh->do( qq(CREATE TABLE $target.transcript (
                 `transcript_id` int(10) unsigned NOT NULL auto_increment,
                 `gene_id` int(10) unsigned NOT NULL default '0',
                 `seq_region_id` int(10) unsigned NOT NULL default '0',
                 `seq_region_start` int(10) unsigned NOT NULL default '0',
                 `seq_region_end` int(10) unsigned NOT NULL default '0',
                 `seq_region_strand` tinyint(2) NOT NULL default '0',
                 `display_xref_id` int(10) unsigned default NULL,
                 PRIMARY KEY  (`transcript_id`),
                 KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
                 KEY `gene_index` (`gene_id`),
                 KEY `xref_id_index` (`display_xref_id`)
                 ) ENGINE=MyISAM DEFAULT CHARSET=latin1
                ));
	$dbh->do( qq(CREATE TABLE $target.gene_description (
                 `gene_id` int(10) unsigned NOT NULL default '0',
                 `description` text,
                 PRIMARY KEY  (`gene_id`)
                 ) ENGINE=MyISAM DEFAULT CHARSET=latin1
                ));
	$dbh->do("DROP TABLE $target.external_db");
	$dbh->do( qq(CREATE TABLE $target.external_db (
                 `external_db_id` int(11) NOT NULL default '0',
                 `db_name` varchar(100) NOT NULL default '',
                 `release` varchar(40) NOT NULL default '',
                 `status` enum('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') NOT NULL default 'KNOWNXREF',
                 PRIMARY KEY  (`external_db_id`)
                 ) ENGINE=MyISAM DEFAULT CHARSET=latin1
                ));
}




