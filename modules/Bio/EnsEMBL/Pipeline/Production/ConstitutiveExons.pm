package Bio::EnsEMBL::Pipeline::Production::ConstitutiveExons;

use strict;
use warnings;

use base qw/Bio::EnsEMBL::Pipeline::Base/;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;

# default run method
# For a given databases of dbtype, flags exons which are constitutive
sub run {
  my ($self)      = @_;
  my $species     = $self->param('species');
  my $dbtype      = $self->param('dbtype');
  my $quick_check = $self->param('quick_check');
  my $ga = Bio::EnsEMBL::Registry->get_adaptor($species, $dbtype, 'Gene');

  if ($quick_check) {
	print "Carrying out quick check...\n";
	my $splice_count = $ga->dbc()->sql_helper()->execute_single_result(
	  -SQL => q/select count(*) from 
  (select distinct(g.gene_id) as id from gene g 
  join transcript t using (gene_id) 
  join seq_region s on (s.seq_region_id=t.seq_region_id) 
  join coord_system using (coord_system_id) 
  where species_id=? group by g.gene_id having count(*)>1) x;/,
	  -PARAMS => [$ga->species_id()]);

	if ($splice_count == 0) {
	  print "No splices detected - setting all exons as constitutive...\n";
	  $ga->dbc()->sql_helper()->execute_update(
		-SQL => q/update exon e, seq_region s, coord_system c 
	set e.is_constitutive=1 
	where e.seq_region_id=s.seq_region_id 
	and c.coord_system_id=s.coord_system_id and c.species_id=?/,
		-PARAMS => [$ga->species_id()]);
	  return;
	}
  }

  my @gene_dbIDs = sort { $a <=> $b } @{$ga->list_dbIDs()};
  my $bin = scalar(@gene_dbIDs)/1000;
  my @gene_list;
  for (my $i = 0; $i < $bin; $i++) {
	push(@gene_list, @{$ga->fetch_all_by_dbID_list([splice(@gene_dbIDs, 0, 1000)])});
  }

  while (my $gene = shift(@gene_list)) {
	my @transcripts      = @{$gene->get_all_Transcripts()};
	my $transcript_count = scalar(@transcripts);
	my %exon_count;
	my %exon_object;

	foreach my $transcript (@transcripts) {
	  my @exons = @{$transcript->get_all_Exons()};
	  while (my $exon = shift(@exons)) {
		my $key = $exon->dbID();
		++$exon_count{$key};
		$exon_object{$key} = $exon;
	  }
	}

	foreach my $exon_key (keys(%exon_count)) {
	  my $exon            = $exon_object{$exon_key};
	  my $is_constitutive = $exon->is_constitutive();
	  if ($exon_count{$exon_key} == $transcript_count) {
		if (!$exon->is_constitutive()) {
		  $is_constitutive = 1;
		}
	  } else {
		if ($exon_object{$exon_key}->is_constitutive()) {
		  $is_constitutive = 0;
		}
	  }
	  $ga->dbc()->sql_helper()->execute_update(-SQL    => 'update exon set is_constitutive = ? where exon_id =?',
											   -PARAMS => [$is_constitutive, $exon->dbID()]);
	}
  } ## end while (my $gene = shift(@gene_list...))
  return;
} ## end sub run

1;
