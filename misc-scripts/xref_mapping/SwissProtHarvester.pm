# Harvester for downloading SwissProt files

package SwissProtHarvester;

use strict;

use vars qw(@ISA);
@ISA = qw(BaseHarvester);

# --------------------------------------------------------------------------------

run() if (!defined(caller()));

# --------------------------------------------------------------------------------

sub run {

  # URL should end with a /
  my $url = "ftp://ftp.ebi.ac.uk/pub/databases/SPproteomes/swissprot_files/proteomes/";
  my @taxons = BaseHarvester->taxonomy_ids();
  my $ext = ".SPC";

  foreach my $taxon (@taxons) {

    my $file = $url . $taxon . $ext;
    my $result = system("wget", "--quiet", "--timestamping", $file);

  }

}

# TODO logging, error handling

# --------------------------------------------------------------------------------

sub new {

  my $self = {};
  bless $self, "SwissProtHarvester";
  return $self;

}

# --------------------------------------------------------------------------------

1;

