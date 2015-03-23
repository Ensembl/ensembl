=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package Bio::Das::ProServer::SourceAdaptor::demo;

use strict;
use vars qw(@ISA);
use Bio::Das::ProServer::SourceAdaptor;
@ISA = qw(Bio::Das::ProServer::SourceAdaptor);

my @id  = qw(CLAT_HUMAN CLAT_HUMAN CLAT_HUMAN 10 AC073366  AC073366);
my @id2  = qw(CLAT_HUMAN CLAT_HUMAN not_same 10 AC073366  AC073366);
my @names = qw(PF00755 PF00755 CLAT_HUMAN 10 AC073366  AC073366); 
my @starts = qw(30 100 50 50150000 1 );
my @ends = qw(748 550 650 50200000 900);
my @method =qw(description PFAM SCOP REPEAT BLAST ANOTH);
my @types  = qw(typ1 typ2 typ3 typ4 typ5 typ6);


sub init {
  my $self                = shift;
  $self->{'dsn'}          = "demo";
  $self->{'capabilities'} = {
			     'features' => '1.0',
 			    };
}

sub length {
  return 0;
}

  
sub build_features {
  my ($self, $opts) = @_;
  my $spid    = $opts->{'segment'};
  my $start   = $opts->{'start'};
  my $end     = $opts->{'end'};

  #  print  $self->config->{'needed_arg'}."\n";
  #
  #Create some dummy notes
  #
  my @notes=();
  for(my $i=0; $i <= $#id; $i++){
    $notes[$i] = "Demo annotation $i ??";
  }

  #create array of fearutes to return;
  my @features = ();
  for(my $i=0; $i <= $#id; $i++){

    next if(($id[$i] ne $spid) ||
	    (defined($start) and ($ends[$i] < $start or $starts[$i] > $end)));
    push @features, {
		     'id'     => $id2[$i],
		     'type'   => $types[$i],  # needed for proteinview protein
		                              # features names in graphical view
		     'feature'=> $names[$i],
		     'method' => $method[$i],
		     'start'  => $starts[$i],
		     'end'    => $ends[$i],
		     'note'  => $notes[$i],
		     'link'   => 'http://www.ensembl.org/Docs/enstour/',
		     'linktxt' => 'Tour',
		    };
    }

  return @features;
}


1;
