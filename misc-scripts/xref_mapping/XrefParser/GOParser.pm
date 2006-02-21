

package XrefParser::GOParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XrefParser::BaseParser;

use vars qw(@ISA);
@ISA = qw(XrefParser::BaseParser);



# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: GoParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $source_id = shift;
  my $species_id = shift;
  my %wrongtype;

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


  my (%swiss) = %{XrefParser::BaseParser->get_valid_codes("uniprot",$species_id)};
  my (%refseq) = %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};

  # complication with GO xrefs from JAX - linked to MGI symbols, which are themselves
  # dependent, so we need to get the MGI->Uniprot mapping and store the *Uniprot*
  # as the master xref
  my (%mgi_to_uniprot) = %{XrefParser::BaseParser->get_existing_mappings("MarkerSymbol", "Uniprot/Swissprot", $species_id)};

  my %worm;
  my %worm_label;
  my $wormset;
  my %fish;
  my $fishset;

  my $count  = 0;

  if(!open(GO,"<".$file)){
    print "ERROR: Could not open $file\n";
    return 1; # 1 error
  }
  my $taxon_line = "taxon:".$species_id;
  my $miss =0;
  while (<GO>) {
    if(/$taxon_line/){
      chomp;
      my @array = split (/\t/,$_);
      $array[9] =~ s/\'/\\\'/g;
      my $master=0;
      if($array[0] =~ /ENSEMBL/){
        #these might be good for a check
        # match GO to Uniprot
        # match Uniprot to ENSEMBL
        # check ENSEMBL's are the same.
      }
      elsif($array[0] =~ /RefSeq/){
        if($refseq{$array[1]}){
          XrefParser::BaseParser->add_to_xrefs($refseq{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
          $count++;
        }
      }
      elsif($array[0] =~ /UniProt/){
        if($swiss{$array[1]}){
          XrefParser::BaseParser->add_to_xrefs($swiss{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
          $count++;
        }
      }
      elsif($array[0] =~ /^WB/){
	#WB      CE20707 ZYG-9           GO:0008017      WB:WBPaper00003099|PMID:9606208 ISS             F                       protein  taxon:6239      20030829        WB
        if(!defined($wormset)){
          $wormset = 1;
          %worm = %{XrefParser::BaseParser->get_valid_xrefs_for_direct_xrefs('worm')};
        }
	my $worm_acc=$array[1];
        if(!defined($worm{$worm_acc})){ 
	  if(defined($worm{$array[10]})){
	    $worm_acc = $array[10];
	  }
	  elsif(defined($worm{$array[2]})){
	    $worm_acc = $array[2];
	  }
	}

        if(defined($worm{$worm_acc})){ 	
	  my ($xref_id, $stable_id, $type, $link) = split(/::/,$worm{$worm_acc});
	  
	  my $new_xref_id=$self->get_xref($array[4],$source_id);
	  
	  if(!defined($new_xref_id)){
	    $new_xref_id = $self->add_xref($array[4],undef,$array[4],"", $source_id, $species_id);
	    $count++;
	  }
	  if(!defined($self->get_direct_xref($stable_id,$type, $array[6]))){
	    $self->add_direct_xref($new_xref_id, $stable_id, $type, $array[6]);
	  }
	}
	else{
	  $miss++;
	}
      }
      elsif($array[0] =~ /^ZFIN/){
	#ZFIN    ZDB-GENE-030131-5418    rfng            GO:0030902      ZFIN:ZDB-PUB-050125-4|PMID:15659486     IMP     ZFIN:ZDB-MRPHLNO-050308-5     radical fringe homolog (Drosophila)              gene    taxon:7955      20050310        ZFIN
        if(!defined($fishset)){
          $fishset = 1;
          %fish = %{XrefParser::BaseParser->get_valid_xrefs_for_dependencies
              ('ZFIN_ID','Uniprot/SPTREMBL','RefSeq_peptide',
               'Uniprot/SWISSPROT')};
        }
        if(defined($fish{$array[1]})){
          XrefParser::BaseParser->add_to_xrefs($fish{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
          $count++;
        }
      }

      elsif($array[0] =~ /MGI/){
        if($mgi_to_uniprot{$array[1]}){
          XrefParser::BaseParser->add_to_xrefs($mgi_to_uniprot{$array[1]},$array[4],'',$array[4],'',$array[6],$source_id,$species_id);
          $count++;
        }
      }

      elsif(!defined($wrongtype{$array[0]})){
        print STDERR "WARNING: unknown type ".$array[0]."\n";
        $wrongtype{$array[0]} = 1;
      }
    }
  }
  print "\t$count GO dependent xrefs added $miss not found\n"; 
  return 0;
}

sub new {

  my $self = {};
  bless $self, "XrefParser::GOParser";
  return $self;

}
 
1;
