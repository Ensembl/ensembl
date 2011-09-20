package XrefParser::GOParser;

use strict;
use warnings;
use POSIX qw(strftime);
use File::Basename;
use Carp;
use base qw( XrefParser::BaseParser );

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $release_file = $ref_arg->{rel_file};
  my $verbose      = $ref_arg->{verbose};

  if((!defined $source_id) or (!defined $species_id) or (!defined $files) ){
    croak "Needs to pass source_id, species_id and  files as pairs";
  }
  $verbose |=0;


  my $file = @{$files}[0];
  my $file_desc = @{$files}[1];

  #
  # Get the descriptions from the desc file.
  #

  my $go_desc_io = $self->get_filehandle($file_desc);

  if ( !defined $go_desc_io ) {
    print STDERR "ERROR: Could not open description file, $file_desc\n";
    return 1;    # 1 error
  }

  my %go_to_desc;
  print "description file for GO\n" if($verbose);
  my $term = undef;
  my $desc = undef;
  while ( $_ = $go_desc_io->getline() ) {
    if(/\<id\>   # start of id tag
       (GO:\d+)  # GO: followed by the id
       \<\/id\>  # end of id tag
       /x){
      $term = $1;
    }
    elsif(/\<name\>   # start of name tag
	   (.*)       # the name we want
	   \<\/name\> # end of name tag
	   /x){
      if(defined($term)){
        $go_to_desc{$term} = $1;
      } 
      $term = undef;
    }
  }
  $go_desc_io->close();


  my %wrongtype;

  #get the "main" GO source id.
  $source_id = $self->get_source_id_for_source_name("GO","main");


  #get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);

  my $swiss_miss=0;
  my (%swiss) = %{$self->get_valid_codes("uniprot/",$species_id)};

  my $refseq_miss=0;
  my (%refseq) = %{$self->get_valid_codes("refseq",$species_id)};

  my %mouse;
  my $mouse_set;

  my %worm;
  my $wormset;

  my %fish;
  my $fishset;

  my %cerevisiae;
  my $cerevisiae_set;

  my $count  = 0;



  my %sp2tax     =  $self->species_id2taxonomy();  #some species have multiple
                                                                    #tax_id i.e. strains

  my @tax_ids = @{$sp2tax{$species_id}};
  foreach my $tax_id ( @tax_ids){

    my $go_io = $self->get_filehandle($file);

    if ( !defined $go_io ) {
      print STDERR "ERROR: Could not open $file\n";
      return 1;    # 1 error
    }

    print "processing for taxon: $tax_id\n" if($verbose);
    my $taxon_line = "taxon:".$tax_id;
    my $miss =0;
    while ( $_ = $go_io->getline() ) {
      next if(!(/$taxon_line/x));

      chomp;
      my @array = split (/\t/x,$_);
	
      # Skip "NOT" terms entirely
      next if ($array[3] eq "NOT");
	
      $array[9] =~ s/\'/\\\'/gx; # replace ' with \'
      my $master=0;
      if($array[0] =~ /ENSEMBL/){
	#these might be good for a check
	# match GO to Uniprot
	# match Uniprot to ENSEMBL
	# check ENSEMBL's are the same.
      }
      elsif($array[0] =~ /RefSeq/){
	if($refseq{$array[1]}) {
	  foreach my $xref_id (@{$refseq{$array[1]}}){
	    $self->add_dependent_xref({ master_xref_id => $xref_id,
					acc            => $array[4],
					label          => $array[4],
					desc           => $go_to_desc{$array[4]} || '',
					linkage        => $array[6],
					source_id      => $source_id,
					species_id     => $species_id} );
	    $count++;
	  }
	}
	else{
	  $refseq_miss++;
	}
      }
      elsif($array[0] =~ /UniProt/){
	if($swiss{$array[1]}){
	  foreach my $xref_id (@{$swiss{$array[1]}}){
	    $self->add_dependent_xref({ master_xref_id => $xref_id,
					acc            => $array[4],
					label          => $array[4],
					desc           => $go_to_desc{$array[4]} || '',
					linkage        => $array[6],
					source_id      => $source_id,
					species_id     => $species_id} );
	    $count++;
	  }
	}
	else{
	  $swiss_miss++;
	}
      }
      elsif($array[0] =~ /^WB/x){
	#WB      CE20707 ZYG-9           GO:0008017      WB:WBPaper00003099|PMID:9606208 ISS             F                       protein  taxon:6239      20030829        WB
	if(!defined($wormset)){
	  $wormset = 1;
	  %worm = %{$self->get_valid_xrefs_for_direct_xrefs('worm')};
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
	  my ($xref_id, $stable_id, $type, $link) = split(/::/x,$worm{$worm_acc});

	  my $new_xref_id = $self->get_xref($array[4],$source_id, $species_id);

	  if(!defined($new_xref_id)){
	    $new_xref_id = 
	      $self->add_xref($array[4],undef,$array[4],"", $source_id, $species_id, "DIRECT");
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
      elsif($array[0] =~ /^ZFIN/x){
	#ZFIN    ZDB-GENE-030131-5418    rfng            GO:0030902      ZFIN:ZDB-PUB-050125-4|PMID:15659486     IMP     ZFIN:ZDB-MRPHLNO-050308-5     radical fringe homolog (Drosophila)              gene    taxon:7955      20050310        ZFIN
	if(!defined($fishset)){
	  $fishset = 1;
	  %fish = %{$self->get_valid_xrefs_for_dependencies
		      ('ZFIN_ID','Uniprot/SPTREMBL','RefSeq_peptide',
		       'Uniprot/SWISSPROT')};
	}
	if(defined($fish{$array[1]})){
	    $self->add_dependent_xref({ master_xref_id => $fish{$array[1]},
					acc            => $array[4],
					label          => $array[4],
					desc           => $go_to_desc{$array[4]} || '',
					linkage        => $array[6],
					source_id      => $source_id,
					species_id     => $species_id} );
	  $count++;
	}
      }

      elsif($array[0] =~ /MGI/x){
	# MGI	MGI:1923501	0610007P08Rik		GO:0004386	MGI:MGI:1354194	IEA
	#  0         1                2         3             4                  5        6
	if(!defined($mouse_set)){
	  $mouse_set = 1;
	  # Todo: Make sure we get this hash populated
	  %mouse = %{$self->get_valid_codes("MGI",$species_id)};
	  
	  print "Got " . keys (%mouse) . " MGI ids\n";
	  
	}
	if ( $mouse{$array[1]} ){
	  foreach my $xref_id ( @{$mouse{$array[1]}} ) {
	    $self->add_dependent_xref({ master_xref_id => $xref_id,
					acc            => $array[4],
					label          => $array[4],
					desc           => $go_to_desc{$array[4]} || '',
					linkage        => $array[6],
					source_id      => $source_id,
					species_id     => $species_id} );
	    $count++;
	  }
	}
      }
      # SGD GO code
      elsif ($array[0] =~ /SGD/x) {
	
	if(!defined($cerevisiae_set)){
	  $cerevisiae_set = 1;
	  # Todo: Make sure we get this hash populated
	  %cerevisiae = %{$self->get_valid_codes("sgd",$species_id)};
	  
	  print "Got " . keys (%cerevisiae) . " cerevisiae ids\n";
	  
	}
	
	if($cerevisiae{$array[1]}){
	  foreach my $xref_id (@{$cerevisiae{$array[1]}}){
	    $self->add_dependent_xref({ master_xref_id => $xref_id,
					acc            => $array[4],
					label          => $array[4],
					desc           => $go_to_desc{$array[4]} || '',
					linkage        => $array[6],
					source_id      => $source_id,
					species_id     => $species_id} );
	    $count++;
	  }
	}
      }
	
      elsif(!defined($wrongtype{$array[0]})){
	print STDERR "WARNING: unknown type ".$array[0]."\n" if($verbose);
	$wrongtype{$array[0]} = 1;
      }
    }

    $go_io->close();

    print "\t$count GO dependent xrefs added $refseq_miss refseq not found and $swiss_miss Swissprot not found \n" if($verbose); 
  }
  if ( defined $release_file ) {
    # Parse and set release information from $release_file.
    my $release_io = $self->get_filehandle($release_file);

    # Slurp mode.
    local $/;
    my $release = <$release_io>;
    $release_io->close();

    $release =~ tr /\n/ /;
    $release =~
          s/.*The following table describes.*?of (GOA.*?)<ul>.*/$1/;
    $release =~ s/<[^>]+>//gx;

    print "GO release: '$release'\n" if($verbose);
    $self->set_release( $source_id, $release );
  }

  return 0;
}

1;
