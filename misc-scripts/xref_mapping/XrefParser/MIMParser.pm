package XrefParser::MIMParser;

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
    print "\nUsage: MIMParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run($ARGV[0]);

}

sub run {

  my $self = shift if (defined(caller(1)));
  my $file = shift;
  my $general_source_id = shift;
  my $species_id = shift;
  my %old_to_new;
  my %removed;
  my $source_id;
  my @sources;
  if(!defined($general_source_id)){
    $general_source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }
  push @sources, $general_source_id;
 

  my %genemap;
  my %morbidmap;
  my $gene_source_id = XrefParser::BaseParser->get_source_id_for_source_name("MIM_GENE");
  push @sources, $gene_source_id;
  my $morbid_source_id =  XrefParser::BaseParser->get_source_id_for_source_name("MIM_MORBID");
  push @sources, $morbid_source_id;

  print "sources are:- ".join(", ",@sources)."\n";
  if(!open(MIM,"<MIM/genemap")){
    print  "ERROR: Could not open MIM/genemap\n";
    return 1; # 1 is an error
  }
    
  #Each entry is a list of fields, separated by the '|' character. 
  #The fields are, in order :
  
  #1  - Numbering system, in the format  Chromosome.Map_Entry_Number
  #2  - Month entered
  #3  - Day     "
  #4  - Year    "
  #5  - Location
  #6  - Gene Symbol(s)
  #7  - Gene Status (see below for codes)
  #8  - Title
  #9  - 
  #10 - MIM Number
  #11 - Method (see below for codes)
  #12 - Comments
  #13 -
  #14 - Disorders
  #15 - Disorders, cont.
  #16 - Disorders, cont
  #17 - Mouse correlate
  #18 - Reference  }

  while(<MIM>){
    my @array = split (/\|/, $_);
    $genemap{$array[9]} = 1; # array starts at 0 so 10 - 1.
  }
  close MIM;
  
  if(!open(MIM,"<MIM/morbidmap")){
    print  "ERROR: Could not open MIM/morbidmap\n";
    return 1; # 1 is an error
  }

  while(<MIM>){
    my @array = split (/\|/, $_);
    # store the description
    $morbidmap{$array[2]} = $array[0];
  }
  close MIM;



  local $/ = "*RECORD*";

  if(!open(MIM,"<".$file)){
    print  "ERROR: Could not open $file\n";
    return 1; # 1 is an error
  }
  
  my $count = 0;
  my $removed_count =0;
  <MIM>; # first record is empty with *RECORD* as the record seperator
  while (<MIM>) {
    #get the MIM number
    my $number = 0;
    my $description = undef;
    my $is_morbid = 0;
    if(/\*FIELD\*\s+NO\n(\d+)/){
      $number = $1;
      $source_id = $gene_source_id;
      if(defined($morbidmap{$number})){
#	$source_id = $morbid_source_id;
	$is_morbid=1;
      }
      elsif(defined($genemap{$number})){
#	$source_id = $gene_source_id;
      }
      else{
	if(/\*FIELD\*\sTI\n([\^\#\%\+\*]*)\d+(.*)\n/){
	  if($1 eq "^"){
	    if(/\*FIELD\*\sTI\n[\^]\d+ MOVED TO (\d+)/){
	      $old_to_new{$number} = $1;
	    }
	    else{
	      $removed{$number} = 1;
	      $removed_count++;
	    }
	    next;
	  }
	}
	next;
#	$source_id = $general_source_id; 
      }
   }
    if($number==0){
      die "ERROR $_";
    }
    if(/\*FIELD\*\sTI\n([\^\#\%\+\*]*)\d+(.*)\n/){
      if($1 eq "^"){
	if(/\*FIELD\*\sTI\n[\^]\d+ MOVED TO (\d+)/){
	  $old_to_new{$number} = $1;
	}
	else{
	  $removed{$number} = 1;
	  $removed_count++;
	}
	next;
      }
      if(!defined($2) or $2 eq ""){
	die "No descripton for $number\n";
      }
      else{
	$description =$2;
	$description =~ s/\;\s[A-Z0-9]+$//; # strip gene name at end
	$count++;
	$self->add_xref($number,"",$number,$description,$gene_source_id,$species_id);
# now also add morbid entry
	if($is_morbid){
	  $self->add_xref($number,"",$number,$morbidmap{$number},$morbid_source_id,$species_id);
	}			}
      #	print $number."\n*".$description."*\n" unless ($count > 100); 
    }
  }
  my $syn_count =0;
  foreach my $mim (keys %old_to_new){
    my $old= $mim;
    my $new= $old_to_new{$old};
    while(defined($old_to_new{$new})){
      $new = $old_to_new{$new};
    }
    if(!defined($removed{$new})){
      $self->add_to_syn_for_mult_sources($new, \@sources, $old);
      $syn_count++;
    }
  }
  print "$count MIM xrefs added\n";
  print "added $syn_count synonyms (defined by MOVED TO)\n";
  return 0; #successful
}



sub new {

  my $self = {};
  bless $self, "XrefParser::MIMParser";
  return $self;

}
 
1;
