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
  my $source_id = shift;
  my $species_id = shift;
  my %old_to_new;
  my %removed;

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


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
    if(/\*FIELD\*\s+NO\n(\d+)/){
      $number = $1;
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
	$self->add_xref($number,"",$number,$description,$source_id,$species_id);
	#	print $number."\n*".$description."*\n" unless ($count > 100); 
      }
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
      $self->add_to_syn($new,$source_id,$old);
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
