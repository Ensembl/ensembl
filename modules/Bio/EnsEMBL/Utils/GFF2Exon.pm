# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

GFF2Exon

=head1 SYNOPSIS

use Bio::EnsEMBL::Utils::GFF2Exon;
my $gff  = GFF2Exon->new(-gff_version => 2);
my $exon = $gff->exon_from_gff( $gff_string );

=cut

# Let the code begin...

package Bio::EnsEMBL::Utils::GFF2Exon;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Root;
@ISA = qw(Bio::EnsEMBL::Root);


sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($gff_version) = $self->_rearrange([qw(GFF_VERSION)],@args);

  if(($gff_version != 1) && ($gff_version != 2)) {
      $self->throw("unknown version: $gff_version");
  }
  $self->gff_version($gff_version);
  return $self;
}

=head2 exon_from_gff

 Usage   : my $exon = $gff->exon_from_gff_string($gff_string);
 Function: Sets properties of an exon object from a GFF-formatted
           string. Interpretation of the string depends on the version
           that has been specified at initialization.
 Args    : The GFF-formatted string

=cut

sub exon_from_gff {
  my ($self, $gff_string, $trans_hash ) = @_;
  
  if($self->gff_version() == 1)  {
    return $self->_from_gff1_string($gff_string, $trans_hash);
  } 
  else{
    return $self->_from_gff2_string($gff_string, $trans_hash);
  }
}

=head2 _from_gff1_string

    Returns : Exon
    Args    : The GFF-formatted string to initialize it from

=cut

sub _from_gff1_string {
  my ($self, $gffstring, $trans_hash) = @_;
    
  #print STDERR "gff_string: $gffstring\n";

  #39896	human_cdna	exon	1030942	1031591	100	+	0	6604
  #39897	human_cdna	exon	1034227	1034285	100	+	0	6604
  #39898	human_cdna	exon	1035280	1035339	100	+	0	6604
  #39899	human_cdna	exon	1035741	1035820	100	+	0	6604
  #39900	human_cdna	exon	1036477	1036629	100	+	0	6604
  #243299	human_cdna	exon	1037161	1037229	100	+	0	40150
  #243301	human_cdna	exon	1043180	1043315	100	+	0	40150
  #243303	human_cdna	exon	1048194	1048340	100	+	0	40150
  #243305	human_cdna	exon	1053724	1053855	100	+	0	40150
  #243307	human_cdna	exon	1054346	1054451	100	+	0	40150
  #243309	human_cdna	exon	1055667	1055833	100	+	0	40150
  #243310	human_cdna	exon	1057048	1057221	100	+	0	40150
  #243311	human_cdna	exon	1058529	1058630	100	+	0	40150
  #243313	human_cdna	exon	1061487	1062080	100	+	0	40150
  my ($seqname, 
      $source, 
      $primary, 
      $start, 
      $end, 
      $score, 
      $strand, 
      $frame, 
      @group) 
    = split(/\s+/, $gffstring);
  
  if ( !defined $frame ) {
    $self->throw("[$gffstring] does not look like GFF to me");
  }
  $frame = 0 unless( $frame =~ /^\d+$/);
    
  my $exon = Bio::EnsEMBL::Exon->new();
  $exon->seqname($seqname);
  #$exon->source_tag($source);
  $exon->start($start);
  $exon->end($end);
  my $phase = ( 3 - $frame )%3;
  $exon->phase($phase);
  $exon->end_phase( ( $exon->phase + $exon->length)%3 );
  if ( $score ){
    $exon->score( $score );
  }
  if ( $strand eq '-' ) { $exon->strand(-1); }
  if ( $strand eq '+' ) { $exon->strand(1); }
  if ( $strand eq '.' ) { $exon->strand(0); }

  ############################################################
  # warning: it parses only the first element of the group
  
  $trans_hash->{ $exon } = $group[0];
  return $exon;
}

=head2 _from_gff2_string
    
    Returns : Exon
    Args    : The GFF2-formatted string to initialize it from

=cut

sub _from_gff2_string {
   my ($self, $gffstring) = @_;
   chomp($gffstring);

   # according to the Sanger website, GFF2 should be single-tab separated elements, and the
   # free-text at the end should contain text-translated tab symbols but no "real" tabs,
   # so splitting on \t is safe, and $attribs gets the 
   # entire attributes field to be parsed later

   my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, @attribs) 
       = split(/\t+/, $gffstring);

   # just in case the rule against tab characters has been broken
   my $attribs = join '', @attribs;  
   
   if ( !defined $frame ) {
       $self->throw("[$gffstring] does not look like GFF2 to me");
   }
   my $exon = Bio::EnsEMBL::Exon->new();
   $exon->seqname($seqname);
   $exon->source_tag($source);
   $exon->analysis->gff_feature($primary);
   $exon->start($start);
   $exon->end($end);
   $exon->phase($frame);
   $exon->end_phase( ( $exon->phase + $exon->length)%3 );
   if ( $score eq '.' ) {
       #$exon->score(undef);
   } else {
       $exon->score($score);
   }
   if ( $strand eq '-' ) { $exon->strand(-1); }
   if ( $strand eq '+' ) { $exon->strand(1); }
   if ( $strand eq '.' ) { $exon->strand(0); }

   #  <Begin Inefficient Code from Mark Wilkinson>
   # this routine is necessay to allow the presence of semicolons in quoted text
   # Semicolons are the delimiting character for new tag/value attributes.
   # it is more or less a "state" machine, with the "quoted" flag going up and down
   # as we pass thorugh quotes to distinguish free-text semicolon and hash 
   # symbols from GFF control characters
   
   # remove comments field (format:  #blah blah blah...  at the end of the GFF line)
   #$attribs =~ s/\#(.*)$//;				 
   #my @att = split //, $attribs;    # split into individual characters
   #my $num = $#att;                 # count them
   #my $flag = 0;
   #my @parsed;		# this is needed to hold the characters that have been parsed	
   
   # run through each character one at a time and check it
   #for (my $a = 0; $a <= $num ; $a +=1){   
   #    if ($att[$a] eq "\""){$flag=($flag==0)?1:0}  # flag up on entering quoted text, down on leaving it
   #    if (($att[$a] eq ";") && $flag){$att[$a] = "INSERT_SEMICOLON_HERE"}  # replace semicolon with an unusual message if the quoted text flag is up
   #    if (($att[$a] eq "#") && !$flag){last}  # an unquoted hash symbol means the beginning of the comments field - discard
   #    push @parsed, $att[$a]                  # take the parsed character and push it onto the parsed list
   #    }
   
   #$attribs = join "", @parsed; # rejoin into a single string
   
   # <End Inefficient Code>   Please feel free to fix this and make it more "perlish"
   
  # my @key_vals = split /;/, $attribs;   # attributes are semicolon-delimited
   
  # foreach my $pair ( @key_vals ) {
   #    $pair =~ s/INSERT_SEMICOLON_HERE/;/g;  # replace semicolons that were removed from free-text above.
   #    my ($blank, $key, $values) = split  /^\s*([\w\d]+)\s/, $pair;	# separate the key from the value
       
    #   my @values;								
       
     #  while ($values =~ s/"(.*?)"//){          # free text is quoted, so match each free-text block
	#   if ($1){push @values, $1};          # and push it on to the list of values (tags may have more than one value...)
       #}
       
      # my @othervals = split /\s+/, $values;  # and what is left over should be space-separated non-free-text values
      # foreach my $othervalue(@othervals){
#	   if (CORE::length($othervalue) > 0){push @values, $othervalue}  # get rid of any empty strings which might result from the split
#       }
       
#       foreach my $value(@values){
#	   $exon->add_tag_value($key, $value);
#       }
 #  }

   return $exon;
}

=head2 gff_string

=cut
    
sub gff_string{
    my ($self, $exon,$transcript) = @_;
    
    if($self->gff_version() == 1) {
	return $self->_gff1_string($exon,$transcript);
    } else {
	return $self->_gff2_string($exon,$transcript);
    }
}

=head2 _gff1_string

=cut

sub _gff1_string{
  my ($self, $exon, $transcript) = @_;
  my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
  
  if( $exon->can('score') ) {
    $score = $exon->score();
  }
  $score = '.' unless defined $score;
  
  if( $exon->can('frame') ) {
    $frame = $exon->frame();
  }
  $frame = '.' unless defined $frame;
  
  $strand = $exon->strand();
  if(! $strand) {
    $strand = ".";
  } elsif( $strand == 1 ) {
    $strand = '+';
  } elsif ( $exon->strand == -1 ) {
    $strand = '-';
  }
  
  $name        = $exon->seqname();
  $source      = "merged";
  $primary_tag = "exon";
  
  $str = join("\t",
	      $name,
	      $source,
	      $primary_tag,
	      $exon->start(),
	      $exon->end(),
	      $score,
	      $strand,
	      $frame);
  
  my $tag_str = $transcript->type;
  $str .= "\t$tag_str";
  return $str;
}

=head2 _gff2_string

=cut

sub _gff2_string{
   my ($gff, $exon) = @_;
   my ($str,$score,$frame,$name,$strand);

   if( $exon->can('score') ) {
       $score = $exon->score();
   }
   $score = '.' unless defined $score;

   if( $exon->can('frame') ) {
       $frame = $exon->frame();
   }
   $frame = '.' unless defined $frame;

   $strand = $exon->strand();
   if(! $strand) {
       $strand = ".";
   } elsif( $strand == 1 ) {
       $strand = '+';
   } elsif ( $exon->strand == -1 ) {
       $strand = '-';
   }

   if( $exon->can('seqname') ) {
       $name = $exon->seqname();
       $name ||= 'SEQ';
   } else {
       $name = 'SEQ';
   }


   $str = join("\t",
                 $name,
		 $exon->source_tag(),
		 $exon->analysis->gff_feature(),
		 $exon->start(),
		 $exon->end(),
		 $score,
		 $strand,
		 $frame);

   # the routine below is the only modification I made to the original
   # ->gff_string routine (above) as on November 17th, 2000, the
   # Sanger webpage describing GFF2 format reads: "From version 2
   # onwards, the attribute field must have a tag value structure
   # following the syntax used within objects in a .ace file,
   # flattened onto one line by semicolon separators. Tags must be
   # standard identifiers ([A-Za-z][A-Za-z0-9_]*).  Free text values
   # must be quoted with double quotes".

   # MW

   #my $valuestr;
   #if ($exon->all_tags){  # only play this game if it is worth playing...
   #     $str .= "\t";     # my interpretation of the GFF2 specification suggests the need for this additional TAB character...??
   #     foreach my $tag ( $exon->all_tags ) {
   #         my $valuestr; # a string which will hold one or more values for this tag, with quoted free text and space-separated individual values.
   #         foreach my $value ( $exon->each_tag_value($tag) ) {
   #      		if ($value =~ /[^A-Za-z0-9_]/){
   #      			$value =~ s/\t/\\t/g;          # substitute tab and newline characters
   #      			$value =~ s/\n/\\n/g;          # to their UNIX equivalents
   #      			$value = '"' . $value . '" '}  # if the value contains anything other than valid tag/value characters, then quote it
   #      		$valuestr .=  $value . " ";								# with a trailing space in case there are multiple values
   #      															# for this tag (allowed in GFF2 and .ace format)		
   #         }
   #         $str .= "$tag $valuestr ; ";                              # semicolon delimited with no '=' sign
   #     }
   #		chop $str; chop $str  # remove the trailing semicolon and space
   # }
   return $str;
}

=head2 gff_version

=cut

sub gff_version {
  my ($self, $value) = @_;
  if(defined $value && (($value == 1) || ($value == 2))) {
      $self->{'GFF_VERSION'} = $value;
  }
  return $self->{'GFF_VERSION'};
}


1;

