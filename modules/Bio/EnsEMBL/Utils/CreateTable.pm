package CreateTable;

sub new {
  my $self = bless {}; 
  return $self;
}

sub get_table_name {
  my $this = shift;
  return $this->{'table_name'};
}

sub set_table_name {
  my $this = shift;
  $this->{'table_name'} = shift;
}


sub get_create_statement {
  my $this = shift;
  return $this->{'create_statement'}; 
}

# returns array by reference
sub get_fields{
  my $this = shift;
  return $this->{'fields'};
}

# pass array by reference
sub set_fields{
  my $this = shift;
  #my $tmp_fieds = shift;
  $this->{'fields'} = shift;
}


sub to_string {
  my $this = shift;
  return "CreateTable[" . $this->get_table_name() . ", " . $this->get_create_statement() . ", (" . scalar(@{$this->get_fields()}) . ")" . join(",", @{$this->get_fields()}) .  "]";
}


sub equals {
  my $this=shift;
  my $other= shift;

  # check names same
  if ( $this->get_table_name() ne $other->get_table_name() ) {return 0;}

  # extract fields in order to compare them
  my @this_fields = @{$this->get_fields()};
  my @other_fields = @{$other->get_fields()};
  #print "THIS " . $this->get_table_name() . "\t" . scalar(@this_fields) . " --> " . join(":", @this_fields) . "\n";
  #print "OTHER " . $other->get_table_name() . "\t" . scalar(@other_fields) . " --> " . join(":",  @other_fields)."\n";

  my $this_fields_len = scalar(@this_fields);
  my $other_fields_len = scalar(@other_fields);
  #check same number of fields
  if (  $this_fields_len != $other_fields_len ) {return 0;}

  #check same field names
  for (my $i=0; $i<$this_fields_len; ++$i) {
    if (@this_fields[$i] ne @other_fields[$i]) { return 0; }
  }

  return 1;
}



# Module function which returns 0 if no entry corresponding to
# $table_name is found, otherwise returns a corresponding CreateTable
# instance.

sub create_from_string {
  my $table_name;
  my $string;

  my @fields;
  my $create_statement;
  
  my @complete_fields;
  my $field_index;

  #shift;
  ($table_name, $string) = @_;

  # extract the sql create table statement
  $string =~ /(create\s+table\s+$table_name\s*\(([^;]*\))\s*\)\s*;)/sig;
  $create_statement = $1;
  my $tmp_fields = $2;
  # remove all (...) sections because these can contain commas which would break our mechanism for sliptting fields at commas
  $tmp_fields =~ s/\([^\)]*\)//g;
  @complete_fields = split(",", $tmp_fields);
  @fields = [];
  $field_index = 0;
  for (my $i=0; $i<scalar(@complete_fields); ++$i) {
    @complete_fields[$i] =~ /^\s*(\w+)\s/;
    my $row = $1;
    unless ($row =~/PRIMARY/i  || $row=~/KEY/i || $row=~/UNIQUE/i || $row=~/^\s*$/) {@fields[$field_index++] = $row;}
  }

  if ( !$create_statement || $create_statement eq "" ) {return 0;}
  else {
    my $v = CreateTable::new();
    $v->{'create_statement'} = $create_statement;
    $v->{'table_name'} = $table_name;
    $v->{'fields'} = \@fields;
    return $v;
  }
}





# returns array of all tables specified in schema string
sub create_all_from_string {
  (my $schema) = @_;

  my @create_tables;
  my $i = 0;
  while($schema =~ m/(create\s*table\s*(\w+)\s*\(([^;]*\))\s*\)\s*;)/sig) {
    my $create_table = CreateTable::create_from_string($2, $1);
    @create_tables[$i++] = $create_table;
  }
  return @create_tables; # return by copying
}

1;
