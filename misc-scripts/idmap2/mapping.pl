# this script produces idmapping for genes, transcripts, translations and exons
# It needs Database orig with stable_ids set
# it needs Database target where the stable ids are to be set
# it produces log output

# it uses direct sql for speed. All queries are in SQL.pm
# it uses LSF for some jobs I guess ...

use DBI;
use SQL;


my $sdbh = DBI->connect( "DBI:mysql:host=geordy;port=3310;database=homo_sapiens_core_new_schema", "ensro" );

# my $tdbh = 

print STDERR "Start: ",scalar(localtime()),"\n";
my $exonInfo = &SQL::orig_exon_information( $sdbh );
print STDERR "Finish: ",scalar(localtime()),"\n";

print "Count: ",scalar( @$exonInfo ),"\n";

direct_mapping( $exonInfo, undef );


# direct mappings
# contig version update



# print_arrayref( $exonInfo );

sub direct_mapping {
  my ( $old, $new ) = @_;
  
  my ( $sold, $snew );

  @{$sold} = sort { $A->{'clone_id'} <=> $B->{'clone_id'} ||
		      $A->{'clone_version'} <=> $B->{'clone_version'} ||
			$A->{'contig_offset'} <=> $B->{'contig_offset'} ||
			  $A->{'exon_start'} <=> $B->{'exon_start'} } @$old;
  
  print STDERR "Sorted: ",scalar(localtime()),"\n";
}




sub print_arrayref {
  my $aref = shift;
  
  for( my $i=1; $i <= @$aref; $i++ ) {
    print "------- new record -------\n";
    for my $key ( keys %{$aref->[$i]} ) {
      print $key,"\t",$aref->[$i]{$key},"\n";
    }
  }
}
