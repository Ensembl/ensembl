#!/ebi/extserv/bin/perl/bin/perl
# an example script demonstrating the use of BioMart webservice


#
# NOTE this could have implemented in the parser itself but the data is needed
# for the simple features so
#
use strict;
use LWP::UserAgent;


 my $xml = (<<XXML);
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
  
  <Dataset name = "dcc" interface = "default" >
    <Attribute name = "mgi_accession_id" />
    <Attribute name = "marker_symbol" />
    <Attribute name = "vector_available" />
    <Attribute name = "escell_available" />
    <Attribute name = "mouse_available" />
    <Attribute name = "ensembl_gene_id" />
    </Dataset>
</Query>
XXML

open (OUT,">ensembl_ikmc_initial.txt");

my $path="http://www.i-dcc.org/biomart/martservice?";
my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
my $ua = LWP::UserAgent->new;

my $response;

$ua->request($request,
	     sub{
		 my($data, $response) = @_;
		 if ($response->is_success) {
		     print OUT "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);

close OUT;

my %symbols;
my %ensembl_ids;
my %status;

open (IN,"ensembl_ikmc_initial.txt");
#nb [9] is now cell_line_bg and [10] is backcross
while (<IN>){
    my @line = split(/\t/,$_);
    chop $line[5];
    my $mgi_id = $line[0];
    $symbols{$mgi_id}=$line[1];
    $ensembl_ids{$mgi_id}=$line[5];
    $status{$mgi_id} = 1 if ($status{$mgi_id} eq '');

    if ($status{$mgi_id} < 4 && $line[4] == 1){
	$status{$mgi_id} = 4;
    }
    elsif ($status{$mgi_id} < 3 && $line[3] == 1){
	$status{$mgi_id} = 3;
    }
    elsif ($status{$mgi_id} < 2 && $line[2] == 1){
	$status{$mgi_id} = 2;
    }

}
close IN;

open (OUT,">ensembl_ikmc_xref.txt");
foreach my $mgi_id(keys %symbols){
    my $description;
    $description = 'No products available yet' if $status{$mgi_id} == 1;
    $description = 'Vector available' if $status{$mgi_id} == 2;
    $description = 'ES cells available' if $status{$mgi_id} == 3;
    $description = 'Mice available' if $status{$mgi_id} == 4;

    print OUT "$mgi_id\t$symbols{$mgi_id}\t$description\t$ensembl_ids{$mgi_id}\n";
}
close OUT;
