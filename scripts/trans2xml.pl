#!/usr/local/bin/perl

=head1 NAME

trans2xml

=head1 SYNOPSIS
 
 trans2xml.pl

=head1 DESCRIPTION

This script produces GAME XML file outputs of the whole database, where each contig produces a new
annotation ID and each transcript being a new feature_set annotation.

=head1 OPTIONS

    -dbtype    Database tpye (only used for TimDB)

    -host      host name for database (gets put as host= in locator)

    -port      For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (dbuser= in locator)

    -dbpass    For RDBs, what password to use (dbpass= in locator)

    -module    Module name to load (Defaults to Bio::EnsEMBL::DBSQL::Obj)

    -help      Displays script documentation with PERLDOC

=cut

use strict;
use Bio::EnsEMBL::DBLoader;
#use Bio::EnsEMBL::TimDB::Obj;
use Bio::SeqIO;
use Getopt::Long;

my $host   = 'localhost';
my $port   = '410000';
my $dbname = 'ensembl07';
my $dbuser = 'root';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::Obj';
my $getall;
my $usefile;
my $help;

&GetOptions( 	     
	     'host:s'   => \$host,
	     'port:n'   => \$port,
	     'dbname:s' => \$dbname,
	     'dbuser:s' => \$dbuser,
	     'dbpass:s' => \$dbpass,
	     'module:s' => \$module,
	     'h|help'   => \$help,
	     'getall'   => \$getall,
	     'usefile'=> \$usefile
	     );

if ($help) {
    exec('perldoc', $0);
}

my $locator = "$module/host=$host;port=$port;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);
my $seqio;
print "<game version=\"1.001\">\n";
my @clones;

if ($help){
    exec('perldoc', $0);
}

if( $usefile == 1 ) {
    while( <> ) {
	my ($en) = split;
	push(@clones,$en);
    }
}
elsif ($getall == 1) {
    @clones=$db->get_all_Clone_id();
}
else {
    @clones = @ARGV;
}


foreach my $clone_id (@clones) {
print STDERR "\nDumping  clone  $clone_id\n";
    
#Get clone, and print header for the whole clone
my $clone = $db->get_Clone($clone_id);
foreach my $contig ( $clone->get_all_Contigs ) {

    #Get supporting features
    my @features;
    push(@features,$contig->get_all_SimilarityFeatures);
    
#Print contig header
    my $length= $contig->length;
    print "   <seq length=\"$length\" type=\"DNA\" id=\"".$contig->id."\">\n";
    print "      <description>\n";
    print "         Contig: ".$contig->id." from Sequence entry: ".$clone->id."\n";
    print "      </description>\n";
    print "      <name>".$contig->id."</name>\n";
    
    #Print contig dna in $nuc blocks
    print "      <residues>\n";
    my $str=$contig->seq;
    my $nuc = 60;
    my $whole = int($length / $nuc) * $nuc;
    my( $i );
    for ($i = 0; $i <= $whole; $i += $nuc) {
	my $blocks = substr($str, $i, $nuc);
	print "         $blocks\n";
    }
    print "      </residues>\n";
    print "   </seq>\n"; 
    
#Print each gene, as an <annotation>
    foreach my $gene ($contig->get_all_Genes('evidence')) {
	print STDERR "         gene   ",$gene->id,"\n";
	
	print "   <annotation>\n";
	print "      <type>gene</type>\n";
	print "      <name>".$gene->id."</name>\n";
	print "      <creation_date>".localtime($gene->created)."</creation_date>\n";
	print "      <modification_date>".localtime($gene->modified)."<modification_date>\n";
	print "      <version>".$gene->version."</version>\n";
	
#Print each transcript as a <feature set>
	foreach my $trans ( $gene->each_Transcript ) {
	    print STDERR"         trans. ",$trans->id,"\n";
	    
	    print "      <feature_set>\n";
	    print "         <type>transcript</type>\n";
	    print "         <name>".$trans->id."</name>\n";
	    print "         <version>".$trans->version."</version>\n";
	    print "         <seq_relationship seq=\"".$contig->id."\">\n";
	    print "            <span>\n";
	    my $switch = 1;
	    my $trans_start;
	    my $trans_end;
	    foreach my $trans_exon ($trans->each_Exon) {
		if ($trans_exon->contig_id eq $contig->id) {
		    if ($switch == 1) {
			$trans_start = $trans_exon->start;
			$switch = 0;
	    }
		    $trans_end = $trans_exon->end;
		}
	    }
	    print "               <start>$trans_start</start>\n";
	    print "               <end>$trans_end</end>\n";
	    print "            </span>\n";
	    print "         </seq_relationship>\n";
	    
	    foreach my $exon ( $trans->each_Exon ) {
		if ($exon->contig_id eq $contig->id) {
		    print "         <feature_span id=\"".$exon->id."\">\n";
		    print "            <type>Exon</type>\n";
		    print "            <name>".$exon->id."</name>\n";
		    print "            <version>".$exon->version."</version>\n";
		    print "            <creation_date>".localtime($exon->created)."</creation_date>\n";
		    print "            <modification_date>".localtime($exon->modified)."</modification_date>\n";
		    print "            <seq_relationship seq=\"",$exon->contig_id,"\">\n";
		    print "               <span>\n";
		    my $strand = $exon->strand;
		    if ($strand == -1) {
			print "                  <start>",$exon->end,"</start>\n";
			print "                  <end>",$exon->start,"</end>\n";
		    }
		    else {
			print "                  <start>",$exon->start,"</start>\n";
			print "                  <end>",$exon->end,"</end>\n";
		    } 
		    print "               </span>\n";
		    print "            </seq_relationship>\n";
		    
		    foreach my $feature ($exon->each_Supporting_Feature) {
			my $source=$feature->source_tag;
			my $type="unknown";
			my $name = $feature->hseqname;
			my $db = "unkown";
			my $acc = "unknown";
			
			if ($source eq "blastp") {
			    $type="homology";
			    $name =~ /(\w+):(\w+)/;
			    if ($1 eq "SW") { $db="SWISSPROT";}
			    elsif ($1 eq "TR") {$db="TREMBL";}
			    elsif ($1 eq "WP") {$db="WORM_PEP";}
			    
			    $acc = $2;
			    if (!$1 || !$2) {
				print STDERR "hseqname $name not understood\n";
			    }
			}
			
			elsif ($source eq "tblastn") {
			    $type="homology";
			    if ($name =~ /gb\|(\w+)/) {
				$db = "GENBANK";
				$acc = $1;
			    }
			    elsif ($name =~ /dbj\|(\w+)/) {
				$db = "DDBJ";
				$acc = $1;
			    }
			    elsif ($name =~ /emb\|(\w+)/) {
				$db = "EMBL";
				$acc = $1;
    }
			    elsif ($name =~ /\|/) {
				print STDERR "hseqname $name not understood\n";
			    }
			    else {
				$db="EMBL";
				$acc=$name;
			    }
			}
			
			elsif ($source eq "hmmpfam") {
			    $type="homology";
			    $db="hmmpfam";
			    $acc=$name;
			}
			else {
			    print STDERR "Source $source not understood!\n";
			}
			
			print "            <evidence type=\"$type\">\n";
			print "               <result_span>\n";
			print "                  <program>".$feature->source_tag."</program>\n";
			print "                  <dbxref>\n";
			print "                     <db>$db</db>\n";
			print "                     <xref_db_id>$acc</xref_db_id>\n";
			print "                  </dbxref>\n";
			print "                  <score>".$feature->score."</score>\n";
			print "                  <seq_relationship type=\"query\" seq=\"".$contig->id."\">\n";
my $strand= $feature->strand;
			if ($strand == -1) {
			    print "                     <start>".$feature->end."</start>\n";
			    print "                     <end>".$feature->start."</end>\n";
			}
			else {
			    print "                     <start>".$feature->start."</start>\n";
			    print "                     <end>".$feature->end."</end>\n";
			}
			print "                  </seq_relationship>\n";
			print "                  <seq_relationship type=\"subject\" seq=\"".$feature->hseqname."\">\n";
			my $strand=$feature->hstrand;
			
			if ($strand == -1) {
			    print "                     <start>".$feature->hend."</start>\n";
			    print "                     <end>".$feature->hstart."</end>\n";
			}
			else {
			    print "                     <start>".$feature->hstart."</start>\n";
			    print "                     <end>".$feature->hend."</end>\n";
			}
			print "                  </seq_relationship>\n";
			print "               </result_span>\n";
			print "            </evidence>\n";
		    }
print "         </feature_span>\n";
		}
	    }
print "      </feature_set>\n";
	    print "      <feature_set id=\"".$trans->translation->id."\">\n";
	    print "         <type>translation</type>\n";
	    print "         <name>".$trans->translation->id."</name>\n";
	    print "         <version>".$trans->translation->version."</version>\n";
	    print "         <seq_relationship seq=\"".$contig->id."\">\n";
	    print "            <span>\n";
	    $switch = 1;
	    my $tr_start;
	    my $tr_end;
	    foreach my $tr_exon ($trans->translateable_exons) {
		if ($tr_exon->contig_id eq $contig->id) {
		    if ($switch == 1) {
			$tr_start = $tr_exon->start;
			$switch = 0;
		    }
		    $tr_end = $tr_exon->end;
		}
	    }
	    print "               <start>$tr_start</start>\n";
	    print "               <end>$tr_end</end>\n";
	    print "            </span>\n";
	    print "         </seq_relationship>\n";
	    print "         <seq length=\"$length\" type=\"AA\" id=\"".$trans->translation->id."\">\n";
	    print "            <description>\n";
	    print "               Translation: ".$trans->translation->id." from Contig: ".$contig->id." from Sequence Entry: ".$clone->id."\n";
	    print "            </description>\n";
	    print "            <name>".$trans->translation->id."</name>\n";
	    
#Print contig dna in $nuc blocks
print "            <residues>\n";
	    my $str=$trans->translate->seq;
	    my $length=$trans->translate->length;
	    my $nuc = 60;
	    my $whole = int($length / $nuc) * $nuc;
	    my( $i );
	    for ($i = 0; $i <= $whole; $i += $nuc) {
		my $blocks = substr($str, $i, $nuc);
		print "               $blocks\n";
	    }
	    print "            </residues>\n";
	    print "         </seq>\n";
	    print "      </feature_set>\n";
	}
	print "   </annotation>\n";
    }
}

#print "   <computational_analysis>\n";
#print "      <date>".localtime($clone->created)."</date>\n";
#print "      <program>assembly</program>\n";
#print "      <version>".$clone->version."</version>\n";
#print "      <result_set seq=\"".$clone->id."\">\n";

#foreach my $contig ($clone->get_all_Contigs) {
#    print "         <result_span id =\"".$contig->id."\">\n";
#    print "            <name>".$contig->id."</name>\n";
#    print "            <seq_relationship seq=\"".$clone->id."\">\n";
#    print "               <span>\n";
#    print "                  <start>".$contig->offset."</start>\n";
#    my $end=$contig->offset+$contig->length;
#    print "                  <end>$end</end>\n";
#    print "               </span>\n";
#    print "            </seq_relationship>\n";
#    print "         </result_span>\n";
#}
#print "      </result_set>\n";
#print "   </computational_analaysis>\n";
}
print "</game>\n";

