## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------

use strict;
use Bio::EnsEMBL::GeneComparison::GeneComparisonStats;
use Bio::EnsEMBL::DBSQL::Obj;
use lib 't';
use EnsTestDB;


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..50\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/gene_comparison.dump");

# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;   
print "ok 2\n";     # 2nd test passes.

my $standard1 = $db->get_Clone('standard1');
my $standard2 = $db->get_Clone('standard2');
my $predictor1 = $db->get_Clone('predictor1');
my $predictor2 = $db->get_Clone('predictor2');
my $predictor3 = $db->get_Clone('predictor3');

die "$0\nError fetching clones : $!" unless defined 
    ($standard1 && $standard2 && $predictor1 && $predictor2 && $predictor3);



#Test 1 compares two identical clones 
my $comparer = new Bio::EnsEMBL::GeneComparison::GeneComparisonStats($standard1, $standard1);
print STDERR "Testing identical clones\n"; 
                                        
if ($comparer->getGeneSpecificity == 1) { 
    print "ok 3\n"; 
}
else {
    print "not ok 3\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getGeneSensitivity == 1) { 
    print "ok 4\n"; 
}
else {
    print "not ok 4\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getExonSpecificity == 1) { 
    print "ok 5\n"; 
}
else {
    print "not ok 5\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getExonSensitivity == 1) { 
    print "ok 6\n"; 
}
else {
    print "not ok 6\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getBaseSpecificity == 1) { 
    print "ok 7\n"; 
}
else {
    print "not ok 7\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getBaseSensitivity == 1) { 
    print "ok 8\n"; 
}
else {
    print "not ok 8\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getMissedGeneScore == 0) { 
    print "ok 9\n"; 
}
else {
    print "not ok 9\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getWrongGeneScore == 0)  { 
    print "ok 10\n"; 
}
else {
    print "not ok 10\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getJoinedGeneScore == 1) { 
    print "ok 11\n"; 
}
else {
    print "not ok 11\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getSplitGeneScore == 1)  { 
    print "ok 12\n"; 
}
else {
    print "not ok 12\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getMissedExonScore == 0) { 
    print "ok 13\n"; 
}
else {
    print "not ok 13\n";
    print STDERR "Error comparing identical clones\n";
}

if ($comparer->getWrongExonScore == 0)  { 
    print "ok 14\n"; 
}
else {
    print "not ok 14\n";
    print STDERR "Error comparing identical clones\n";
}


#Test 2 compares two clones that have no genes in common
$comparer = new Bio::EnsEMBL::GeneComparison::GeneComparisonStats($standard1, $predictor1);
print STDERR  "Testing clones with no genes in common\n";
                                        
if ($comparer->getGeneSpecificity == 0) { 
    print "ok 15\n"; 
}
else {
    print "not ok 15\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getGeneSensitivity == 0) { 
    print "ok 16\n"; 
}
else {
    print "not ok 16\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getExonSpecificity == 0) { 
    print "ok 17\n"; 
}
else {
    print "not ok 17\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getExonSensitivity == 0) { 
    print "ok 18\n"; 
}
else {
    print "not ok 18\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getBaseSpecificity == 0) { 
    print "ok 19\n"; 
}
else {
    print "not ok 19\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getBaseSensitivity == 0) { 
    print "ok 20\n"; 
}
else {
    print "not ok 20\n";
    print STDERR "Error comparing clones with no genes in common\n";
}


if ($comparer->getMissedGeneScore == 1) { 
    print "ok 21\n"; 
}
else {
    print "not ok 21\n";
    print STDERR "Error comparing clones with no genes in common\n";
}

if ($comparer->getWrongGeneScore == 1)  { 
    print "ok 22\n"; 
}
else {
    print "not ok 22\n";
    print STDERR "Error comparing clones with no genes in common\n";
}

if ($comparer->getJoinedGeneScore == 0) { 
    print "ok 23\n"; 
}
else {
    print "not ok 23\n";
    print STDERR "Error comparing clones with no genes in common\n";
}

if ($comparer->getSplitGeneScore == 0)  { 
    print "ok 24\n"; 
}
else {
    print "not ok 24\n";
    print STDERR "Error comparing clones with no genes in common\n";
}

if ($comparer->getMissedExonScore == 1) { 
    print "ok 25\n"; 
}
else {
    print "not ok 25\n";
    print STDERR "Error comparing clones with no genes in common\n";
}

if ($comparer->getWrongExonScore == 1)  { 
    print "ok 26\n"; 
}
else {
    print "not ok 26\n";
    print STDERR "Error comparing clones with no genes in common\n";
}



#Test 3 compares two clones that have 1 gene in common and 1 unique gene each
$comparer = new Bio::EnsEMBL::GeneComparison::GeneComparisonStats($standard1, $predictor2);
print STDERR "Testing clones with 1 gene in common and 1 unique gene each\n";

if ($comparer->getGeneSpecificity == 0.5) { 
    print "ok 27\n"; 
}
else {
    print "not ok 27\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getGeneSensitivity == 0.5) { 
    print "ok 28\n"; 
}
else {
    print "not ok 28\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if (($comparer->getExonSpecificity > 0.33) && ($comparer->getExonSpecificity < 0.34)) { 
    print "ok 29\n"; 
}
else {
    print "not ok 29\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getExonSensitivity == 0.4) { 
    print "ok 30\n"; 
}   
else {
    print "not ok 30\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if (($comparer->getBaseSpecificity > 0.49) && ($comparer->getBaseSpecificity < 0.5)) { 
    print "ok 31\n"; 
}
else {
    print "not ok 31\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if (($comparer->getBaseSensitivity > 0.4) && ($comparer->getBaseSensitivity < 0.41)) { 
    print "ok 32\n"; 
}
else {
    print "not ok 32\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getMissedGeneScore == 0.5) { 
    print "ok 33\n"; 
}
else {
    print "not ok 33\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getWrongGeneScore == 0.5)  { 
    print "ok 34\n"; 
}
else {
    print "not ok 34\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getJoinedGeneScore == 1) { 
    print "ok 35\n"; 
}
else {
    print "not ok 35\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getSplitGeneScore == 1)  { 
    print "ok 36\n"; 
}
else {
    print "not ok 36\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if ($comparer->getMissedExonScore == 0.6) { 
    print "ok 37\n"; 
}
else {
    print "not ok 37\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}

if (($comparer->getWrongExonScore > 0.66) && ($comparer->getWrongExonScore < 0.67)) { 
    print "ok 38\n";
}
else {
    print "not ok 38\n";
    print STDERR "Error comparing clones with 1 gene in common and 1 unique gene each\n";
}


#Test 4 compares two clones that have have one joined and one split gene
#Compare the clones 
$comparer = new Bio::EnsEMBL::GeneComparison::GeneComparisonStats($standard2, $predictor3);
print STDERR "Testing clones that have 1 joined and 1 split gene\n";
                                           
if ($comparer->getGeneSpecificity == 0.75) { 
    print "ok 39\n"; 
}
else {
    print "not ok 39\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getGeneSensitivity == 0.75) { 
    print "ok 40\n"; 
}
else {
    print "not ok 40\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getExonSpecificity == 1)    { 
    print "ok 41\n"; 
}
else {
    print "not ok 41\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getExonSensitivity == 1)    { 
    print "ok 42\n"; 
}
else {
    print "not ok 42\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getBaseSpecificity == 1)    { 
    print "ok 43\n"; 
}
else {
    print "not ok 43\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getBaseSensitivity == 1)    { 
    print "ok 44\n"; 
}
else {
    print "not ok 44\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getMissedGeneScore == 0)    { 
    print "ok 45\n"; 
}
else {
    print "not ok 45\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getWrongGeneScore == 0)     { 
    print "ok 46\n"; 
}
else {
    print "not ok 46\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getJoinedGeneScore == 1.2)  { 
    print "ok 47\n"; 
}
else {
    print "not ok 47\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getSplitGeneScore == 1.2)   { 
    print "ok 48\n"; 
}
else {
    print "not ok 48\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getMissedExonScore == 0)    { 
    print "ok 49\n"; 
}
else {
    print "not ok 49\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

if ($comparer->getWrongExonScore == 0)     { 
    print "ok 50\n"; 
}
else {
    print "not ok 50\n";
    print STDERR "Error comparing clones that have 1 joined and 1 split gene\n";
}

$db = undef;

