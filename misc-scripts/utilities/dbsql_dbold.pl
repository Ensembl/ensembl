#USe with:
#mkdir temp/
#cd temp/
#ls ../*.pm | xargs -i perl /work2/elia/src/ensembl/misc-scripts/utilities/dbsql_dbold.pl {}
#mv -f *.pm ../
#rm -rf temp/

my $file=shift(@ARGV);
open (FILE,"<$file");
$file =~ /(\w+\.pm)/;
my $new = $1;
open (NEW,">$new");
while (<FILE>) {
    $line =$_; 
    $line =~ s/DBSQL/DBOLD/g; 
    print NEW $line;
}
