use strict;
use GD;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use SeqContigDraw;
use Parameters;


my $image_param_ref=&Parameters::gene_image_par;
my $im=new GD::Image($image_param_ref->{x_img_len},$image_param_ref->{y_img_len});


my $dbuser = 'ensembl';
my $dbname = 'db11';
my $host = 'sol28';
my $dbpass = undef; 
my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my $gene_id=shift @ARGV;$gene_id || die 'I need contig id e.g. ENSG00000021276 ';
my $gene=$db->get_Gene($gene_id);



print "<body bgcolor=\"#FFFFFF\">\n";
print  "<IMG SRC=\"contig.gif\" BORDER=\"0\" USEMAP=\"#Ngs22\"><BR>\n";
print  "<MAP Name=\"Ngs22\">\n";


&SeqContigDraw::draw_gene($im,$gene);


open(GIF,">contig.gif");
print GIF $im->gif;
close(GIF);


print  "</map>\n";


