BEGIN{

    unshift @INC, "../../modules";
    unshift @INC,"../../../bioperl-ensembl-06/";
}



use strict;
use GD;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use SeqContigDraw;
use Parameters;




my $image_param_ref=&Parameters::contig_image_par;
my $im=new GD::Image($image_param_ref->{x_img_len},$image_param_ref->{y_img_len});

my $dbuser = 'ensro';
my $dbname = 'ensembl';
my $host = 'ensrv3.sanger.ac.uk';
my $dbpass = undef; 
my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

my $contig_id=shift @ARGV; $contig_id || die 'I need contig id e.g. AC004471.00001  ';
my $contig=$db->get_Contig($contig_id);


print "<body bgcolor=\"#FFFFFF\">\n";
print  "<IMG SRC=\"contig.gif\" BORDER=\"0\" USEMAP=\"#Ngs22\"><BR>\n";
print  "<MAP Name=\"Ngs22\">\n";


&SeqContigDraw::draw_contig_image($im,$contig);

open(GIF,">contig.gif");
print GIF $im->gif;
close(GIF);

print  "</map>\n";



