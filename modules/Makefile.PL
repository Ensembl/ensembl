# Let the code begin...

require 5.003;
use ExtUtils::MakeMaker;

$NAME     = 'Bio';
$DISTNAME = "ensembl";
$VERSION  = "0.05";

WriteMakefile(
    NAME	=> $NAME,
    DISTNAME    => $DISTNAME,
    VERSION     => $VERSION,
    'dist'      => {  COMPRESS     => 'gzip -9f',
                      SUFFIX       => '.gz', 
                      DIST_DEFAULT => 'all tardist',
                   },

   'realclean'  => {  FILES        => join(' ',@BACKUP_FILES) }
);
