use CPAN;

my $FTP = find_path('ftp');
my $GZIP = find_path('gzip');
my $LYNX = find_path('lynx');
my $MAKE = find_path('make');
my $NCFTPGET = find_path('ncftpget');
my $NCFTP = find_path('ncftp');
my $LESS = find_path('less');
my $TAR = find_path('tar');
my $UNZIP = find_path('unzip');
my $WGET = find_path('wget');

my $conf = bless {}, CPAN::Config;
$CPAN::Config->{'build_cache'} = 10;
$CPAN::Config->{'build_dir'} = "$ENV{'HOME'}/.cpan/build";
$CPAN::Config->{'cache_metadata'} = 1;
$CPAN::Config->{'cpan_home'} = "$ENV{'HOME'}/.cpan";
$CPAN::Config->{'ftp'} = $FTP;
$CPAN::Config->{'ftp_proxy'} = q[];
$CPAN::Config->{'getcwd'} = q[cwd];
$CPAN::Config->{'gzip'} = $GZIP;
$CPAN::Config->{'http_proxy'} = q[];
$CPAN::Config->{'inactivity_timeout'} = q[0];
$CPAN::Config->{'index_expire'} = q[1];
$CPAN::Config->{'inhibit_startup_message'} = q[0];
$CPAN::Config->{'keep_source_where'} = "$ENV{'HOME'}/.cpan/sources";
$CPAN::Config->{'lynx'} = $LYNX;
$CPAN::Config->{'make'} = $MAKE;
$CPAN::Config->{'make_arg'} = q[];
$CPAN::Config->{'make_install_arg'} = q[];
$CPAN::Config->{'makepl_arg'} = q[];
$CPAN::Config->{'ncftpget'} = $NCFTPGET;
$CPAN::Config->{'ncftp'} = $NCFTP;
$CPAN::Config->{'no_proxy'} = q[];
$CPAN::Config->{'pager'} = $LESS;
$CPAN::Config->{'prerequisites_policy'} = q[follow];
$CPAN::Config->{'scan_cache'} = q[atstart];
$CPAN::Config->{'shell'} = $ENV{'SHELL'};
$CPAN::Config->{'tar'} = $TAR;
$CPAN::Config->{'term_is_latin'} = 1;
$CPAN::Config->{'unzip'} = $UNZIP;
$CPAN::Config->{'urllist'} = [];
$CPAN::Config->{'wait_list'} = [q[wait://ls6.informatik.uni-dortmund.de:1404]];
$CPAN::Config->{'wget'} = $WGET;

$conf->commit($ARGV[0]);


sub find_path {
    my $prog = shift;
    my $path = `which $prog`;
    chomp $path;
    $path =~ s/\033.*\007//;
    if ($path =~ /(no $prog in|Comman not found)/) {
	return;
    } else {
	return $path;
    }
}
