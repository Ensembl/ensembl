#!/usr/bin/perl

# $Id$
#
# GUZZLE
#
# A stand-alone perl-based web DAS-client
#
#
# Author:   Andreas Kähäri (EMBL-EBI/ensembl),
#	    andreas.kahari@ebi.ac.uk

use strict;
use warnings;

# Where you keep non-standard modules (including Bioperl)
use lib qw(/opt/local/libdata/perl5/i386-openbsd/5.8.0
	   /opt/local/libdata/perl5/i386-openbsd
	   /opt/local/libdata/perl5/site_perl/i386-openbsd
	   /opt/local/libdata/perl5/site_perl
	   /opt/local/libdata/perl5);

use Bio::Das;
use CGI::Pretty qw(:standard -compile);

use Data::Dumper;

#---------------------------------------------------------------
# Configurable things (change these):

# $contact:
# Who's responsible for this thing (webmaster or similar).
#
my $contact = 'Yours Truly<br>yt@example.com';

# $htdocs:
# Where the web server looks when someone requests the web
# server root ('/') URL.
#
my $htdocs  = '/var/www/htdocs';

# $tmpurl:
# The URL to the directory which hold the generated images.
#
my $tmpurl  = '/guzzle_tmp';

# $tmpdir:
# Where the directory of $tmpurl is located on the system.
# MAKE SURE THIS DIRECTORY EXISTS AND THAT IT IS WRITABLE AND
# EXECUTABLE BY THE WEB SERVER USER!
#
my $tmpdir  = $htdocs . $tmpurl;
if (! (-d $tmpdir && -w $tmpdir && -x $tmpdir)) {
    die "Check value of \$tmpdir, '$tmpdir' is not a writable directory\n";
}

# $query_page_title and $result_page_title:
# These holds the title stings for the two pages.
#
my $query_page_title	= 'Guzzle query page';
my $result_page_title	= 'Guzzle result page';

# @presets:
# The DSN's to present by default on the query page.  A DSN with
# a non-zero CHECKED field will be checked (used/queried) by
# default.
#
my @presets = (
    {	NAME	=> 'Sprot (using Ensembl peptide IDs and ' .
		   'simple ID mapping) [test]',
	DSN	=> 'http://uhuru/cgi-bin/das/sprot',
	MAPFILE	=> '/home/ak/ensembl-cvs/ensembl/misc-scripts/' .
		   'das_client/ensembl_sprot.simple_map',
	MAPTYPE	=> 'simple',
	CHECKED	=> 1 },
    {	NAME	=> 'Sprot1 (using Ensembl peptide IDs and ' .
		   'alignment mapping) [test]',
	DSN	=> 'http://uhuru/cgi-bin/das/sprot1',
	MAPFILE	=> '/home/ak/ensembl-cvs/ensembl/misc-scripts/' .
		   'das_client/ensembl_sprot.align_map',
	MAPTYPE	=> 'align',
	CHECKED	=> 1 },
    {	NAME	=> 'Sprot [test]',
	DSN	=> 'http://uhuru/cgi-bin/das/sprot',
	MAPTYPE	=> 'none',
	CHECKED	=> 0 },
    {	NAME	=> 'Pfam-A [test]',
	DSN	=> 'http://uhuru/cgi-bin/das/pfam',
	MAPTYPE	=> 'none',
	CHECKED	=> 0 } );

# $nblanks:
# The number of blank/empty fields to present on the query page
# by default.
#
my $nblanks = 1;

# @colours:
# The colours that the user should be able to
# choose from.  Must be "named colors", see e.g.
# "http://webdesign.about.com/library/bl_namedcolors.htm".
#
my @colours = qw( red green blue magenta cyan yellow gray black );

# $stylesheet:
# A CSS stylesheet that will be included in every generated web
# page.
#
my $stylesheet = <<EOT;
body {
    color:	    #000;
    background:	    #6a6;
    margin:	    5%;
    font-family:    helvetica, arial, sans-serif;
}
td {
    font-size:	    x-small;
}
h1 {
    font-family:    helvetica, arial, sans-serif;
}
address {
    font-family:    times, serif;
    color:	    #666;
    font-size:	    small;
}
tt {
    font-family:    courier, monospace;
}
.thetable {
    background:	    #9c9;
}
EOT

#
# OPTIONAL FEATURES
#

# $use_stylesheets
# Whether to try to make use of DAS stylesheets or not
# (non-zero is "yes").
#
my $use_stylesheets = 1;

# $be_nice
# Whether to be nice to DAS servers or not.
#
my $be_nice = 1;

# $use_graphics
# Whether to generate and use graphics or not.
#
my $use_graphics = 1;

# $use_*_mapping
# Whether to make available and use mappings between
# e.g. Ensembl and Swissprot.
#
my $use_simple_mapping = 1;
my $use_align_mapping = 1;

# No servicable parts inside...
#---------------------------------------------------------------

if ($use_stylesheets) {
    require LWP::UserAgent;
}

if ($be_nice) {
    require Data::Serializer;
}

if ($use_graphics) {
    require Bio::Graphics;
    require Bio::SeqFeature::Generic;

    require File::Temp;
    import File::Temp 'tempfile';
}

if ($use_align_mapping) {
    require Storable;
    import Storable 'dclone';

    require Guzzle::Mapper;
}

my $use_mapping = ($use_align_mapping || $use_simple_mapping);

#---------------------------------------------------------------

# do_query():
# Will perform the query.
#
sub do_query
{
    my $cgi = shift;
    my $sources = shift;

    my $start = $cgi->param('START') || undef;
    my $stop  = $cgi->param('STOP')  || undef;
    my $segment;

    my $seqid = $cgi->param('ID');

    my $range;
    if (defined $start && defined $stop) {
	$range = ':' . $start . ',' . $stop;
    } else {
	$range = '';
    }

    my @replies;
    foreach my $source (@{ $sources }) {
	my %query;

	if ($source->{MAPTYPE} eq 'none') {
	    # Don't do mapping for this source.
	    $query{$seqid}{SEGMENT} = $seqid . $range;
	    push(@{ $query{$seqid}{DSN} }, $source->{DSN});
	} elsif ($source->{MAPTYPE} eq 'simple') {

	    # Simple 1:N ID mapping:
	    #
            # Notation: The ID requested by the user is the
            #           "query ID".  The ID actually used by
            #           the client to perform the query is the
            #           "target ID".
	    #
            # 1. Find the target ID(s) that corresponds to
            #    the requested query ID in the map file.  The
            #    query ID should be found in column 1, and the
            #    corresponding target ID(s) will then be in
            #    column 2.  The map file is a comma-delimited
            #    file and the query IDs in column 1 need not be
            #    unique if one query ID corresponds to more than
            #    one target ID.
	    #
	    # 2. Perform the query.
	    #
            # 3. Tag all occurances of the target ID in the
            #    reply with the query ID.
	    #

	    open(IN, $source->{MAPFILE}) or
		die "Can't open map data file '" .
		    $source->{MAPFILE} ."': " .  $!;

	    while (defined(my $line = <IN>)) {
		chomp $line;
		my ($qi, $ti) = split /,/, $line;

		last if ($qi gt $seqid);
		next if ($qi ne $seqid);

		$query{$ti}{SEGMENT} = $ti . $range;
		push(@{ $query{$ti}{DSN} }, $source->{DSN});
	    }

	    close IN;
	} elsif ($source->{MAPTYPE} eq 'align') {

	    # Alignment mapping:
	    #
            # Notation: The ID requested by the user is the
            #           "query ID".  The ID actually used by
            #           the client to perform the query is the
            #           "target ID".
	    #
            # 1. Find the target ID(s) that corresponds to
            #    the requested query ID in the map file.  The
            #    query ID should be found in column 1, and the
            #    corresponding target ID(s) will then be in
            #    column 2.  The map file is a comma-delimited
            #    file and the query IDs in column 1 need not be
            #    unique if one query ID corresponds to more than
            #    one target ID.
	    #
            #    Column 3 and 4 should contain the start
            #    of the alignment in the query and target
            #    coordinate systems respecitvely.  The last
            #    column should contain the Ensembl style cigar
            #    string (e.g. "20MD340M2I20M" or somesuch
            #    thing).
	    #
            # 2. Parse the cigar line and build a
            #    Guzzle::Mapper object from it.
	    #
            # 3. If a range was specified, map it into the
            #    target coordinate system.
	    #
	    # 4. Perform the query.
	    #
            # 5. Tag all occurances of the target ID in the
            #    reply with the query ID.
	    #
            # 6. Map the features back to the query coordinate
            #    system.
            #

	    open(IN, $source->{MAPFILE}) or
		die "Can't open map data file '" .
		    $source->{MAPFILE} ."': " .  $!;

	    while (defined(my $line = <IN>)) {
		chomp $line;
		my ($qi, $ti, $qab, $tab, $C) = split /,/, $line;

		last if ($qi gt $seqid);
		next if ($qi ne $seqid);

		my $map = new Guzzle::Mapper('queryCOORD', 'targetCOORD');

		my ($qpos, $tpos) = ($qab, $tab);

		# Parse the cigar string:
		while ($C =~ /(\d*)([MID])/g) {
		    my ($len, $op) = ($1, $2);
		    my ($qadd, $tadd) = (0, 0);

		    if ($len eq '') {
			$len = 1;
		    }

		    if ($op eq 'M') {
			$map->add_map_coordinates(
			    'queryID', $qpos, $qpos + $len, 1,
			    'targetID', $tpos, $tpos + $len);

			$qpos += $len;
			$tpos += $len;

		    } elsif ($op eq 'D') {
			$qpos += $len;
		    } elsif ($op eq 'I') {
			$tpos += $len;
		    } else {
			die "Unknown cigar string operation '$op'\n";
		    }
		}

		if ($range ne '') {
		    # Map the requested range.
		    my @mapped = $map->map_coordinates(
			'queryID', $start, $stop, 1,
			'queryCOORD');

		    foreach my $mapped (@mapped) {
			next if ($mapped->isa('Guzzle::Mapper::Gap'));
			$range = ':' . $mapped->start . ',' .  $mapped->end;
			push(@{ $query{$ti}{SEGMENT} }, $ti . $range);
		    }
		} else {
		    $query{$ti}{SEGMENT} = $ti . $range;
		}

		$query{$ti}{MAPPER} = $map;
		push(@{ $query{$ti}{DSN} }, $source->{DSN});
	    }

	    close(IN);

	} else {
	    die "Unknown mapping type: " .
		$source->{MAPTYPE} . "\n";
	}

	my $das = new Bio::Das(15);

	foreach my $query (values %query) {
	    foreach my $dsn (@{ $query->{DSN} }) {
		my $reply = $das->features(
		    -dsn	    => $dsn,
		    -segment    => $query->{SEGMENT});

		next if (!$reply->is_success);

		#print $cgi->pre(Dumper($reply));

		if (exists $query->{MAPPER} ) {
		    # Map results using align mapper.
		    # Scary stuff.  If this works, I deserve a beer.

		    my @results = $reply->results;
		    foreach my $result (@results) {
			# Tag the string that we later use in the tabla
			# and graphics.

			# Do the reverse mapping from the target
			# coodinate system into the query coordinate
			# system.  This will yield an array of "gaps"
			# and "coordinates".  If a range is mapped as a
			# whole (no holes) then just change the start
			# and stop coodinates.  Tag the result with
			# "mapped", "fragmented", and "unmappable" as
			# neccesary.

			my @mapped = $query->{MAPPER}->map_coordinates(
			    'targetID', $result->start, $result->stop, 1,
			    'targetCOORD');

			if (scalar @mapped > 1) {
			    for my $i (0 .. (scalar @mapped - 1)) {

				my $resultcopy = dclone($result);

				$resultcopy->start($mapped[$i]->start);
				$resultcopy->stop($mapped[$i]->end);

				$resultcopy->group($result->group .
				    " [fragmented in $seqid]");

				# The following if-statment and
				# push-call breaks the OO badly since
				# it assumes that the underlying
				# representation of the objects are
				# known, but it was the only way I could
				# make it work.

				$resultcopy->{segment} = $result->{segment};

				if ($mapped[$i]->
					isa('Guzzle::Mapper::Gap')) {
				    $resultcopy->{type}{label} .=
					" [fragment " .  (1 + $i) .
					" (not mapped)]";
				} else {
				    $resultcopy->{type}{label} .=
					" [fragment " .  (1 + $i) . "]";
				}

				push(@{ $reply->{results}[1] }, $resultcopy);
			    }
			    # Will take it off the display
			    $result->group('NOSHOW');

			} elsif ($mapped[0]->isa('Guzzle::Mapper::Gap')) {
			    $result->group($result->group .
				" [unmappable in $seqid]");
			} else {
			    $result->group($result->group .
				" [mapped to $seqid]");
			    $result->start($mapped[0]->start);
			    $result->stop($mapped[0]->end);
			}
		    }
		}
		push(@replies, $reply);
	    }
	}
    }

    return @replies;
}

#---------------------------------------------------------------

# page_start_and_head():
# Will create the head of a HTML page with the CSS stylesheet
# from above ($stylesheet) and a H1 header with the page title.
#
sub page_start_and_head
{
    my $cgi = shift;
    my $title = shift;

    print $cgi->header( -expires => 'now' );

    print $cgi->start_html(
	-title => $title,
	-style => { -code => $stylesheet });

    print $cgi->h1({ style => 'text-align:center' }, $title);
}

# page_foot_and_end():
# Will put the contact address in the page footer and then
# close/finish the page.
#
sub page_foot_and_end
{
    my $cgi = shift;

    print $cgi->p({ style => 'text-align:right' }, 
	$cgi->small($cgi->tt('$Revision$')));

    print $cgi->hr, $cgi->address($contact), $cgi->end_html;
}

#---------------------------------------------------------------

# create_graphics():
# Given the array of structures created in result_page(), will
# create a PNG image and return its URL.
#
sub create_graphics
{
    my $table = shift;	# "Table" since it's the data that goes
			# into the table on the result page.

    my $do_descr = shift;   # Whether to include the descriptions
			    # on each feature in the graph or not.

    my $minpos = $table->[0]{FEATURE}->start;
    my $maxpos = $table->[0]{FEATURE}->stop;

    my $full_seq;
    my %all;
    foreach my $table_row (@{ $table }) {
	my $feature = $table_row->{FEATURE};

	if ($minpos > $feature->start) {
	    $minpos = $feature->start;
	}
	if ($maxpos < $feature->stop) {
	    $maxpos = $feature->stop;
	}

	my $tag	    = $table_row->{COLUMNS}[0];
	my $source  = $table_row->{COLUMNS}[1];

	if ($feature->group ne 'NOSHOW') {
	    if ($feature->group =~ /^Sequence:/) {
		if (!defined $full_seq) {
		    $full_seq = new Bio::SeqFeature::Generic(
			-start	=> $feature->start,
			-end	=> $feature->end );
		}
	    } else {
		push(@{ $all{$source}->{$tag} },
		    new Bio::SeqFeature::Generic(
			-start	    => $feature->start,
			-end	    => $feature->end,
			-primary    => $feature->group,
			-source_tag => $feature->type->label,
			-tag	    => {
			    colour  => $table_row->{COLOUR}
			} ));
	    }
	}
    }

    my $panel = new Bio::Graphics::Panel(
	-start	    => $minpos,
	-stop	    => $maxpos,
	-grid	    => 1,
	-key_style  => 'between',
	-pad_left   => 10,
	-pad_right  => 100,
	-pad_top    => 10,
	-pad_bottom => 10);

    my $full_length = new Bio::SeqFeature::Generic(
	-start	=> $minpos,
	-end	=> $maxpos);

    # Add a nice arrow at the top of the image.
    $panel->add_track($full_length,
	-glyph	    => 'arrow',
	-tick	    => 2,
	-fgcolor    => 'black',
	-double	    => 1);

    # Populate the panel.
    foreach my $source (values %all) {
	foreach my $tag (keys %{ $source }) {
	    if (defined $do_descr && $do_descr == 1) {
		$panel->add_track($source->{$tag},
		    -key	    => $tag,
		    -description    => sub { $_[0]->source_tag },
		    -bgcolor	    => $source->{$tag}[0]->_tag_value('colour'),
		    -font2color	    => 'black' );
	    } else {
		$panel->add_track($source->{$tag},
		    -key	=> $tag,
		    -bgcolor	=> $source->{$tag}[0]->_tag_value('colour'));
	    }
	}
    }

    # Get a handle of a temporary (persistent) file.
    my ($tmpfh, $tmpname) = tempfile('pdas_XXXXX',
	UNLINK	=> 0,
	SUFFIX	=> '.png',
	DIR	=> $tmpdir);

    # Dump the image data into the file.
    syswrite($tmpfh, $panel->png);

    # Make the image file readable for web clients.
    chmod(0644, $tmpname);

    $tmpname =~ s#.*/##;    # Does a 'basename' operation.

    return $tmpurl . '/' . $tmpname;
}

#---------------------------------------------------------------

# query_page():
# Creates the page with the query form.
#
sub query_page
{
    my $cgi = shift;

    if (defined $cgi->param('COUNT')) {
	# We've been here before, so figure out how many blank
	# lines we should display.
	$nblanks = $cgi->param('COUNT') - scalar(@presets);
    }

    if (defined $cgi->param('INCR')) {
	# User wants to insert more blank lines.
	++$nblanks;
    } elsif (defined $cgi->param('DECR')) {
	# User wants to remove blank lines.
	--$nblanks;
    }

    $nblanks = 0 if ($nblanks < 0);

    page_start_and_head($cgi, $query_page_title);

    $cgi->delete_all;

    print $cgi->start_form;

    print $cgi->start_table({
	-class		=> 'thetable',
	-border		=> '0',
	-cellpadding	=> '2',
	-cellspacing	=> '2',
	-align		=> 'center' });
    if ($use_stylesheets) {
	print $cgi->Tr($cgi->th( [ qw( Use Name URL Colour ),
	    '... or Stylesheet' ] ));
    }

    my $count = 0;
    foreach my $preset (@presets) {
	print $cgi->Tr($cgi->td( { -align => 'center' },
	    $cgi->checkbox(
		-name	    => 'USE',
		-value	    => $count,
		-checked    => $preset->{CHECKED},
		-label	    => '')),
	    $cgi->td( [
		$preset->{NAME},
		$preset->{DSN},
		$cgi->popup_menu(
		    -name	=> 'COLOUR',
		    -values	=> \@colours,
		    -default	=> $colours[$count % scalar @colours])
		]),
	    ($use_stylesheets ?
	    $cgi->td( { -align => 'center' },
		$cgi->checkbox(
		    -name	=> 'STYLE',
		    -value	=> $count,
		    -checked	=> '1',
		    -label	=> '' )) : '' ) );  # That wasn't pretty...
	print $cgi->hidden(
	    -name   => 'NAME',
	    -value  => $preset->{NAME});
	print $cgi->hidden(
	    -name   => 'DSN',
	    -value  => $preset->{DSN});
	print $cgi->hidden(
	    -name   => 'MAPTYPE',
	    -value  => $preset->{MAPTYPE}) if exists($preset->{MAPTYPE});
	print $cgi->hidden(
	    -name   => 'MAPFILE',
	    -value  => $preset->{MAPFILE}) if exists($preset->{MAPFILE});
	++$count;
    }

    for (my $i = 0; $i < $nblanks; ++$i) {
	print $cgi->Tr($cgi->td( { -align => 'center' },
	    $cgi->checkbox(
		-name	    => 'USE',
		-value	    => $count,
		-checked    => '0',
		-label	    => '')),
	    $cgi->td( [
		$cgi->textfield(
		    -name	=> 'NAME',
		    -size	=> '10',
		    -default	=> 'Source #' . (1 + $count)),
		$cgi->textfield(
		    -name	=> 'DSN',
		    -size	=> '50'),
		$cgi->popup_menu(
		    -name	=> 'COLOUR',
		    -values	=> \@colours,
		    -default	=> $colours[$count % scalar @colours])
		] ),
	    ($use_stylesheets ?
	    $cgi->td( { -align => 'center' },
		$cgi->checkbox(
		    -name	=> 'STYLE',
		    -value	=> $count,
		    -checked	=> '0',
		    -label	=> '' )) : '' ) );  # Stylish? No...
	++$count;
    }
    print $cgi->Tr($cgi->td( [ '&nbsp;', '&nbsp;',
	$cgi->submit(
	    -name   => 'INCR',
	    -value  => 'More blank fields') . '&nbsp;' .
	$cgi->submit(
	    -name   => 'DECR',
	    -value  => 'Fewer blank fields'), '&nbsp;',
	    ($use_stylesheets ? '&nbsp;' : '') ] ));

    print $cgi->hidden(
	-name	=> 'COUNT',
	-value	=> $count);
    print $cgi->hidden(
	-name	=> 'DESCR',
	-value	=> 1);

    print $cgi->Tr($cgi->th(
	{ -colspan => '2', -align => 'right' }, 'Sequence ID'),
	$cgi->td({ -colspan => '2' }, $cgi->textfield(
	    -name	=> 'ID',
	    -size	=> '25',
	    -default	=> 'ENSP00000186985') ),
	($use_stylesheets ? $cgi->td('&nbsp;') : ''));


    print $cgi->Tr($cgi->th(
	{ -colspan => '2', -align => 'right' }, 'Range'),
	$cgi->td($cgi->textfield(
	    -name	=> 'START',
	    -size	=> '5') . ' to ' .
	$cgi->textfield(
	    -name	=> 'STOP',
	    -size	=> '5')),
	$cgi->td({ -align => 'center',
	    -colspan => ($use_stylesheets ? '2' : '1') },
	    $cgi->submit( -name => 'GO', -value => '   Go   ' ) ));

    print $cgi->end_table, $cgi->end_form;

    page_foot_and_end($cgi);
}

#---------------------------------------------------------------

# result_page():
# Queries the DAS server, generates the graphics, and creates
# the result page.
#
sub result_page
{
    my $cgi = shift;

    my $cereal;
    if ($be_nice) {
        # Create a serialization/deserialization object
        # ($cereal) that we use to embed the results from the
        # query in a variable (RESULT) in the form on the result
        # page.  This means that if the user wants to re-sort
        # the result, we won't need to query the DAS server
        # again.

	$cereal = new Data::Serializer(
	    serializer  => 'Storable',
	    compress    => 1,
	    portable    => 1);
    }

    page_start_and_head($cgi, $result_page_title);

    my $table = [ ];
    if (!defined $cgi->param('RESULTS') || !$be_nice) {
	# The results are fetched from the DAS server.

	my @sources;
	foreach my $idx ($cgi->param('USE')) {
	    my %source = (
		NAME	    => [ $cgi->param('NAME')    ]->[$idx],
		DSN	    => [ $cgi->param('DSN')	]->[$idx],
		COLOUR	    => [ $cgi->param('COLOUR')  ]->[$idx],
		USESTYLE    => [ $cgi->param('STYLE')	]->[$idx],
		MAPTYPE	    => [ $cgi->param('MAPTYPE')	]->[$idx],
		MAPFILE	    => [ $cgi->param('MAPFILE')	]->[$idx] );

	    $source{DSN} =~ s/\s//g;
	    next unless (length $source{DSN} > 0);

	    push(@sources, \%source);
	}

	my @replies = do_query($cgi, \@sources);

	foreach my $reply (@replies) {
	    next unless ($reply->is_success && defined $reply->results);

	    my $source_url = $reply->dsn->base;
	    if (defined $reply->dsn->id) {
		$source_url .= '/' . $reply->dsn->id;
	    }
	    $source_url =~ s/\/features\?.*//;

	    foreach my $source (@sources) {
		next if ($source_url ne $source->{DSN});

                # FIXME:  This is a quick'n dirty hack to work
                # around the non-existant Bio::Das::Stylesheet
                # module (in the CPAN distribution of the
                # Bio::Das modules, as of 2003-09-22).  It will
                # pick up and use the *first* FGCOLOR from the
                # stylesheet, if one exists.

		if ($use_stylesheets && defined $source->{USESTYLE}) {
		    my $ua		= new LWP::UserAgent;
		    my $response	= $ua->get($source_url . '/stylesheet');
		    if ($response->is_success) {
			my $headers = $response->headers;
			if ($headers->header('x-das-status') == 200) {
			    $response->content =~ m#<FGCOLOR>(\w+)</FGCOLOR>#;
			    $source->{COLOUR} = $1 if (defined $1);
			}
		    }
		}

		foreach my $feature ($reply->results) {
		    my %table_row = (
			COLUMNS => [
			    $feature->group,	    # Label
			    $source->{NAME},	    # Source
			    $feature->type->label,  # Description
			    $feature->start,	    # Start
			    $feature->stop,	    # Stop
			    $feature->score ],	    # Score
			FEATURE => $feature,	# Yes, this will duplicate
						# some of the data...
			COLOUR	=> $source->{COLOUR} );
		    push(@{ $table }, \%table_row);
		}
		last;
	    }
	}

	if (scalar @replies == 0 || scalar @{ $table } == 0) {
	    print $cgi->b("Sorry, no features were found for your query.");
	    page_foot_and_end($cgi);
	    exit;
	}

    } else {
        # The results are picked up from the encoded string,
        # not from any DAS server.  Hopefully this works
        # everywhere...

	$table = $cereal->deserialize($cgi->param('RESULTS'));
    }

    print $cgi->start_form;

    my $sort1;
    my $sort2;

    if (defined $cgi->param('SORT1')) {
	$sort1 = $cgi->param('SORT1');
	$sort2 = $cgi->param('SORT2');
    } else {
	$sort1 = 1;
	$sort2 = 3;
    }

    # Sorting.  Must distinguish between the numerical and the
    # non-numerical columns (the last few ones and the first few
    # ones, respectively).

    # Column 5 (score) is a special case.

    if ($sort1 == 5) {
	$table = [ sort {
	    substr($a->{COLUMNS}[2], 0, 1 + index($a->{COLUMNS}[2], ':')) cmp
	    substr($b->{COLUMNS}[2], 0, 1 + index($b->{COLUMNS}[2], ':')) ||
	    $a->{COLUMNS}[5] <=> $b->{COLUMNS}[5] } @{ $table } ];
    } else {
	if ($sort1 <= 2) {
	    if ($sort2 <= 2) {
		$table = [ sort {
		    $a->{COLUMNS}[$sort1] cmp
		    $b->{COLUMNS}[$sort1] ||
		    $a->{COLUMNS}[$sort2] cmp
		    $b->{COLUMNS}[$sort2] } @{ $table } ];
	    } else {
		$table = [ sort {
		    $a->{COLUMNS}[$sort1] cmp
		    $b->{COLUMNS}[$sort1] ||
		    $a->{COLUMNS}[$sort2] <=>
		    $b->{COLUMNS}[$sort2] } @{ $table } ];
	    }
	} else {
	    if ($sort2 <= 2) {
		$table = [ sort {
		    $a->{COLUMNS}[$sort1] <=>
		    $b->{COLUMNS}[$sort1] ||
		    $a->{COLUMNS}[$sort2] cmp
		    $b->{COLUMNS}[$sort2] } @{ $table } ];
	    } else {
		$table = [ sort {
		    $a->{COLUMNS}[$sort1] <=>
		    $b->{COLUMNS}[$sort1] ||
		    $a->{COLUMNS}[$sort2] <=>
		    $b->{COLUMNS}[$sort2] } @{ $table } ];
	    }
	}
    }

    if ($be_nice) {
	# Embed the results for fast re-sort.

	$cgi->autoEscape(0);
	print $cgi->hidden(
	    -name	=> 'RESULTS',
	    -value	=> $cereal->serialize($table));
	$cgi->autoEscape(1);
    } else {
	# Save the state from the query page.
	print $cgi->hidden(
	    -name => 'USE',
	    -value => $cgi->param('USE') );
	print $cgi->hidden(
	    -name => 'NAME',
	    -value => $cgi->param('NAME') );
	print $cgi->hidden(
	    -name => 'DSN',
	    -value => $cgi->param('DSN') );
	print $cgi->hidden(
	    -name => 'COLOUR',
	    -value => $cgi->param('COLOUR') );
	print $cgi->hidden(
	    -name => 'STYLE',
	    -value => $cgi->param('STYLE') );
	print $cgi->hidden(
	    -name => 'START',
	    -value => $cgi->param('START') );
	print $cgi->hidden(
	    -name => 'STOP',
	    -value => $cgi->param('STOP') );
	print $cgi->hidden(
	    -name => 'ID',
	    -value => $cgi->param('ID') );
    }

    print $cgi->start_table({
	-class		=> 'thetable',
	-border		=> '0',
	-cellpadding	=> '2',
	-cellspacing	=> '2',
	-align		=> 'center' });

    if ($use_graphics) {
	my $imgsrc = create_graphics($table, $cgi->param('DESCR'));

	print $cgi->Tr($cgi->td( { -colspan => '6', -align => 'center' },
	    $cgi->img({ -src => $imgsrc, -alt => 'Result (PNG)' }),
	    $cgi->p({ -style => 'text-align:right' },
	    $cgi->checkbox(
		-name	    => 'DESCR',
		-value	    => 1,
		-checked    => 1,
		-label	    => ' Show descriptions on tracks' ),
	    $cgi->br,
	    $cgi->submit( -name => 'GO', -value => 'Update display' )),
	    $cgi->p('&nbsp;')
	     ));
    }

    print $cgi->Tr($cgi->th( [
	qw( &nbsp; Label Source Description Start Stop Score ) ] ));

    foreach my $table_row (@{ $table }) {
	# Don't display the full sequence reply (good/bad?)
	next if ($table_row->{COLUMNS}[0] =~ /^Sequence:/);
	# ... or if it is set to 'NOSHOW'
	next if ($table_row->{COLUMNS}[0] eq 'NOSHOW');

	# Make the first column the link.
	$table_row->{COLUMNS}[0] = $cgi->a({
	    -href => $table_row->{FEATURE}->link },
	    $table_row->{COLUMNS}[0]);

	print $cgi->Tr($cgi->td( {
	    -style => 'background:' . $table_row->{COLOUR} }, '&nbsp;' ),
	    $cgi->td($table_row->{COLUMNS}));
    }

    print $cgi->Tr($cgi->td({ -colspan => 6 }, '&nbsp;' ));
    print $cgi->Tr($cgi->th({ -colspan => 6 }, 'Sorting' ));

    $cgi->autoEscape(0);
    print $cgi->Tr($cgi->td({ -colspan => 2, -align => 'right' },
	'First sort on '),
	$cgi->td({ -colspan => 2 },
	$cgi->radio_group(
	    -name	=> 'SORT1',
	    -values	=> [ 0, 1, 2, 3, 4, 5 ],
	    -labels	=> {
		0	=> ' Label',
		1	=> ' Source',
		2	=> ' Description',
		3	=> ' Start',
		4	=> ' Stop',
		5	=> ' Score' },
	    -default	=> 1 ) ),
	$cgi->td({
	    -rowspan	=> 2,
	    -colspan	=> 2,
	    -align	=> 'center',
	    -valign	=> 'center' },
	$cgi->submit( -name => 'GO', -value => 'Update display' ) ));

    print $cgi->Tr($cgi->td({ -colspan => 2, -align => 'right' },
	'Then sort on '),
	$cgi->td({ -colspan => 2 },
	$cgi->radio_group(
	    -name	=> 'SORT2',
	    -values	=> [ 0, 1, 2, 3, 4, 5 ],
	    -labels	=> {
		0	=> ' Label',
		1	=> ' Source',
		2	=> ' Description',
		3	=> ' Start',
		4	=> ' Stop',
		5	=> ' Score' },
	    -default	=> 3 ) ));
    $cgi->autoEscape(1);
    print $cgi->end_table, $cgi->end_form;

    page_foot_and_end($cgi);
}

#---------------------------------------------------------------

my $cgi = new CGI;

if (defined $cgi->param('GO')) {
    result_page($cgi);
} else {
    query_page($cgi);
}
