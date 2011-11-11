#!/bin/perl
# 1st Feb 2011
# Generate an HTML documentation page from an SQL file.
#
# It needs to have a "javascript like" documentation above each table. e.g.:
####################################################################################
#/**
#@table variation

#@desc This is the schema's generic representation of a variation.

#@colour #FF0000

#@column variation_id				Primary key, internal identifier.
#@column source_id					Foreign key references to the @link source table.
#@column name								Name of the variation. e.g. "rs1333049".
#@column validation_status	Variant discovery method and validation from dbSNP.
#@column ancestral_allele		Taken from dbSNP to show ancestral allele for the variation.
#@column flipped						This is set to 1 if the variant is flipped.
#@column class_so_id				Class of the variation, based on the Sequence Ontology.

#@see variation_synonym
#@see flanking_sequence
#@see failed_variation
#@see variation_feature
#@see variation_group_variation
#@see allele
#@see allele_group_allele
#@see individual_genotype_multiple_bp
#*/
#
#
#create table variation (
#		variation_id int(10) unsigned not null auto_increment, # PK
#		source_id int(10) unsigned not null, 
#		name varchar(255),
#		validation_status SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious'),
#		ancestral_allele text,
#		flipped tinyint(1) unsigned NULL DEFAULT NULL,
#		class_so_id ENUM('SO:0001483','SO:1000002','SO:0000667') DEFAULT 'SO:0001059', # default to sequence_alteration
#
#		primary key( variation_id ),
#		unique ( name ),
#		key source_idx (source_id)
#);
#
#/**
#@legend #FF0000 Table storing variation data
#*/
#
########################################################################################
# Tags description:
# /** and */ : begin and end of the document block
# @header    : tag to create a group of tables
# @table     : name of the sql table
# @desc      : description of the role/content of the table, set or info tags
# @colour    : tag to colour the header of the table (e.g. if the tables are coloured in the graphic SQL schema and you want to reproduce it in the HTML version)
# @column    : column_name [tab(s)] Column description. Note: 1 ligne = 1 column
# @see       : tables names linked to the described table
# @link      : Internal link to an other table description. The format is ... @link table_name ...
# @info      : tag to describe additional information about a table or a set of tables
# @legend    : tag to fill the colour legend table at the end of the HTML page

use strict;
use POSIX;
use Getopt::Long;

###############
### Options ###
###############
my ($sql_file,$html_file,$db_team,$show_colour,$header_flag,$sort_headers,$sort_tables,$help);

usage() if (!scalar(@ARGV));
 
GetOptions(
    'i=s' => \$sql_file,
    'o=s' => \$html_file,
		'd=s' => \$db_team,
		'c=i' => \$show_colour,
    'show_header=i' => \$header_flag,
		'sort_headers=i' => \$sort_headers,
    'sort_tables=i' => \$sort_tables,
		'help!' => \$help
);

usage() if ($help);

if (!$sql_file) {
	print "> Error! Please give a sql file using the option '-i' \n";
	usage();
}
if (!$html_file) {
	print "> Error! Please give an output file using the option '-o'\n";
	usage();
}

$show_colour  = 1 if (!defined($show_colour));
$header_flag  = 1 if (!defined($header_flag));
$sort_headers = 1 if (!defined($sort_headers));
$sort_tables  = 1 if (!defined($sort_tables));


##############
### Header ###
##############
my $html_header = qq{
<html>
<head>
<meta http-equiv="CONTENT-TYPE" content="text/html; charset=utf-8" />
<title>$db_team Schema Documentation</title>

<script language="Javascript" type="text/javascript">
	// Function to show/hide the columns table
	function show_hide (param) {
		div   = document.getElementById('div_'+param);
		alink = document.getElementById('a_'+param);
		if (div.style.display=='inline') {
			div.style.display='none';
			alink.innerHTML='Show';
		}
		else {
			if (div.style.display=='none') {
				div.style.display='inline';
				alink.innerHTML='Hide';
			}
		}
	}
</script>

</head>

<body>
<h1>Ensembl $db_team Schema Documentation</h1>

<h2>Introduction</h2>

<p><i>
&lt;please, insert your introduction here&gt;
</i></p>
<br />
};


##############
### Footer  ##
##############
my $html_footer = qq{
</body>
</html>};



################
### Settings  ##
################
my $default_colour = '#000000'; # Black

my %display_col = ('Show' => 'none', 'Hide' => 'inline');
my $documentation = {};
my $tables_names = {'default' => []};
my @header_names = ('default');
my @colours = ($default_colour);
my %legend;


my $in_doc = 0;
my $in_table = 0;
my $header = 'default';
my $table = '';
my $info = '';
my $nb_by_col = 15;
my $count_sql_col = 0;
my $tag_content = '';
my $tag = '';
my $display = 'Show';
my $parenth_count = 0;
my $header_colour;

#############
## Parser  ##
#############

open SQLFILE, "< $sql_file" or die "Can't open $sql_file : $!";
while (<SQLFILE>) {
	chomp $_;
	next if ($_ eq '');
	next if ($_ =~ /^\s*(DROP|PARTITION)/i);
	
	# Verifications
	if ($_ =~ /^\/\*\*/)  { $in_doc=1; next; }  # start of a table documentation
	if ($_ =~ /^\s*create\s+table\s+(if\s+not\s+exists\s+)?(\w+)/i) { # start to parse the content of the table
		my $sql_t_name = remove_char($2);
		if ($sql_t_name eq $table) { 
			$in_table=1;
			$parenth_count++;
		}
		else { 
			print STDERR "The documentation of the table $sql_t_name has not be found!\n";
		}
		next;
	}	
	next if ($in_doc==0 and $in_table==0);
	
	my $doc = $_;
	
	## Parsing of the documentation ##
	if ($in_doc==1) {
		# Header name
		if ($doc =~ /^\@header\s*(.+)$/i and $header_flag == 1) {
			$header = $1;
			push (@header_names,$header);
			$tables_names->{$header} = [];
			next;
		}		
		# Table name
		elsif ($doc =~ /^\@table\s*(\w+)/i) {
			$table = $1;
			push(@{$tables_names->{$header}},$table);
			$documentation->{$header}{'tables'}{$table} = { 'desc' => '', 'colour' => '', 'column' => [], 'see' => [], 'info' => [] };
			$tag = $tag_content = '';		
		}
		# Description (used for both set, table and info tags)
		elsif ($doc =~ /^\@(desc)\s*(.+)$/i) {
			fill_documentation ($1,$2);
		}
		# Colour of the table header (used for both set, table) (optional)
		elsif ($doc =~ /^\@(colour)\s*(.+)$/i and $show_colour) {
			fill_documentation ($1,$2);
		}
		# Column
		elsif ($doc =~ /^\@(column)\s*(.+)$/i) {
			fill_documentation ($1,$2);
		}
		# See other tables
		elsif ($doc =~ /^\@(see)\s*(.+)$/i) {
			fill_documentation ($1,$2);	
		}
		# Addtional information block
		elsif ($doc =~ /^\@(info)\s*(.+)$/i) {
			fill_documentation ();
			$info = $2;
			next;
		}
		# End of documentation
		elsif ($doc =~ /^\*\//) { # End of the documentation block
			fill_documentation (); # Add the last tag content to the documentation hash
			$in_doc=0;
			next; 
		}
		# Add legend colour description
		elsif ($doc =~ /^\@(legend)\s*(#\S+)\s+(.+)$/i) {
			$legend{$2} = $3;
		}
		elsif ($doc =~ /^\s*(.+)$/) { # If a tag content is split in several lines
			$tag_content .= " $1";
		}
	}
	
	## Parsing of the SQL table to fetch the columns types ##
	elsif ($in_table==1) {
	
	  #END OF TABLE DEFINITION
	  #Can't do this easily with a simply regex as there are varying valid formats
	  #The end of the table definition is actually defined by 2nd enclosing bracket
	  
	  #Regex counting VOODOO!
	  #This basically puts the regex in a list context
	  #before inc/dec'ing with it in a scalar context.
	  $parenth_count +=()= $doc =~ /\(/gi;
	  $parenth_count -=()= $doc =~ /\)/gi;


	 # warn "parenth count is $parenth_count";

	  if ($parenth_count == 0) { # End of the sql table definition
		if (scalar @{$documentation->{$header}{'tables'}{$table}{column}} > $count_sql_col) {
		  
		  #use Data::Dumper;
		  #warn "col count $count_sql_col";
		  #warn Data::Dumper::Dumper(\$documentation);

		  print STDERR "Description of a non existant column in the table $table!\n";
		}

		$in_table=0;
		$count_sql_col = 0;
		$table='';
		$parenth_count = 0;
	  }
	  else{

	#	warn "Processing table";

		## INDEXES ##
		if ($doc =~ /^\s*(primary\s+key)\s*\((.+)\)/i or $doc =~ /^\s*(unique)\s*\((.+)\)/i){ # Primary or unique
			my $icol = remove_char($2);
			add_column_index($1,$icol);
			next;
		}
		elsif ($doc =~ /^\s*(unique\s+)?(key|index)\s+(\S+)\s*\((.+)\)/i) { # Keys and indexes
			my $icol = remove_char($4);
			add_column_index("$1$2",$icol,$3);
			next;
		}
		elsif ($doc =~ /^\s*(key)\s+\((.+)\)/i) { # Keys
			my $icol = remove_char($2);
			add_column_index("$1",$icol,'');
			next;
		}
		
		## TYPES & DEFAULT VALUES ##
		my $col_name = '';
		my $col_type = '';
		my $col_def  = '';
		
		
		# All the type is contained in the same line (type followed by parenthesis)
		if ($doc =~ /^\W*(\w+)\W+(\w+\s?\(.*\))/ ){
			$col_name = remove_char($1);
			$col_type = $2;
			if ($doc =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
		
		# The type is written in several lines
		elsif ($doc =~ /^\W*(\w+)\W+(enum|set)(\s?\(.*)/i){ # The content is split in several lines
			$col_name= remove_char($1);
			$col_type="$2$3<br />";
			my $end_type = 0;
			while ($end_type != 1){
				my $line = <SQLFILE>;
				chomp $line;

				#Regex counting VOODOO again
				$parenth_count +=()= $line =~ /\(/gi;
				$parenth_count -=()= $line =~ /\)/gi;
				
				if ($line =~ /\)/) { # Close parenthesis
					$end_type=1; 
					$line =~ /^\s*(.+)\)/;
					$col_type .= "$1)"; 
				}
				else { # Add the content of the line
					$line =~ /^\s*(.+)/;
					$col_type .= $1.'<br />';
				}
				if ($line =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
			}
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
		
		# All the type is contained in the same line (type without parenthesis)
		elsif ($doc =~ /^\s*\W*(\w+)\W+(\w+)/ ){
			$col_name = remove_char($1);
			$col_type = $2;
			if ($doc =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1;} # Default value
			add_column_type_and_default_value($col_name,$col_type,$col_def);
		}
	  }
	}
}
close(SQLFILE);


############
### Core  ##
############

my $html_content;

# Sort the headers names by alphabetic order
if ($sort_headers == 1) {
	@header_names = sort(@header_names);
}
# Sort the tables names by alphabetic order
if ($sort_tables == 1) {
	while ( my($header_name,$tables) = each (%{$tables_names})) {
		@{$tables} = sort(@{$tables});
	}
}

# Legend link
if ($show_colour and scalar @colours > 1) {
	$html_content .= qq{A colour legend is available at the <a href="#legend">bottom of the page</a>.};
}

# List of tables by header
$html_content .= display_tables_list();
my $table_count = 1;
my $col_count = 1;

# Display the detailled tables by header
foreach my $header_name (@header_names) {
	my $tables = $tables_names->{$header_name};
	my $hcolour = $documentation->{$header_name}{'colour'};
		
	# Header display	
	if ($header_flag == 1 and $header_name ne 'default') {
		$html_content .= qq{<br />\n<hr />\n<h2>$header_name</h2>\n};
		my $header_desc = $documentation->{$header_name}{'desc'};		
		$html_content .= qq{<p>$header_desc</p>} if (defined($header_desc));
	}
	elsif ($header_flag == 0 or scalar @header_names == 1) {
		$html_content .= qq{<hr />\n};
	}
	
	# Additional information
	if ($header_name eq 'default' and defined($documentation->{$header_name}{'info'})) {
		$html_content .= qq{<h2>Additional information about the schema</h2>\n};
	}	
	$html_content .= add_info($documentation->{$header_name}{'info'});
	
	# Tables display
	foreach my $t_name (@{$tables}) {
		my $data = $documentation->{$header_name}{'tables'}{$t_name};
		my $colour;
		
		if ($header_flag && $hcolour) {
			$colour = $hcolour;
		} else {
			$colour = $data->{colour};
		}
		
		$html_content .= add_table_name($t_name,$colour);
		$html_content .= add_description($data->{desc});
		$html_content .= add_info($data->{info});	
		$html_content .= add_columns($t_name,$data);
		$html_content .= add_see($data->{see});
	}
}
$html_content .= add_legend();


## HTML/output file ##
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML $html_content."\n";
print HTML $html_footer."\n";
close(HTML);




###############
##  Methods  ##
###############

sub display_tables_list {
	my $html = qq{<h3 id="top">List of the tables:</h3>\n};
	foreach my $header_name (@header_names) {
		
		my $tables = $tables_names->{$header_name};
		my $count = scalar @{$tables};
		next if ($count == 0);
		
		my $nb_col = ceil($count/$nb_by_col);
		my $table_count = 0;
		my $col_count = 1;
	
		if ($nb_col>3) { 
			while ($nb_col>3) {
				$nb_by_col += 5;
				$nb_col = ceil($count/$nb_by_col);
			}
			$nb_col = 3;
		}
		
		# Header
		if ($header_flag == 1 and $header_name ne 'default') {
			$html .= display_header($header_name);
		}
		$html .= qq{<table><tr><td>\n <ul>\n};

		# Table
		foreach my $t_name (@{$tables}) {
			if ($table_count == $nb_by_col and $col_count<$nb_col and $nb_col>1){
				$html .= qq{	</ul></td><td><ul>\n};
				$table_count = 0;
			}
			$html .= add_table_name_to_list($t_name);
			$table_count ++;
		}
		$html .= qq{		</ul>\n</td></tr></table>\n};
	}
	return $html;
}


sub display_header {
	my $header_name = shift;
	my $html;
	
	if ($show_colour && $header_colour) {
		my $hcolour = $documentation->{$header_name}{colour};
		$hcolour = $default_colour if (!$hcolour);
		$html .= qq{
		  <table style="border: 1px solid #CCCCCC;padding:0px;margin:0px;background-color:#FAFAFF"><tr>
		    <td style="background-color:$hcolour;width:10px"></td>
        <td><h2 style="margin:5px;">$header_name</h2></td>
		  </tr></table>\n
		};
	}
	else {
		$html .= qq{<h2>$header_name</h2>\n};
	}
	return $html;
}


# If the line starts by a @<tag>, the previous tag content is added to the documentation hash.
# This method allows to describe the content of a tag in several lines.
sub fill_documentation {
	my $t1 = shift;
	my $t2 = shift;
	
	if ($tag ne '') {
		# Description tag (info, table or header)
		if ($tag eq 'desc') {
			# Additional description
			if ($info ne '') {
				$tag_content = $info.'@info@'.$tag_content;
				# Table: additional information				
				if ($table ne '') {
					push(@{$documentation->{$header}{'tables'}{$table}{'info'}},$tag_content);
				}
				# Header: additional information
				else {
					if (!$documentation->{$header}{'info'}) {
						$documentation->{$header}{'info'} = [];
					}
					push(@{$documentation->{$header}{'info'}},$tag_content);
				}
				$info = '';
			}
			# Table description
			elsif ($table ne '') {
				$documentation->{$header}{'tables'}{$table}{$tag} = $tag_content;
			}
			# Header description
			else {
				$documentation->{$header}{'desc'} = $tag_content;
			}
		}
		elsif ($tag eq 'colour') {
			if ($table ne '') {
				$documentation->{$header}{'tables'}{$table}{$tag} = $tag_content;
				if (! grep {$tag_content eq $_} @colours) {
					push (@colours,$tag_content);
				}
			}
			else {
				$documentation->{$header}{'colour'} = $tag_content;
				$header_colour = 1;
			}
		}
		elsif ($tag eq 'column') {
			$tag_content =~ /(\w+)[\s\t]+(.+)/;
			
			my $column = { 'name'    => $1,
								     'type'    => '',
			               'default' => '',
								     'index'   => '',
						         'desc'    => $2
							     };
			push(@{$documentation->{$header}{'tables'}{$table}{$tag}},$column);
		}
		else{
			push(@{$documentation->{$header}{'tables'}{$table}{$tag}},$tag_content);
		}
	}
	# New tag initialized
	if ($t1) {
		$tag = $t1;
		$tag_content = $t2;
	}
	else {
		$tag = $tag_content = '';	
	}
}
 

sub add_table_name_to_list {
	my $t_name = shift;
	
	my $html = qq{		<li><a href="#$t_name"><b>$t_name</b></a></li>\n};
	return $html;
}

sub add_table_name {
	my $t_name = shift;
	my $colour = shift || '#000000';
	
	my $c_box = '';
	if ($show_colour) {
		$c_box = qq{
      <td style="padding:0px;width:10px;background-color:$colour"></td>
      <td style="width:2px"></td>};
	}
	
	my $html = qq{\n<br />
  <table style="border: 2px groove #CCCCCC;background-color:#FAFAFF">
    <tr style="vertical-align:middle">$c_box
      <td style="width:420px;text-align:left;height:10px"><span id="$t_name" style="font-size:11pt;font-weight:bold">$t_name</span></td>
      <td style="width:180px;text-align:right"><a id="a_$t_name" style="cursor:pointer;text-decoration:underline" onclick="show_hide('$t_name')">Show</a> columns
				&nbsp;<b>|</b>&nbsp;<a href="#top">[Back to top]</a></td>
    </tr>
  </table>\n};
	
	return $html;
}


sub add_description {
	my $desc = shift;
	return qq{<p>$desc</p>\n};
}

sub add_info {
	my $infos = shift;
	my $html = '';
	
	foreach my $inf (@{$infos}) {
		my ($title,$content) = split('@info@', $inf);
		$html .= qq{
    <table>
      <tr class="bg3"><th>$title</th></tr>
      <tr class="bg1"><td>$content</td></tr>
    </table>\n};
	}
	
	return $html;
}

sub add_columns {
	my $table = shift;
	my $data  = shift;
	my $cols  = $data->{column};
	my $display_style = $display_col{$display};
	
	my $html = qq{\n	<div id="div_$table" style="display:$display_style">\n
	<table style="border:1px outset #222222">
		<tr class="bg3 center"><th style="width:180px">Column</th><th style="width:150px">Type</th><th style="width:100px">Default value</th><th style="width:400px">Description</th><th style="width:150px">Index</th></tr>\n};
	my $bg = 1;
	foreach my $col (@$cols) {
		my $name    = $col->{name};
		my $type    = $col->{type};
		my $default = $col->{default};
		my $desc    = $col->{desc};
		my $index   = $col->{index};
		
		# links
		$desc =~ /\@link\s?(\w+)/;
		my $link = $1;
		if (!grep {$link eq $_} @{$data->{see}} and defined($link)) {
			push @{$data->{see}}, $link;
		}
		my $table_to_link = qq{<a href="#$link">$link</a>};
		$desc =~ s/\@link\s?\w+/$table_to_link/;
		
		#$col =~ /^\s*(\w+)[\s\t]+(.+)\t+(.+)\t(.*)/;
		$html .= qq{		<tr class="bg$bg"><td><b>$name</b></td><td>$type</td><td>$default</td><td>$desc</td><td>$index</td></tr>\n};
		if ($bg==1) { $bg=2; }
		else { $bg=1; }
	}
	$html .= qq {</table>\n</div>\n};
	
	return $html;
}

sub add_see {
	my $sees = shift;
	my $html = '';
	
	if (scalar @$sees) {
		$html .= qq{<p><b>See also:</b></p>\n<ul>\n};
		foreach my $see (@$sees) {
			$html .= qq{	<li><a href="#$see">$see</a></li>\n};
		}
		$html .= qq{</ul>\n};
	}
	
	return $html;
}


sub add_column_index {
	my $idx_type = shift;
	my $idx_col  = shift;
	my $idx_name = shift;
	
	my $index = $idx_type;
	$idx_name = remove_char($idx_name);
	if (defined($idx_name)) {
		$index .= ": $idx_name";
	}
	
	my @idx_cols = split(',',$idx_col); # The index can involve several columns
	
	my %is_found = ();
	foreach my $i_col (@idx_cols) {
		$i_col =~ s/\s//g; # Remove white spaces
		$i_col = remove_char($i_col); # Remove the ` character
		# In case of index using a number characters for a column
		if ($i_col =~ /(.+)\(\d+\)/) {
			$i_col = $1;
		} 
		$is_found{$i_col} = 0;
		foreach my $col (@{$documentation->{$header}{tables}{$table}{column}}) {
			if ($col->{name} eq $i_col) {
				if ($col->{index} ne '') {
					$col->{index} .= '<br />';
				}
				$col->{index} .= lc($index);
				$is_found{$i_col} = 1;
				last;
			}
		}
	}
	# Description missing
	while (my ($k,$v) = each(%is_found)) {
		if ($v==0) {
			print STDERR "INDEX: The description of the column '$k' is missing in the table $table!\n";
		}
	}
}


sub add_column_type_and_default_value {
	my $c_name    = shift;
	my $c_type    = shift;
	my $c_default = shift;
	$count_sql_col ++;
	
	my $is_found = 0;
	foreach my $col (@{$documentation->{$header}{'tables'}{$table}{column}}) {
		if ($col->{name} eq $c_name) {
			$col->{type} = $c_type;
			$col->{default} = $c_default if ($c_default ne '');
			$is_found = 1;
			last;
		}
	}
	# Description missing
	if ($is_found==0) {
		print STDERR "COLUMN: The description of the column '$c_name' is missing in the table $table!\n";
	}
}


sub add_legend {
	my $html = '';
	my $default = 'Other tables';
	
	return $html if (scalar @colours == 1 or $header_colour);
	
	$html .= qq{<br />\n<hr />\n<h3 id="legend">Colour legend</h3>\n<table>};
	
	foreach my $c (@colours) {
		my $desc = '';
		if ($c eq $default_colour && !$legend{$default_colour}) {
			$desc = $default;
		}
		else {
			$desc = $legend{$c};
		}
		$html .= qq{<tr><td style="width:25px;height:15px;background-color:$c"></td><td>$desc</td></tr>\n};
	}
	$html .= '</table>';
	
	return $html;
}


sub remove_char {
	my $text = shift;
	$text =~ s/`//g;
	return $text;
}


sub usage {
	
  print qq{
  Usage: perl sql2html.pl [OPTION]
  
  Convert the SQL documentation into an HTML document.
	
  Options:

    -help           Print this message
      
    An input file must be specified. This file must be a SQL file, with the "Java-doc like documentation". 
    For more information, please visit the following page: 
    http://www.ebi.ac.uk/seqdb/confluence/display/EV/SQL+documentation

    -i              A SQL file name (Required)
    -o              An HTML output file name (Required)
    -d              The name of the database (e.g Core, Variation, Functional Genomics, ...)
    -c              A flag to display the colours associated with the tables (1) or not (0).
                    By default, the value is set to 1.
    -show_header    A flag to display headers for a group of tables (1) or not (0). 
                    By default, the value is set to 1.
    -sort_headers   A flag to sort (1) or not (0) the headers by alphabetic order.
                    By default, the value is set to 1.
    -sort_tables    A flag to sort (1) or not (0) the tables by alphabetic order.
                    By default, the value is set to 1.							 
  } . "\n";
  exit(0);
}
