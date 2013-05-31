#!/usr/bin/env perl
# 1st Feb 2011
# Generate an HTML documentation page from an SQL file.
#
# It needs to have a "javascript like" documentation above each table.
# See the content of the method sql_documentation_format();
####################################################################################


use strict;
use POSIX;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


###############
### Options ###
###############

my ($sql_file,$html_file,$db_team,$show_colour,$version,$header_flag,$format_headers,$sort_headers,$sort_tables,$intro_file,$help,$help_format);
my ($host,$port,$dbname,$user,$pass,$skip_conn,$db_handle);

usage() if (!scalar(@ARGV));
 
GetOptions(
    'i=s' => \$sql_file,
    'o=s' => \$html_file,
    'd=s' => \$db_team,
    'c=i' => \$show_colour,
    'v=i' => \$version,
    'show_header=i'    => \$header_flag,
    'format_headers=i' => \$format_headers,
    'sort_headers=i'   => \$sort_headers,
    'sort_tables=i'    => \$sort_tables,
    'host=s'           => \$host,
    'port=i'           => \$port,
    'dbname=s'         => \$dbname,
    'user=s'           => \$user,
    'pass=s'           => \$pass,
    'skip_connection'  => \$skip_conn,
    'intro=s'          => \$intro_file,
    'help!'            => \$help,
    'help_format'      => \$help_format,
);

usage() if ($help);
sql_documentation_format() if ($help_format);


if (!$sql_file) {
  print "> Error! Please give a sql file using the option '-i' \n";
  usage();
}
if (!$html_file) {
  print "> Error! Please give an output file using the option '-o'\n";
  usage();
}

$show_colour    = 1 if (!defined($show_colour));
$header_flag    = 1 if (!defined($header_flag));
$format_headers = 1 if (!defined($format_headers));
$sort_headers   = 1 if (!defined($sort_headers));
$sort_tables    = 1 if (!defined($sort_tables));

$skip_conn      = undef if ($skip_conn == 0);


# Dababase connection (optional)
if (defined($host) && !defined($skip_conn)) {
  my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host => $host,
    -user => $user,
    -pass => $pass,
    -port => $port,
    -dbname => $dbname
  ) or die("DATABASE CONNECTION ERROR: Could not get a database adaptor for $dbname on $host:$port\n");
  $db_handle = $db_adaptor->dbc->db_handle;
}




################
### Settings  ##
################

my $default_colour = '#000'; # Black
my $list_bg = "background-color:#F2F2F2";

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

my $SQL_LIMIT = 50;
my $img_plus  = qq{<img src="/i/16/plus-button.png" style="width:12px;height:12px;vertical-align:middle" alt="show"/>};
my $img_minus = qq{<img src="/i/16/minus-button.png" style="width:12px;height:12px;vertical-align:middle;margin-right:6px" alt="hide"/>};
my $link_text = 'columns';




##############
### Header ###
##############

my $html_header = qq{
<html>
<head>
<meta http-equiv="CONTENT-TYPE" content="text/html; charset=utf-8" />
<title>$db_team Schema Documentation</title>

<script language="Javascript" type="text/javascript">
  var img_plus   = '$img_plus';
  var img_minus  = '$img_minus';
  var link_text  = 'columns';
  var span_open  = ' <span style="vertical-align:middle">';
  var span_close = '</span>';

  // Function to show/hide the columns table
  function show_hide (param, type) {
    // Example tables
    if (type === 'example') {
      div    = document.getElementById('ex_'+param);
      alink  = document.getElementById('e_'+param);
      a_text = 'query results'; 
    } 
    // Schema tables
    else { 
      div    = document.getElementById('div_'+param);
      alink  = document.getElementById('a_'+param);
      a_text = link_text;
    }  
    
    if (div.style.display=='inline') {
      div.style.display='none';
      alink.innerHTML=img_plus+span_open+'Show '+a_text+span_close;
    }
    else {
      if (div.style.display=='none') {
        div.style.display='inline';
        alink.innerHTML=img_minus+span_open+'Hide '+a_text+span_close;
      }
    }
  }
  
  // Function to show/hide all the tables
  function show_hide_all () {
    expand_flag = document.getElementById('expand');
    divs = document.getElementsByTagName('div');
    for(var i=0; i < divs.length; i++){
      param = divs[i].id.substring(4);
      div   = document.getElementById('div_'+param);
      alink = document.getElementById('a_'+param);
      
      if (div && alink) {
        if (expand_flag.value=='0') {
          div.style.display='inline';
          alink.innerHTML=img_minus+span_open+'Hide '+link_text+span_close;
        }
        else {
          div.style.display='none';
          alink.innerHTML=img_plus+span_open+'Show '+link_text+span_close;
        }
      }
    }
    if (expand_flag.value=='0') {
      expand_flag.value='1';
    }
    else {
      expand_flag.value='0';
    }
  }
</script>

</head>
<body>
};


##############
### Footer  ##
##############

my $html_footer = qq{
</body>
</html>};




#############
## Parser  ##
#############

# Create a complex hash "%$documentation" to store all the documentation content

open SQLFILE, "< $sql_file" or die "Can't open $sql_file : $!";
while (<SQLFILE>) {
  chomp $_;
  next if ($_ eq '');
  next if ($_ =~ /^\s*(DROP|PARTITION)/i);
  next if ($_ =~ /^\s*(#|--)/); # Escape characters
  
  # Verifications
  if ($_ =~ /^\/\*\*/)  { $in_doc=1; next; }  # start of a table documentation
  if ($_ =~ /^\s*create\s+table\s+(if\s+not\s+exists\s+)?`?(\w+)`?/i) { # start to parse the content of the table
    if ($table eq $2) { 
      $in_table=1;
      $parenth_count++;
    }
    else { 
      print STDERR "The documentation of the table $2 has not be found!\n";
    }
    next;
  }  
  next if ($in_doc==0 and $in_table==0);
  
  my $doc = remove_char($_);
  
  #================================================#
  # Parsing the documentation part of the SQL file #
  #================================================#
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
      $documentation->{$header}{'tables'}{$table} = { 'desc' => '', 'colour' => '', 'column' => [], 'example' => [], 'see' => [], 'info' => [] };
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
    # Example 
    elsif ($doc =~ /^\@(example)\s*(.+)$/i) {
      fill_documentation ($1,$2);
    }
    # See other tables
    elsif ($doc =~ /^\@(see)\s*(\w+)\s*$/i) {
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
  
  #=====================================================#
  # Parsing of the SQL table to fetch the columns types #
  #=====================================================#
  elsif ($in_table==1) {
    
    # END OF TABLE DEFINITION
    # Can't do this easily with a simply regex as there are varying valid formats
    # The end of the table definition is actually defined by 2nd enclosing bracket
    
    # Regex counting VOODOO!
    # This basically puts the regex in a list context
    # before inc/dec'ing with it in a scalar context.
    $parenth_count +=()= $doc =~ /\(/gi;
    $parenth_count -=()= $doc =~ /\)/gi;


    if ($parenth_count == 0) { # End of the sql table definition
    if (scalar @{$documentation->{$header}{'tables'}{$table}{column}} > $count_sql_col) {
      print STDERR "Description of a non existant column in the table $table!\n";
    }

    $in_table=0;
    $count_sql_col = 0;
    $table='';
    $parenth_count = 0;
    }
    else {

      #---------#
      # INDEXES #
      #---------#
      
      # Skip the FOREIGN KEY
      next if ($doc =~ /^\s*foreign\s+key/i);
      
      if ($doc =~ /^\s*(primary\s+key)\s*\w*\s*\((.+)\)/i or $doc =~ /^\s*(unique)\s*\((.+)\)/i){ # Primary or unique;
        add_column_index($1,$2);
        next;
      }
      elsif ($doc =~ /^\s*(unique\s+)?(key|index)\s+([^\s\(]+)\s*\((.+)\)/i) { # Keys and indexes
        add_column_index("$1$2",$4,$3);
        next;
      }
      elsif ($doc =~ /^\s*(unique)\s+(\S*)\s*\((.+)\)/i) { # Unique
        add_column_index("$1",$3,$2);
        next;
      }
      elsif ($doc =~ /^\s*(key)\s+\((.+)\)/i) { # Keys
        add_column_index("$1",$2);
        next;
      }
      
      #----------------------------------#
      # COLUMNS & TYPES & DEFAULT VALUES #
      #----------------------------------#
      my $col_name = '';
      my $col_type = '';
      my $col_def  = '';
      
      # All the type is contained in the same line (type followed by parenthesis)
      if ($doc =~ /^\W*(\w+)\W+(\w+\s?\(.*\))/ ){
        $col_name = $1;
        $col_type = $2;
        if ($doc =~ /default\s+([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
        add_column_type_and_default_value($col_name,$col_type,$col_def);
      }
    
      # The type is written in several lines
      elsif ($doc =~ /^\W*(\w+)\W+(enum|set)(\s?\(.*)/i){ # The content is split in several lines
        $col_name= $1;
        $col_type="$2$3<br />";
        my $end_type = 0;
        while ($end_type != 1){
          my $line = <SQLFILE>;
          chomp $line;
          $line = remove_char($line);

          # Regex counting VOODOO again
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
          if ($line =~ /default\s+([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1; } # Default value
        }
        add_column_type_and_default_value($col_name,$col_type,$col_def);
      }
    
      # All the type is contained in the same line (type without parenthesis)
      elsif ($doc =~ /^\s*\W*(\w+)\W+(\w+)/ ){
        $col_name = $1;
        $col_type = $2;
        if ($doc =~ /default\s*([^,\s]+)\s*.*(,|#).*/i) { $col_def = $1;} # Default value
        add_column_type_and_default_value($col_name,$col_type,$col_def);
      }
    }
  }
}
close(SQLFILE);




###########
## Core  ##
###########

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
if ($show_colour and scalar @colours > 1 and $header_flag != 1) {
  $html_content .= qq{A colour legend is available at the <a href="#legend">bottom of the page</a>.\n<br /><br />};
}

#=============================================#
# List of tables by header (i.e. like a menu) #
#=============================================#
$html_content .= display_tables_list();
my $table_count = 1;
my $col_count = 1;

#========================================#
# Display the detailled tables by header #
#========================================#
foreach my $header_name (@header_names) {
  my $tables = $tables_names->{$header_name};
  my $hcolour = $documentation->{$header_name}{'colour'};
   
  #----------------#  
  # Header display #
  #----------------#  
  if ($header_flag == 1 and $header_name ne 'default') {
    $html_content .= qq{\n<br /><br />
<div style="$list_bg;padding:5px;margin:5px 0px;border-top:2px solid $hcolour">
  <h2 style="display:inline;color:#000">$header_name</h2>
</div>\n};
    my $header_desc = $documentation->{$header_name}{'desc'};    
    $html_content .= qq{<p style="width:800px">$header_desc</p>} if (defined($header_desc));
  }
  
  #------------------------#
  # Additional information #
  #------------------------#
  if ($header_name eq 'default' and defined($documentation->{$header_name}{'info'})) {
    $html_content .= qq{<h2>Additional information about the schema</h2>\n};
  }  
  $html_content .= add_info($documentation->{$header_name}{'info'});
  
  #----------------#
  # Tables display #
  #----------------#
  foreach my $t_name (@{$tables}) {
    my $data = $documentation->{$header_name}{'tables'}{$t_name};
    my $colour = ($header_flag && $hcolour) ? $hcolour : $data->{colour};
    
    $html_content .= add_table_name($t_name,$colour);
    $html_content .= add_description($data);
    $html_content .= add_info($data->{info},$data);  
    $html_content .= add_columns($t_name,$data);
    $html_content .= add_examples($t_name,$data);
    $html_content .= add_see($data->{see});
  }
}
$html_content .= add_legend();




######################
## HTML/output file ##
######################
open  HTML, "> $html_file" or die "Can't open $html_file : $!";
print HTML $html_header."\n";
print HTML slurp_intro($intro_file)."\n";
print HTML $html_content."\n";
print HTML $html_footer."\n";
close(HTML);




###############
##  Methods  ##
###############

# List the table names. 
# Group them if the header option "-format_headers" is selected.
# By default there is one group, named "default" and it contains all the tables.
sub display_tables_list {

  my $html; 
  
  $header_flag = 0 if (scalar @header_names == 1);
  
  if ($header_flag == 1) {
    $html .= qq{\n<h3 id="top">List of the tables:</h3>\n};
    $html .= qq{<div>\n} if ($format_headers == 1);
  } 
  else {
    my $list_width;
    if (scalar @header_names == 1) {
      my $list_count = scalar @{$tables_names->{'default'}};
      my $list_nb_col = ceil($list_count/$nb_by_col);
      $list_width = length_names($tables_names->{'default'},$list_nb_col);
    }
    $html .= qq{
<div>      
  <div id="top" style="$list_bg;border-radius:5px;margin-bottom:20px;float:left;$list_width">
    <div style="padding:5px;background-color:#336;border-top-left-radius:5px;border-top-right-radius:5px">
      <img src="/i/16/rev/info.png" style="vertical-align:top" />
      <h3 style="display:inline;color:#FFF">List of the tables:</h3>
    </div>};
  }
  
  my $has_header = 0;
  my $nb_col_line = 0;
  
  foreach my $header_name (@header_names) {
    
    my $tables = $tables_names->{$header_name};
    my $count = scalar @{$tables};
    next if ($count == 0);
    
    # Number of columns needed to display the tables of the group
    my $nb_col = ceil($count/$nb_by_col);
    my $nbc = $nb_col;
    my $table_count = 0;
    my $col_count = 1;
  
    if ($nb_col>3) { 
      while ($nb_col>3) {
        $nb_by_col += 5;
        $nb_col = ceil($count/$nb_by_col);
      }
      $nb_col = 3;
    }
    
    
    # Header #
    if ($header_flag == 1) {
      if ($header_name ne 'default') {
        if ($nb_col_line+$nbc > 4 and $format_headers == 1) {
          $html .= qq{  <div style="clear:both" />\n</div>\n\n<div>};
          $nb_col_line = 0;
        }
      
        $html .= display_header($header_name,$nbc);
        $nb_col_line += $nbc;
        $has_header = 1;
      }
      
      # List of tables #
      $html .= qq{\n      <div style="float:left">} if ($count > $nb_by_col);
      $html .= qq{\n      <ul style="padding:0px 4px 0px 22px;margin-bottom:2px">\n};
      my $t_count = 0;
      foreach my $t_name (@{$tables}) {
        my $t_colour;
        if ($has_header == 0 && $show_colour) {
          $t_colour = $documentation->{$header_name}{'tables'}{$t_name}{'colour'};
        }
        $html .= add_table_name_to_list($t_name,$t_colour);
        $t_count++;
        if ($t_count>=$nb_by_col) {
          $html .= qq{\n      </ul>\n      </div>};
          $html .= qq{\n      <div style="float:left">};
          $html .= qq{\n      <ul style="padding:0px 4px 0px 22px;margin-bottom:2px">\n};
          $t_count = 0;
        }
      }
      $html .= qq{\n      </ul>};
      $html .= qq{\n      </div>} if ($count > $nb_by_col);
      $html .= qq{\n    </div>\n} if ($format_headers == 1);   
    }
    else {
      $html .= qq{\n    <table style="padding:0px 2px"><tr><td>\n      <ul style="padding-left:20px">\n};

      # List of tables #
      foreach my $t_name (@{$tables}) {
        if ($table_count == $nb_by_col and $col_count<$nb_col and $nb_col>1){
          $html .= qq{      </ul>\n    </td><td>\n      <ul style="padding-left:20px">\n};
          $table_count = 0;
        }
        my $t_colour;
        if ($has_header == 0 && $show_colour) {
          $t_colour = $documentation->{$header_name}{'tables'}{$t_name}{'colour'};
        }
        $html .= add_table_name_to_list($t_name,$t_colour);
        $table_count ++;
      }
      $html .= qq{      </ul>\n    </td></tr></table>\n};
    }
  }
  
  my $input_margin;
  if ($header_flag == 1 and $format_headers == 1){
    $html .= qq{\n  <div style="clear:both" />\n</div>};
  } else {
    $input_margin = qq{ style="margin-left:10px;margin-bottom:5px"};
  }
  $html .= qq{
  <input type="button" onclick="show_hide_all()" class="fbutton" value="Show/hide all"$input_margin/>
  <input type="hidden" id="expand" value="0" />
  };
  
  $html .= qq{\n  </div>\n  <div style="clear:both" />\n</div>} if ($header_flag!=1 and $format_headers == 1);
  
  return $html;
}


# If the option "-show_header" is selected, the tables will be displayed by group ("Header") in the HTML page.
# This method generates the HTML code to display the group names & descriptions.
sub display_header {
  my $header_name = shift;
  my $nb_col = shift;
  
  my $html;
  
  if ($format_headers == 1) {
    my $width = length_names($tables_names->{$header_name},$nb_col);
    
    $html .= qq{\n  <div style="$list_bg;border-radius:5px;margin-bottom:15px;float:left;margin-right:20px$width">
    <div style="padding:2px 5px;background-color:#336;border-top-left-radius:5px;border-top-right-radius:5px">};
    
    my $text_colour = '#FFF';
    
    if ($show_colour && $header_colour) {
      my $hcolour = $documentation->{$header_name}{colour};
         $hcolour = $default_colour if (!defined($hcolour));
      
      $html .= qq{
      <div style="background-color:$hcolour;border:1px solid #FFF;padding:0px 8px;display:inline;vertical-align:middle"></div>
      <h2 style="margin-left:8px;display:inline;color:$text_colour;vertical-align:middle">$header_name</h2>\n};
    }
    else {
      $html .= qq{<h2 style="display:inline;color:$text_colour">$header_name</h2>\n};
    }
    $html .= qq{    </div>};
  } 
  else {
    $html .= qq{    <h2 >$header_name</h2>};
  }
  return $html;
}


# Method to pick up the documentation information contained in the SQL file.
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
      # Header description
      elsif(!$documentation->{$header}{'tables'}) {
        $documentation->{$header}{'desc'} = $tag_content;
      }
      # Table description
      else {
        $documentation->{$header}{'tables'}{$table}{$tag} = $tag_content;
      }
    }
    elsif ($tag eq 'colour') {
      if(!$documentation->{$header}{'tables'}) {
        $documentation->{$header}{'colour'} = $tag_content;
        $header_colour = 1;
      }
      elsif ($table ne '') {
        $documentation->{$header}{'tables'}{$table}{$tag} = $tag_content;
        if (! grep {$tag_content eq $_} @colours) {
          push (@colours,$tag_content);
        }
      }
    }
    elsif ($tag eq 'column') {
      $tag_content =~ /(\w+)[\s\t]+(.*)/;
      
      my $column = { 'name'    => $1,
                     'type'    => '',
                     'default' => '',
                     'index'   => '',
                     'desc'    => $2
                   };
      if ($2 eq '') {
        print STDERR "COLUMN: The description content of the column '$1' is missing in the table $table!\n";
      }
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
 

# Method generating the HTML code to display the table name into the top menu.
sub add_table_name_to_list {
  my $t_name = shift;
  my $t_colour = shift;
  my $c = $t_colour;
  if (defined($t_colour)) {
    $t_colour = ($t_colour ne '') ? qq{;background-color:$t_colour} : '';
    $t_colour = qq{<div style="padding:0px;margin-left:0px$t_colour;display:inline">&nbsp;</div> };
  }
  my $html = qq{        <li>$t_colour<a href="#$t_name"><b>$t_name</b></a></li>\n};
  return $html;
}


# Method generating the HTML code to display the title/header of the table description block
sub add_table_name {
  my $t_name = shift;
  my $colour = shift || '#000';
  
  my $c_box = '';
  if ($show_colour) {
    $c_box = qq{    <div style="float:left;padding:0px;height:20px;width:10px;background-color:$colour;margin-right:10px"></div>};
  }

  my $html = qq{
  <div id="$t_name" style="width:820px;height:20px;border: 2px groove #CCC;background-color:#FAFAFF;padding:2px;margin-top:35px;margin-bottom:2px">
    $c_box
    <div style="float:left;text-align:left;font-size:11pt;font-weight:bold">$t_name</div>
    <div style="float:right;text-align:right;padding-right:1px">
      <a id="a_$t_name" style="cursor:pointer;text-decoration:none" onclick="show_hide('$t_name')">
        $img_plus
        <span style="vertical-align:middle">Show $link_text</span>
      </a> 
      <b> | </b> <a href="#top">[Back to top]</a>
    </div>
  </div>\n};  
  
  return $html;
}


# Method generating the HTML code to display the description content
sub add_description {
  my $data = shift;
  
  # Search if there are some @link tags in the description text.
  my $desc = add_internal_link($data->{desc},$data);
  
  return qq{  <p style="padding:5px 0px;margin-bottom:0px;width:800px">$desc</p>\n};
}


# Method generating the HTML code to display additional information contained in the tags @info
sub add_info {
  my $infos = shift;
  my $data  = shift;
  my $html  = '';
  
  foreach my $inf (@{$infos}) {
    my ($title,$content) = split('@info@', $inf);
    $content = add_internal_link($content,$data) if (defined($data));
    
    $html .= qq{
    <table>
      <tr class="bg3"><th>$title</th></tr>
      <tr class="bg1"><td>$content</td></tr>
    </table>\n};
  }
  
  return $html;
}


# Method generating the HTML code of the table listing the columns of the given SQL table.
sub add_columns {
  my $table = shift;
  my $data  = shift;
  my $cols  = $data->{column};
  my $display_style = $display_col{$display};
  
  my $html = qq{\n  <div id="div_$table" style="display:$display_style">
    <table style="border:1px solid #667aa6;padding:0px;min-width:1000px;max-width:1200px">
      <tr class="center" style="color:#FFFFFF;background-color:#667aa6"><th style="color:#FFF;padding:2px">Column</th><th style="color:#FFF;padding:2px">Type</th><th style="color:#FFF;padding:2px;min-width:80px">Default value</th><th style="color:#FFF;padding:2px;min-width:500px">Description</th><th style="color:#FFF;padding:2px;min-width:100px">Index</th></tr>\n};
  my $bg = 1;
  
  foreach my $col (@$cols) {
    my $name    = $col->{name};
    my $type    = $col->{type};
    my $default = $col->{default};
    my $desc    = $col->{desc};
    my $index   = $col->{index};
    
    # links
    $desc = add_internal_link($desc,$data);
    
    $html .= qq{      <tr class="bg$bg"><td><b>$name</b></td><td>$type</td><td>$default</td><td>$desc</td><td>$index</td></tr>\n};
    if ($bg==1) { $bg=2; }
    else { $bg=1; }
  }
  $html .= qq {    </table>\n  </div>\n};
  
  return $html;
}


# Method generating the HTML code to display the content of the tags @example (description + SQL query + Table of SQL query results)
sub add_examples {
  my $table = shift;
  my $data  = shift;
  my $examples  = $data->{example};
  my $html;

  my $nb = (scalar(@$examples) > 1) ? 1 : '';

  foreach my $ex (@$examples) {
    my @lines = split("\n",$ex);
    my $nb_display = ($nb ne '') ? " $nb" : $nb;
    $html .= qq{<div style="margin: 10px 5px"><p style="font-weight:bold">Example$nb_display:</p>};
    my $has_desc = 0;
    my $sql;
    
    # Parse the example lines
    foreach my $line (@lines) {
      chomp($line);
      
      # Pick up the SQl query if it exists
      if ($line =~ /(.*)\s*\@sql\s*(.+)/) {
        $html .= ($has_desc == 1) ? $1 : qq{<p style="padding-left:10px">$1};
        $sql = $2;
      } elsif (!defined($sql)){
        $html .= qq{<p style="padding-left:10px">} if ($has_desc == 0);
        $html .= $line;
        $has_desc = 1;
      }
      else {
        $sql .= $line;
      }
    }
    $html .= qq{</p>};
    
    # Search if there are some @link tags in the example description.
    $html = add_internal_link($html,$data);
    
    # Add a table of examples
    if (defined($sql)) {
      my $show_hide = '';
      my $sql_table = '';
      if (!defined($skip_conn) && defined($host)) {
        $show_hide .= qq{<a id="e_$table$nb" style="cursor:pointer;text-decoration:none" onclick="show_hide('$table$nb','example')">$img_plus<span style="vertical-align:middle"> Show query results</span></a>};
        $sql_table = get_example_table($sql,$table,$nb);
      }
      if (defined($sql)) {
        foreach my $word (qw(SELECT DISTINCT CONCAT FROM LEFT JOIN USING WHERE LIMIT DESC ORDER GROUP)) {
          $sql =~ s/$word /$word /i;
        }
      }
      $html .= qq{<pre style="display:inline;border:1px solid #555;padding:2px;margin-right:15px;margin-left:10px">$sql</pre> $show_hide\n$sql_table};
    }
    $html .= qq{</div>};
    $nb ++;
  }
  
  return $html;
}


# Method generating the HTML code to display the content of the tags @see
sub add_see {
  my $sees = shift;
  my $html = '';
  
  if (scalar @$sees) {
    $html .= qq{  <p style="margin-top:10px"><b>See also:</b></p>\n  <ul>\n};
    foreach my $see (@$sees) {
      $html .= qq{    <li><a href="#$see">$see</a></li>\n};
    }
    $html .= qq{  </ul>\n\n};
  }
  
  return $html;
}


# Method searching the tag @link into the string given as argument and replace it by an internal HTML link 
sub add_internal_link {
  my $desc = shift;
  my $data = shift;
  while ($desc =~ /\@link\s?(\w+)/) {
    my $link = $1;
    if ((!grep {$link eq $_} @{$data->{see}}) and defined($link)) {
      push @{$data->{see}}, $link;
    }
    my $table_to_link = qq{<a href="#$link">$link</a>};
    $desc =~ s/\@link\s?\w+/$table_to_link/;
  }
  return $desc;
}


# Method parsing the index information from the SQL table description in order to display it in the
# HTML table listing the columns of the corresponding SQL table.
sub add_column_index {
  my $idx_type = shift;
  my $idx_col  = shift;
  my $idx_name = shift;
  
  my $index = $idx_type;
  if (!defined($idx_name)) {
    $idx_name = $idx_col;
  }
  if ($idx_type !~ /primary/i) {
    $index .= ": $idx_name";
  }
  my @idx_cols = split(',',$idx_col); # The index can involve several columns
  
  my %is_found = ();
  foreach my $i_col (@idx_cols) {
    $i_col =~ s/\s//g; # Remove white spaces
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


# Method parsing the column type and default value from the SQL table description, in order to display them in the
# HTML table listing the columns of the corresponding SQL table.
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


# Method to query the database with the SQL query example, get the result and display it
# in an HTML table.
sub get_example_table {
  my $sql   = shift;
  my $table = shift;
  my $nb    = shift;
  my $html;
  
  $sql =~ /select\s+(.+)\s+from/i;
  my $cols = $1;
  my @tcols;
     
  foreach my $col (split(',',$cols)) {
  
    # Columns selection like the expressions "table.*" or "*"
    my $table_name;
    $table_name = $table if ($cols eq '*');
    $table_name = $1 if ($col =~ /(\S+)\.\*/ and !defined($table_name));
    if (defined($table_name)) {
      my $table_cols = $db_handle->selectall_arrayref(qq{SHOW COLUMNS FROM $table_name});
      foreach my $col (@$table_cols) {
        push(@tcols,$col->[0]);
      }
      next;
    }
     
    # Check for alias
    $col = $1 if ($col =~ /\s+as\s+(\S+)$/i);
       
    $col =~ s/ //g;
    push(@tcols,$col);
  }
  
  my $results = $db_handle->selectall_arrayref($sql);
  if (scalar(@$results)) {
    $html .= qq{
  <div id="ex_$table$nb" style="display:none;">
    <table class="ss" style="width:50%;margin-top:20px;margin-left:10px;border-spacing:2px">\n      <tr><th>};
    $html .= join("</th><th>",@tcols);
    $html .= qq{</th></tr>};
    
    my $bg = '';
    
    my $count = 0;
    foreach my $result (@$results) {
      last if ($count >= $SQL_LIMIT);
      $html .= qq{      <tr$bg><td>};
      $html .= join("</td><td>", @$result);
      $html .= qq{</td></tr>};
      
      $bg = ($bg eq '')  ? ' class="bg2"' : '';  
      $count ++;
    }
    $html .= qq{    </table>\n  </div>};
  } else {
    my $msg = qq{ERROR: the SQL query displayed above returned no results!};
    $html .= qq{<div style="padding:5px;margin:10px;width:500px;font-weight:bold;border:2px solid red;color:red">$msg</div>};
    print STDERR qq{SQL: $sql\n$msg\n};
  }
  
  return $html;
}


# Method generating a "Colour legend" paragraph, based on the header colours.
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
    $html .= qq{  <tr><td style="width:25px;height:15px;background-color:$c"></td><td>$desc</td></tr>\n};
  }
  $html .= qq{</table>};
  
  return $html;
}


# Removed the character(s) ` from the read line.
sub remove_char {
  my $text = shift;
  $text =~ s/`//g;
  return $text;
}


# Count largest length in the table names
sub length_names {
  my $list = shift;
  my $nb_col = shift;
  
  my $max = 0;
  foreach my $name (@$list) {
    my $length = length($name);
    $max = $length if ($length>$max);
  }
  return '' if ($nb_col>1 || $max>25);
  return ";width:200px" if ($max<=25);
}

# Insert the introduction text of the web page
sub slurp_intro {
  my ($intro_file) = @_;
  if (!defined $intro_file) {
    return qq{<h1>Ensembl $db_team Schema Documentation</h1>\n<h2>Introduction</h2>\n<p><i>please, insert your introduction here</i><p><br />};
  }

  local $/=undef;
  open my $fh, "< $intro_file" or die "Can't open $intro_file: $!";
  my $intro_html = <$fh>;
  close $fh;
  
  $intro_html =~ s/####DB_VERSION####/$version/g if (defined($version));
  
  return $intro_html;
}


##################
## Help methods ##
##################

sub sql_documentation_format {
  print q{
  
#--------------------------#
# Example of documentation #  
#--------------------------#
  
/**
@table variation

@desc This is the schema's generic representation of a variation.

@colour #FF0000

@column variation_id       Primary key, internal identifier.
@column source_id          Foreign key references to the \@link source table.
@column name               Name of the variation. e.g. "rs1333049".
@column validation_status  Variant discovery method and validation from dbSNP.
@column ancestral_allele   Taken from dbSNP to show ancestral allele for the variation.
@column flipped            This is set to 1 if the variant is flipped.
@column class_so_id        Class of the variation, based on the Sequence Ontology.

@example Example of SQL query for to retrieve data from this table:
         @sql SELECT * FROM variation WHERE source_id=1 LIMIT 10;

@see variation_synonym
@see flanking_sequence
@see failed_variation
@see variation_feature
@see variation_group_variation
@see allele
@see allele_group_allele
@see individual_genotype_multiple_bp
*/


create table variation (
    variation_id int(10) unsigned not null auto_increment, # PK
    source_id int(10) unsigned not null, 
    name varchar(255),
    validation_status SET('cluster','freq','submitter','doublehit','hapmap','1000Genome','failed','precious'),
    ancestral_allele text,
    flipped tinyint(1) unsigned NULL DEFAULT NULL,
    class_so_id ENUM('SO:0001483','SO:1000002','SO:0000667') DEFAULT 'SO:0001059', # default to sequence_alteration

    primary key( variation_id ),
    unique ( name ),
    key source_idx (source_id)
);

/**
@legend #FF0000 Table storing variation data
*/


#========================================================================================================================#


#------------------#
# Tags description #
#------------------#

 /** and */ : begin and end of the document block
 @header    : tag to create a group of tables
 @table     : name of the sql table
 @desc      : description of the role/content of the table, set or info tags
 @colour    : tag to colour the header of the table (e.g. if the tables are coloured in the graphic SQL schema and you want to reproduce it in the HTML version)
 @column    : column_name [tab(s)] Column description. Note: 1 ligne = 1 column
 @see       : tables names linked to the described table
 @link      : Internal link to an other table description. The format is ... @link table_name ...
 @info      : tag to describe additional information about a table or a set of tables
 @legend    : tag to fill the colour legend table at the end of the HTML page
 @example   : tag to add some examples, like examples of SQL queries
 @sql       : tag inside the @example tag, used to delimit a SQL query

};
  exit(0);
}


sub usage {
  
  print q{
  Usage: perl sql2html.pl [OPTION]
  
  Convert the SQL documentation into an HTML document.
  
  Options:

    -help            Print this message
    -help_format     Print the description of the documentation format in the SQL files
      
    An input file must be specified. This file must be a SQL file, with the "Java-doc like documentation". 
    For more information, please visit the following page: 
    http://www.ebi.ac.uk/seqdb/confluence/display/EV/SQL+documentation

    -i                A SQL file name (Required)
    -o                An HTML output file name (Required)
    -d                The name of the database (e.g Core, Variation, Functional Genomics, ...)
    -c                A flag to display the colours associated with the tables (1) or not (0). By default, the value is set to 1.
    -v                Version of the schema. Replace the string ####DB_VERSION#### by the value of the parameter "-v", in the introduction text. (Optional)
    -intro            A html/text file to include in the Introduction section (Optional. If not provided a default text will be inserted)
    -show_header      A flag to display headers for a group of tables (1) or not (0). By default, the value is set to 1.
    -format_headers   A flag to display formatted headers for a group of tables (1) or not (0) in the top menu list. By default, the value is set to 1.                
    -sort_headers     A flag to sort (1) or not (0) the headers by alphabetic order. By default, the value is set to 1.
    -sort_tables      A flag to sort (1) or not (0) the tables by alphabetic order. By default, the value is set to 1.
                     
    Other optional options - if you want to add some SQL query results as examples
    
    -host             Host name of the MySQL server
    -port             Port of the MySQL server
    -dbname           Database name
    -user             MySQL user name
    -pass             MySQL password (not always required)
    -skip_connection  Avoid to run the MySQL queries contained in the "@example" tags.
    
  } . "\n";
  exit(0);
}
