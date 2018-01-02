=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 LICENSE

  Doxygen Pre-Processor for Perl
  Copyright (C) 2002  Bart Schuller
  Copyright (C) 2006  Phinex Informatik AG
  All Rights Reserved
  
  Doxygen Filter is free software; you can redistribute it and/or modify
  it under the same terms as Perl itself.
  
  Larry Wall's 'Artistic License' for perl can be found in
  http://www.perl.com/pub/a/language/misc/Artistic.html
  ------------------------------------------------
  Author: Aeby Thomas, Phinex Informatik AG,
 	  Based on DoxygenFilter from Bart Schuller
  E-Mail: tom.aeby@phinex.ch
 
  Phinex Informatik AG
  Thomas Aeby
  Kirchweg 52
  1735 Giffers
 
  ------------------------------------------------
  This completely rewritten version of Doxygen::PerlFilter 

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

  EnsEMBL::PerlFilter

=head1 DESCRIPTION

  Implementation of EnsEMBL::PerlFilter.
  Derived from http://www.bigsister.ch/doxygenfilter by Bart Schuller and Thomas Aeby
  Original distributed under Perl artistic license, see: http://www.bigsister.ch/doxygenfilter/license.html

  Adaptations by Kieron Taylor (ktaylor@ebi.ac.uk), 2011

  State Machine rewrite of existing filter.
  Was going to use DFA::Command to handle the logic, but actually it won't work well for parsing Perl, 
  hence own simplified state machine. Intentionally not using PPI package to parse Perl, too complex.
  This is a 80/20 EnsEMBL specific POD->Doxygen converter, however it should work somewhat with other code.
 

=cut 

package EnsEMBL::PerlFilter;

use warnings;
use strict;

use base qw(EnsEMBL::Filter);

# Possible states
use constant {
    NORMAL      => 0,
    INHERIT     => 1,
    PODTOP      => 2,
    PODSECTION  => 3,
    PODMETHOD   => 4,
    SEEALSO	    => 5,
    TERMINAL    => 6, #would have been END, but that's a reserved word
    CODE        => 7,
};

# State-determining functions
my @parse = ( \&normal_parser, \&inheritance_parser, \&pod_top_parser, \&pod_section_parser, \&pod_method_parser, \&see_also_parser, \&finish, \&code_parser );
# State-reactive functions
my @act = ( \&normal_action, \&inheritance_action, \&pod_top_action, \&pod_section_action, \&pod_method_action, \&see_also_action, \&finish, \&code_action );

my @buffer;
my $state; # state of state machine, see?

my @big_buffer; # to absorb everything we want to print right up until we know whom we inherit from.
my $class_declaration;
my @inheritance;
my @leading_text;
my $method_description;
my $previous_doc_header;
my $brackets =0;

my $id = __PACKAGE__;

sub filter {
    my($self, $infile) = @_;
    open(my $infh, $infile);
    my $current_class = "";
    $state = NORMAL;
    my $line;
    # Read file, using lookup table to run correct parser on each line.
    # Parsing is done in two stages:
    #   1) Parsing based on the previous mode. The methods are called "_parse"
    #   2) Parsing after any mode switch in stage 1. The methods are called "_action"
    # Strictly speaking actions are also performed in "parse" code, but that's what was needed in the end.
    while( defined($line = <$infh>) && $state != TERMINAL) {
		my $sub_ref = $parse[$state];
		$parse[$state]->([$self,$line]);
		$act[$state]->([$self,$line]);
    }
    # Create the filtered file:
    # beware, #include declarations are coming from elsewhere ( inheritance_action() ).
    
    my @namespaces;
    my $class_name;
    if (defined($class_declaration)) {
		@namespaces = split(/::/,$class_declaration);
		$class_name = pop @namespaces;
		foreach (@namespaces) {
	    	$self->print("namespace ".$_." {\n");
		}
		$self->more(@leading_text);
		$self->print("class ".$class_name);
    }
    else {
		$self->print("# No class definition in this file.");
		warn "No package line found in $infile\n";
    }
    if (scalar @inheritance > 0) {
		my $string = shift @inheritance;
		$self->print(" : public ".$string);
		foreach my $parent (@inheritance) {
	    	$self->print(", public ".$parent);
		}
    }
    $self->print(" {\n");
    $self->print("public: \n");
    $self->more(@big_buffer);
    $self->print("};\n");
    foreach (@namespaces) {
		print("}\n");
    }
}


my $packaged;
sub normal_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    
    chomp($line);
	
	if ($line =~ /^=head1/) { 
	    $state = PODTOP;
	}
	elsif ($line =~/^=head2/) { # head2 usually signifies a doc-block just before a method
	    $state = PODMETHOD;
	}
	elsif ($line =~/^1;/) {
	    $state = TERMINAL;
	    warn "Reached end of code: 1;\n";
	}
	elsif ($line =~ /^\s*package\s+(.*);/) { 
        if ($packaged) {$state = TERMINAL} else {
            $class_declaration = $1;
            $packaged = 1;
        }
	}
	elsif ($line =~/^\s*use\s/ || $line =~/^(our|my)?\s*\@ISA/ || $line =~(/^\s*.*::ISA/) ){
	    $state = INHERIT;
	}
    
}

sub normal_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    # for catching undocumented subroutines and adding code blocks to documentation
    if ($line =~ /^\s*sub\s+([\w:]+)/) {
		$state = CODE;
		my $method_name = $1;
		if ($line =~ /{/) {$brackets = 1;} else {$brackets = 0;}
		#warn "Previous: $previous_doc_header. Present: $method_name\n";
		if (defined($previous_doc_header) && $previous_doc_header =~ /$method_name/){
		    # We've found the corresponding sub to go with the documentation.
		    $previous_doc_header = "";
		}
		else {
		    # Create an undocumented entry
		    my $scope = "public";
		    if ($method_name =~ /^_/) {$scope = "protected";}
		    push @big_buffer,"/** \@fn $scope $method_name( ) \n Undocumented method\n\n";
	        $method_description = $scope." ".$method_name;
		    warn "Found undocumented method $method_name\n";
		}
		my $html_lump = "
<div id='codesection-$method_name' class='dynheader closed' style='cursor:pointer;' onclick='return toggleVisibility(this)'>
\@htmlonly 
	<img id='codesection-$method_name-trigger' src='closed.png' style='display:inline'><b>Code:</b>
</div>
<div id='codesection-$method_name-summary' class='dyncontent' style='display:block;font-size:small;'>click to view</div>
<div id='codesection-$method_name-content' class='dyncontent' style='display: none;'> 
\@endhtmlonly
\@code\n";
        #push @big_buffer,"\@par Code:\n\@code\n";
        push @big_buffer,$html_lump;
        push @big_buffer,$line;
        if ($line =~ /sub.*{.*}/ || $line =~ /^\s*1;/) {$brackets = 0;} #one-line subroutines must be catered for, but only after the magic <Div> is created.
    }
}

sub inheritance_parser {
    # We only get here if dealing with multiple inheritance over comma-terminated multiple lines
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    
    chomp($line);
    $line =~ s/;//;
    my @parents = $line =~ /Bio::EnsEMBL::[\w:]+/g;
    push @inheritance,@parents;
}
sub inheritance_action {
    my ($include, $inherit);
    
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    chomp($line);
    $line =~ s/;//;
    # simple inheritance suited only to Ensembl code. Multiple inheritance from one line possible
    # There are a few ignored cases of Bio::PrimarySeqI and other things from BioPerl(?)
    if ($line =~ /\@ISA/) {
	   my @parents = $line =~ /Bio::EnsEMBL::[\w:]+/g;
	   push @inheritance,@parents;
    }
    elsif ($line =~ /use base/ or $line =~ /use parent/) {
        my @parents = $line =~ /Bio::EnsEMBL::[\w:]+/g;
        push @inheritance,@parents;
    }
    elsif ($line =~ /use/) {
        $line =~ /use\s+([\w:]+)/;
        $include = $1;
        if (defined($include)) {
            unless ($include eq "strict" || $include eq "warnings" || $include eq "vars" || $include eq "Exporter" || $include eq "base") {
                $include =~ s/::/\//g; # allows doxygen to know where to look for other packages
                $self->print("#include \"".$include.".pm\"\n");
            }
        }
        else {
            warn "Inheritance issue with: $line";   
        }
    }
    
    $line =~ s/\s*#.*$//; # strip trailing comments
    if ($line =~ /,$/) { 
        $state = INHERIT; # technically redundant, but easier to follow state 
    } else {
        $state = NORMAL;
    }
}

sub pod_top_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    #secondary switching logic found elsewhere in pod_section_parser, due to unterminated pod sections. 
    if ($line =~ /DESCRIPTION|SYNOPSIS/) {$state = PODSECTION;} 
    elsif ($line =~ /SEE ALSO|OPTIONS/) {$state = SEEALSO;} 
    elsif ($line =~ /^=cut/) {$state = NORMAL;}
}
sub pod_top_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];

}

sub pod_section_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    if ($line =~ /^=head1\s+(.+)|^(=cut)/) {
        #end of section. Flush out, otherwise keep on slurping through pod_section_action
        my $header = $1;	
        if ($buffer[0] =~ /DESCRIPTION/) {
            push @leading_text,"/**  \@section Description\n\@htmlonly<pre>";
            shift @buffer; #discard the description pod header
		    foreach (@buffer) {
                $_ =~ s/\@/\\@/g; # escape @array references but only in descriptions.
                $_ =~ s/=head(\d)\s*(.*)/<\/pre>\n<h$1>$2<\/h$1>\n<pre>/;   # replace in-block head commands with formatting
            }
            push @leading_text,@buffer; 
            push @leading_text,"</pre>\@endhtmlonly*/ \n";
            @buffer = ();
        }
        elsif ($buffer[0] =~ /SYNOPSIS/) { 
            push @leading_text,"/**  \@section Synopsis\n\@code\n";
            shift @buffer;
            push @leading_text,@buffer; 
            push @leading_text,"\@endcode */ \n";
            @buffer = ();
        }
        if (defined($header) && ( $header eq "DESCRIPTION" || $header eq "SYNOPSIS") ) {
            $state = PODSECTION;
        }
        elsif (defined($header) && $header eq "SEE ALSO") {
            $state = SEEALSO;
        }
        elsif (not defined($header)) {
            $state = NORMAL; # this fires when the =cut pattern matches.
        }
        else {
            $state = PODTOP;
        }
    }
}
sub pod_section_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    $line =~ s/[BICLFS]<(.+?)>/$1/g; # remove POD formatting commands
    #$line =~ s/(<|>)/\\$1/g; #protect HTML-like stuff that isn't HTML
    push @buffer,$line;
}

sub pod_method_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    # optional extraction of parameters by guesswork?
    my $proto = "";
    if ($line =~ /^=cut/) {
		#flush out method
		$method_description = shift @buffer;
		chomp $method_description;
		$method_description =~ s/retval//; # remove any still unassigned return types
		# trim trailing brackety stuff off method header. It is upsetting Doxygen
		$method_description =~ s/\s*\(.+\).*$//;
		# Protect the code examples with a @code block
	    my $code;
		foreach (@buffer) {
		    if ($code and /\w+\s*:(?!:)/) { #closes code section at end of Example: 
                $_ = '@endcode <pre>'.$_;
                $code = undef;
            }
		    if (/Example\s+:(.*)/) { #starts code section at Example : in POD
		        $_ =~ s#(Examples?\s+:)(.*)#$1</pre>\@code$2#;
		        $code = 1;
		    }
		}
		if ($code) {push @buffer,'\@endcode'};
		push @big_buffer,"\n /** \@fn ".$method_description."\n<pre>";
		push @big_buffer,@buffer;
		push @big_buffer," </pre>\n "; # this comment block is still open. To be finished in code_parser
		@buffer = ();
		$previous_doc_header = $method_description;
		$state = NORMAL;
    }
    my $return_type;
    if ($line =~ /^\s*returns?\s*(type)?\s*:/i) { #picking up "Returns : ", "Return type:"
		$return_type = $self->sanitize_return_values($line);
		if (not defined($return_type)) {$return_type = "";}
		$buffer[0] =~ s/retval/$return_type/;
    }
    
}

sub pod_method_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    my ($protection,$method_name);
    # extract method name from header, including methods with spaces in names, unnecessary brackets on the ends and so on
    if ($line =~ /^=head2\s+([\w_\-\&\s]+)/) {
		$method_name = $1;
        chomp $method_name;
        # use _method coding convention to identify scope of method
        if( substr( $method_name, 0, 1 ) eq "_" ) {
            $protection = "protected";
        }
        else {
            $protection = "public";
        }
        $method_description = "$protection retval $method_name";
        push @buffer,$method_description."( );\n";
    }
    else {
        $line =~ s/[BICLFS]<(.+?)>/$1/g; # remove POD formatting commands
        #$line =~ s/(\@|&|<|>|\\|\%|#)/\\$1/g; #sanitising the oddities that will bewilder Doxygen 
		$line =~ s/(?<!isn't\s)DEPRECATED/\@deprecated/i; #make use of Doxygen's deprecated list features
        push @buffer,$line;
    }
}

sub see_also_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    
    if ($line =~ /^=cut/ || $line =~ /^=head1/) {
        $buffer[0] =~ s/=head1/\/**\@section/;
        $buffer[0] =~ s/SEE\sALSO/See-also/; #join up SEE ALSO so that doxygen doesn't chop off the SEE
        push @buffer,'*/';
        push @leading_text,@buffer;
        @buffer = ();
        $state = NORMAL;
    }
}
sub see_also_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    push @buffer,$line;
    
}

sub finish {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
}


sub code_parser {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    
    #count the brackets in the code to know when we've run out of subroutine
    $line =~ s/\s#.*?$//; #rips comments off the end of a line prior to counting
    $line =~ s#/.*/##g; #remove any conventional regexp from the line, as they can contain brackets
    my $open_brackets = () = $line =~ /{/g;
    my $close_brackets = () = $line =~ /}/g;
    $brackets = $brackets + $open_brackets - $close_brackets;
}

sub code_action {
    my $args = $_[0];
    my $self = $args->[0];
    my $line = $args->[1];
    
    push @big_buffer,$line; 
    # When we run out of open brackets, or we hit weird unpaired brackets in strings or comments
    if ($brackets <=0 || $line =~ /^=/ || $line =~ /^\s*sub[\s{]/ || $line =~ /^\s*1;/) {
        $state = NORMAL;
        push @big_buffer,"\@endcode\n </div>*/\n";
    	#Add fake function for doxygen to find after the comment.
		push @big_buffer,$method_description."( );\n";
    }
}

# Methods sometimes have well-specified return values, these must be dressed up for Doxygen.
sub sanitize_return_values {
    my ( $self, $line ) = @_;
    my $return_value;
    my $type;
    $line =~ /^\s*returns?\s*(type)?\s*:\s*([\w:]+)/i;
    if ($2) {
		my $type = $2;
		if ($type =~ /^int/i) {$return_value = "Int";}
		elsif ($type eq "undef") {$return_value = "Undef";}
		elsif ($type =~ /^tri/i) {$return_value = "Boolean Or Undef";}
		elsif ($type =~ /none/i) {$return_value = "void";}
		elsif ($type eq "1" || $type eq "0" || $type =~ /TRUE/i || $type =~ /FALSE/i) {$return_value = "Boolean";}
		elsif ($line =~ /\sSQL\s/) {$return_value = "SQLStatement";}
		elsif ($line =~ /subclass type/) {$return_value = "\$this";}
		else {$return_value = ucfirst($type);}
    }
    else {
		$type = "void";
    }
    return $return_value;
}

1;
