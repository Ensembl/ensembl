#!/usr/bin/env perl


use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use FindBin;
use Getopt::Long;
use Config::Tiny;

use JSON;
use HTTP::Request;
use LWP::UserAgent;
use Term::ReadKey;
use iCal::Parser;

use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::ApiVersion;

# ----------------
# global variables
# ----------------
my ( $dryrun, $xref_species );

main();

sub main {

  # -----------------
  # initialize logger
  # -----------------
  my $logger = Bio::EnsEMBL::Utils::Logger->new();

  # $logger->init_log();

  # ----------------------------
  # read command line parameters
  # ----------------------------
  my ( $relco, $release, $password, $help, $tickets_tsv, $ical_url, $config );

  GetOptions(
	     'relco=s'    => \$relco,
	     'release=i'  => \$release,
	     'password=s' => \$password,
	     'p=s'        => \$password,
	     'tickets=s'  => \$tickets_tsv,
	     'config=s'   => \$config,
	     'c=s'        => \$config,
             'ical=s'     => \$ical_url,
             'species=s'  => \$xref_species,
             'dryrun'     => \$dryrun,
	     'help'       => \$help,
	     'h'          => \$help,
	    );

  # ------------
  # display help
  # ------------
  if ($help) {
    usage();
  }

  # ---------------------------------
  # deal with command line parameters
  # ---------------------------------
  ( $relco, $release, $password, $tickets_tsv, $config )
    = set_parameters( $relco, $release, $password, $tickets_tsv, $config, $logger );

  # ---------------------------
  # read config file parameters
  # ---------------------------
  my $parameters = Config::Tiny->read($config);

  # check_dates($parameters);

  # integrate command line parameters to parameters object
  $parameters->{relco}       = $relco;
  $parameters->{password}    = $password;
  $parameters->{release}     = $release;
  $parameters->{tickets_tsv} = $tickets_tsv;
  $parameters->{config}      = $config;

  # If we have an ical url, try and fetch dates from there
  if($ical_url) {
      fetch_dates($ical_url, $parameters);
  }

  validate_parameters( $parameters, $logger );

  # ------------------
  # parse tickets file
  # ------------------
  my $tickets = parse_tickets_file( $parameters, $tickets_tsv, $logger );
  
  # --------------------------------
  # get existing tickets for current
  # release from the JIRA server
  # --------------------------------
  my $existing_tickets_response
    = post_request( 'rest/api/latest/search',
		    {
		     "jql" => "fixVersion = " . $parameters->{release} },
		    $parameters, $logger );
  my $existing_tickets
    = decode_json( $existing_tickets_response->content() );
  #use Data::Dumper; print Dumper $existing_tickets; exit;
  
  # --------------------
  # check for duplicates
  # --------------------
  my %tickets_to_skip;
  for my $ticket ( @{$tickets} ) {
    my $duplicate = check_for_duplicate( $ticket, $existing_tickets );

    if ($duplicate) {
      $tickets_to_skip{ $ticket->{summary} } = $duplicate;
    }
  }

  # --------------------
  # validate JIRA fields
  # --------------------
  # for my $ticket ( @{$tickets} ) {

  # $logger->info( 'Validating' . ' "' . $ticket->{summary} . '" ... ' );

  # validate_fields( $ticket, $parameters, $logger );
  # $logger->info("Done\n");

  # }
  
  # -----------------------
  # create new JIRA tickets
  # -----------------------
  for my $ticket ( @{$tickets} ) {
    $logger->info( 'Creating' . ' "' . $ticket->{summary} . '" ... ' );

    # if the ticket to be submitted is a subtask then fetch the parent key and
    # replace the parent summary with the parent key
    if ( $ticket->{'issuetype'}->{'name'} eq 'Sub-task' ) {
      my $parent_key
  	= get_parent_key( $ticket->{'parent'}, $parameters, $logger );
      $ticket->{'parent'} = { 'key'  => $parent_key };
    }

    if ( $tickets_to_skip{ $ticket->{summary} } ) {
      $logger->info(
  		    'Skipped: This seems to be a duplicate of https://www.ebi.ac.uk/panda/jira/browse/'
                    . $tickets_to_skip{ $ticket->{summary} }
                    . "\n" );
    } else {
      my $ticket_key = create_ticket( $ticket, $parameters, $logger );
      $logger->info( "Done\t" . $ticket_key . "\n" );
    }

  }

}

=head2 set_parameters

  Arg[1]      : String $relco - a Regulation team member name or JIRA username
  Arg[2]      : Integer $release - the EnsEMBL release version
  Arg[3]      : String $password - user's JIRA password
  Arg[4]      : String $tickets_tsv - path to the tsv file that holds the input
  Arg[5]      : String $config - path to the config file holding handover dates
  Arg[6]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Description : Makes sure that the parameters provided through the command line
                are valid and assigns default values to the ones which where not 
                supplied
  Return type : Listref
  Exceptions  : none

=cut

sub set_parameters {
    my ( $relco, $release, $password, $tickets_tsv, $config, $logger ) = @_;

    $relco = $ENV{'USER'} unless $relco; 
    validate_user_name( $relco, $logger );

    $release = Bio::EnsEMBL::ApiVersion->software_version() if !$release;

    $tickets_tsv = $FindBin::Bin . '/jira_recurrent_tickets.tsv'
        if !$tickets_tsv;

    if ( !-e $tickets_tsv ) {
        $logger->error(
            'Tickets file '
                . $tickets_tsv
                . ' not found! Please specify one using the -tickets option!',
            0, 0
        );
    }

    $config = $FindBin::Bin . '/jira.conf' if !$config;

    if ( !-e $config ) {
        $logger->error(
            'Config file '
                . $config
                . ' not found! Please specify one using the -config option!',
            0, 0
        );
    }

    printf( "\trelco: %s\n\trelease: %i\n\ttickets: %s\n\tconfig: %s\n",
        $relco, $release, $tickets_tsv, $config );
#    print "Are the above parameters correct? (y,N) : ";
#    my $response = readline();
#    chomp $response;
#    if ( $response ne 'y' ) {
#        $logger->error(
#            'Aborted by user. Please rerun with correct parameters.',
#            0, 0 );
#    }

    if ( !$password ) {
        print 'Please type your JIRA password:';

        ReadMode('noecho');    # make password invisible on terminal
        $password = ReadLine(0);
        chomp $password;
        ReadMode(0);           # restore typing visibility on terminal
        print "\n";
    }

    return ( $relco, $release, $password, $tickets_tsv, $config );
}

=head2 validate_parameters

  Arg[1]      : Hashref $parameters - the configuration parameters
  Arg[2]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : validate_parameters( $parameters )
  Description : Asks the user to look over the settings before proceeding
  Return type : none
  Exceptions  : Dies if the user rejects the settings

=cut

sub validate_parameters {
  my ( $parameters, $logger ) = @_;

  print "Configuration:\n";
  printf( "\trelco: %s\n\trelease: %i\n\ttickets: %s\n\tconfig: %s\n",
	  $parameters->{relco},
	  $parameters->{release},
	  $parameters->{tickets_tsv},
	  $parameters->{config} );

  print "Dates:\n";
  foreach my $date_label (keys $parameters->{dates}) {
      printf( "%33s: %s\n",
	      $date_label,
	      $parameters->{dates}->{$date_label} );
  }

  if($xref_species) {
      my @species_list = split /,/, $xref_species;

      print "Xref species:\n";

      foreach my $species (@species_list) {
	  my $valid = validate_species($species, $parameters);

	  if(! $valid ) {
	      die "Species $species doesn't seem to be valid in the metadata registry";
	  }

	  print "\t$species\n"
      }
  }

  if($dryrun) {
      print "DRY RUN.\n";
  }

  print "Are the above parameters correct? (y,N) : ";
  my $response = readline();
  chomp $response;
  if ( $response ne 'y' ) {
      $logger->error(
	  'Aborted by user. Please rerun with correct parameters.',
	  0, 0 );
  }

}

=head2 validate_species

  Arg[1]      : String $species - a latin species name
  Example     : my $valid = validate_species($species)
  Description : Validate a species via the metadata registry
  Return type : Boolean
  Exceptions  : none

=cut

sub validate_species {
  my ($species, $parameters) = @_;

  my $content = fetch_url($parameters->{_}->{metadata_registry_url} . $species,
			  'application/json');

  my $count = decode_json( $content )->{'count'};

  return ($count > 0);
}

=head2 validate_user_name

  Arg[1]      : String $user - a Core team member JIRA username
  Arg[2]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : validate_user_name($user, $logger)
  Description : Checks if the provided user name is valid
  Return type : none
  Exceptions  : none

=cut

sub validate_user_name {
  my ( $user, $logger ) = @_;

  my %valid_user_names =
    (
     'mr6' => 1,
     'avullo' => 1,
     'ktaylor' => 1,
     'lairdm' => 1,
     'prem' => 1,
     'killm9m1' => 1
    );

  return if exists $valid_user_names{$user};
    
  my $valid_names = join( "\n", sort keys %valid_user_names );
  $logger->error("User name $user not valid! Here is a list of valid names:\n" . $valid_names, 0, 0);  
}

=head2 parse_tickets_file

  Arg[1]      : Hashref $parameters - parameters from command line and config
  Arg[2]      : String $tickets_tsv - path to the tsv file that holds the input
  Arg[3]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $tickets = parse_tickets_file( $parameters, $tickets_tsv, $logger );
  Description : Reads the tsv input file, replaces placeholder tags and returns
                a listref of tickets to be submitted
  Return type : Listref
  Exceptions  : none

=cut

sub parse_tickets_file {
  my ( $parameters, $tickets_tsv, $logger ) = @_;
  my @tickets;

  # create main ticket for the given release, will have tickets read from the file
  # as sub-tasks except Xrefs, which will have its own main ticket since it contains
  # the individual species xrefs as sub-tasks and we cannot create sub-sub-tasks
  my $project = 'ENSCORESW';
  my $priority = 'Major';
  
  push @tickets,
    {
     'project'     => { 'key'  => $project },
     'issuetype'   => { 'name' => 'Task' },
     'summary'     => sprintf("Release %d", $parameters->{release}),
     # 'reporter'    => { 'name' => $reporter },
     'assignee'    => { 'name' => $parameters->{relco} },
     'priority'    => { 'name' => $priority },
     'fixVersions' => [ { 'name' => sprintf("%d", $parameters->{release}) } ],
     'duedate'     => $parameters->{dates}{release},
     # 'components'  => \@components,
     'description' => sprintf "Core activities for release %s", $parameters->{release},
    };
  
  open my $tsv, '<', $tickets_tsv;

  my $header = readline $tsv;
  chomp $header;

  while ( readline $tsv ) {
    my $line = $_;
    chomp $line;

    $line = replace_placeholders( $line, $parameters );

    my ( $summary, $reporter, $assignee,
	 $due_date, $component_string, $description,
	 $species_list ) = split /\t/, $line;

    # print "Summary: $summary\nReporter: $reporter\nAssignee: $assignee\nPriority: $priority\nDue date: $due_date\nComponents: $component_string\nDescription: $description\n\n";
    
    if($reporter) {
      validate_user_name( $reporter, $logger );
    } else {
      $reporter = $parameters->{relco};  
    }
    if ($assignee) {
      validate_user_name( $assignee, $logger );
    } else {
      $assignee = $parameters->{relco};
    }
    
    my @components;
    map { push @components, { 'name' => $_ } } split /,/, $component_string;

    my %ticket;
    
    # Xrefs is a task with its own sub-tasks, one for each species

    # If we have a global xref_species, ie we've been given the species on the command
    # line, override the one in the ticket config
    if ($summary =~ /xrefs/i) {
      if($xref_species) {
	$species_list = $xref_species;
      }

      $species_list or $logger->error("Empty list of xref species", 0, 0);
      my @species = split /,/, $species_list;

      # this is the main Xref issue
      push @tickets,
	{
	 'project'     => { 'key'  => $project },
	 'issuetype'   => { 'name' => 'Task' },
	 'summary'     => $summary,
	 # 'reporter'    => { 'name' => $reporter },
	 'assignee'    => { 'name' => $assignee },
	 'priority'    => { 'name' => $priority },
	 'fixVersions' => [ { 'name' => sprintf("%d", $parameters->{release}) } ],
	 'duedate'     => $due_date,
	 'components'  => \@components,
	 'description' => $description . "\nRun xrefs for the following species:\n\n" . join("\n", @species)
	};

      # now proceed with its sub-tasks, one for every species
      foreach my $species (@species) {
	push @tickets,
	  {
	   'project'     => { 'key'  => $project },
	   'issuetype'   => { 'name' => 'Sub-task' },
	   'summary'     => sprintf("Xrefs, %s %d", $species, $parameters->{release}),
	   'parent'      => $summary,
	   # 'reporter'    => { 'name' => $reporter },
	   'assignee'    => { 'name' => $assignee },
	   'priority'    => { 'name' => $priority },
	   'fixVersions' => [ { 'name' => sprintf("%d", $parameters->{release}) } ],
	   'duedate'     => $due_date,
	   'components'  => \@components
	  };
      }
    } else {
      # create an issue as a sub-task of the main release one
      %ticket = (
		 'project'     => { 'key'  => $project },
		 'issuetype'   => { 'name' => 'Sub-task' },
		 'summary'     => $summary,
		 # the parent summary is replaced by the parent key
		 # just before the ticket submission
		 'parent'      => sprintf("Release %d", $parameters->{release}),
		 # 'reporter'    => { 'name' => $reporter },
		 'assignee'    => { 'name' => $assignee },
		 'priority'    => { 'name' => $priority },
		 'fixVersions' => [ { 'name' => sprintf("%d", $parameters->{release}) } ],
		 'duedate'     => $due_date,
		 'components'  => \@components,
		 'description' => $description,
		);

      # delete empty fields from ticket
      for my $key ( keys %ticket ) {
	if ( !$ticket{$key} ) {
	  delete $ticket{$key};
	}
      }
      
      push @tickets, \%ticket;

    }
    
  }

  # use Data::Dumper; print Dumper \@tickets; <STDIN>;
  return \@tickets;
}

=head2 replace_placeholders

  Arg[1]      : String $line - One line from the tsv input file
  Arg[2]      : Hashref $parameters - parameters from command line and config
  Example     : $line = replace_placeholders( $line, $parameters );
  Description : Replaces the placeholder tags with valid values and returns a
                a new string
  Return type : String
  Exceptions  : none

=cut

sub replace_placeholders {
  my ( $line, $parameters ) = @_;

  for my $param (keys %$parameters) {
      next if(ref($parameters->{$param}) eq 'HASH');

      if($line =~ /<$param>/) {
	  $line =~ s/<$param>/$parameters->{$param}/eg;
      }
  }

  for my $param (keys $parameters->{dates}) {
      if($line =~ /<$param>/) {
	  $line =~ s/<$param>/$parameters->{dates}->{$param}/eg;
      }
  }

  return $line;

}

=head2 get_parent_key

  Arg[1]      : String $summary - Summary of the parent ticket
  Arg[2]      : Hashref $parameters - parameters from command line and config
  Arg[3]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $parent_key
                = get_parent_key( $ticket->{'parent'}, $parameters, $logger );
  Description : Gets the ticket key of the parent task
  Return type : String
  Exceptions  : none

=cut

sub get_parent_key {
  my ( $summary, $parameters, $logger ) = @_;

  # jql=summary ~ "Update declarations" AND fixVersion %3D release-88
  my $content
    = {   "jql" => 'fixVersion = ' . $parameters->{release}
	  . ' and summary ~ "'
	  . $summary
	  . '"' };

  my $response = post_request( 'rest/api/latest/search',
			       $content, $parameters, $logger );

  my $parent = decode_json( $response->content() )->{'issues'}->[0];

  return $parent->{'key'};
}

=head2 create_ticket

  Arg[1]      : Hashref $line - Holds the ticket data
  Arg[2]      : Hashref $parameters - parameters from command line and config
  Arg[3]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $ticket_key = create_ticket( $ticket, $parameters, $logger );
  Description : Submits a post request to the JIRA server that creates a new
                ticket. Returns the key of the created ticket
  Return type : String
  Exceptions  : none

=cut

sub create_ticket {
  my ( $ticket, $parameters, $logger ) = @_;
  my $endpoint = 'rest/api/latest/issue';

  # If we're in dry run mode, don't actually post the ticket to jira
  if($dryrun) {
      use Data::Dumper;
      print "DRY RUN, ticket:\n";
      print Dumper($ticket);
      print "\n";
      return "DRYRUN";
  }

  my $content = { 'fields' => $ticket };
  my $response = post_request( $endpoint, $content, $parameters, $logger );

  return decode_json( $response->content() )->{'key'};
}

=head2 fetch_url

  Arg[1]      : String $url - the url to request
  Arg[2]      : String $accept - the accept time (optional)
  Example     : my $response_body = fetch_url( $url );
  Description : Fetch a url's content and return as a string
  Return type : String
  Exceptions  : Dies upon fetch error

=cut

sub fetch_url {
  my ( $url, $accept ) = @_;

  my $ua = LWP::UserAgent->new;

  if($accept) {
      $ua->default_header('Accept' => $accept);
  }

  my $response = $ua->get($url);

  if($response->is_success) {
      return $response->decoded_content;
  }

  die $response->status_line;
}

=head2 post_request

  Arg[1]      : String $endpoint - the request's endpoint
  Arg[2]      : Hashref $content - the request's content
  Arg[3]      : Hashref $parameters - parameters used for authorization
  Arg[4]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $response = post_request( $endpoint, $content, $parameters, $logger )
  Description : Sends a POST request to the JIRA server
  Return type : HTTP::Response object
  Exceptions  : none

=cut

sub post_request {
    my ( $endpoint, $content, $parameters, $logger ) = @_;

    my $host = 'https://www.ebi.ac.uk/panda/jira/';
    my $url  = $host . $endpoint;

    my $json_content = encode_json($content);

    my $request = HTTP::Request->new( 'POST', $url );

    $request->authorization_basic( $parameters->{relco},
        $parameters->{password} );
    $request->header( 'Content-Type' => 'application/json' );
    $request->content($json_content);

    my $agent    = LWP::UserAgent->new();
    my $response = $agent->request($request);

    if ( $response->code() == 401 ) {
        $logger->error( 'Your JIRA password is not correct. Please try again',
            0, 0 );
    }

    if ( $response->code() == 403 ) {
        $logger->error(
            'Your do not have permission to submit JIRA tickets programmatically',
            0, 0
        );
    }

    if ( !$response->is_success() ) {
        my $error_message = $response->as_string();

        $logger->error( $error_message, 0, 0 );
    }

    return $response;
}

=head2 fetch_dates

  Arg[1]      : String $ical_url - Hold the url to fetch the ical release
                calendar
  Arg[2]      : Hashref $parameters - parameters from config file
  Example     : fetch_dates($ical_url, $parameters);
  Description : Attempts to fetch an ical calendar containing the dates
                of the release cycle
  Return type : void
  Exceptions  : none

=cut

sub fetch_dates {
  my ( $ical_url, $parameters ) = @_;

  print "Fetching ical...\n";
  my $ical_str = fetch_url( $ical_url );

  print "\n";
  my $parser=iCal::Parser->new();

  my $hash = $parser->parse_strings($ical_str);

  my $events = $hash->{events};

  foreach my $milestone_name (keys $parameters->{ical}) {
      my $calendar_label = $parameters->{ical}->{$milestone_name};
      $calendar_label = replace_placeholders($calendar_label, $parameters);
      print "Looking for '$calendar_label' in ical...";
      my $event = find_event_by_name($calendar_label, $events);

      # We've found the event, extract the start date and substitute it
      # in our dates configuration section
      if($event) {
	  $parameters->{dates}->{$milestone_name} = sprintf("%04s-%02s-%02s",
						       $event->{DTSTART}->year(),
						       $event->{DTSTART}->month(),
						       $event->{DTSTART}->day());
	  print "found.\n";
      } else {
	  print "not found.\n";
      }
  }

  print "\n";
}

=head2 find_event_by_name

  Arg[1]      : String $label - label to find in the calendar
  Arg[2]      : Hashref $events - the hash of events from the calendar
  Example     : my $event = find_event_by_name($label, $events);
  Description : Tries to find the event matching the label
  Return type : Hashref or undef
  Exceptions  : none

=cut

sub find_event_by_name {
  my ( $label, $events ) = @_;

  foreach my $year (keys %$events) {
      foreach my $month (keys $events->{$year}) {
	  foreach my $day (keys $events->{$year}->{$month}) {
	      foreach my $event (keys $events->{$year}->{$month}->{$day} ) {
		  if($events->{$year}->{$month}->{$day}->{$event}->{SUMMARY} eq $label) {
		      return $events->{$year}->{$month}->{$day}->{$event};
		  }
	      }
	  }
      }
  }

  return;
}

=head2 check_for_duplicate

  Arg[1]      : Hashref $ticket - holds the data for the ticket which is about
                to be submitted
  Arg[2]      : Hashref $existing_tickets - holds the data for all tickets that
                already exist on the JIRA server for the current EnsEMBL release
  Example     : my $duplicate = check_for_duplicate($ticket, $existing_tickets);
  Description : Checks whether the ticket which is about to be submitted exists
                already on the JIRA server and returns the relevant key if this
                is true
  Return type : String
  Exceptions  : none

=cut

sub check_for_duplicate {
    my ( $ticket, $existing_tickets ) = @_;
    my $duplicate;

    for my $existing_ticket ( @{ $existing_tickets->{issues} } ) {
        if ( $ticket->{summary} eq $existing_ticket->{fields}->{summary} ) {
            $duplicate = $existing_ticket->{key};
            last;
        }
    }

    return $duplicate;
}

sub usage {
    print <<EOF;
=head1 NAME

create_release_tasks.pl

=head1 SYNOPSIS

create_release_tasks.pl -relco <string> -password <string> -release <integer> -tickets <file> -config <file> 

-relco               JIRA username. Optional, will be inferred from current system user if not supplied.
-password | -p       JIRA password. Will need to be typed in standard input if not supplied.
-release             EnsEMBL Release. Optional, will be inferred from EnsEMBL API if not supplied.
-tickets             File that holds the input data for creating the JIRA tickets in tab separated format.
                     If not supplied, the script is looking for the default one 'jira_recurrent_tickets.tsv'
                     in the same directory as the executable.
-config              Configuration parameters file, currently holds the handover deadlines, may be expanded
                     in the future. If not supplied, the script is looking for the default one 'jira.conf'
                     in the same directory as the executable.
-ical                The URL to fetch the Release calendar from Google in iCal format
-species             A comma separated list of species to make xref tickets for, this overrides the
                     species in the jira_recurrent_tickets.tsv file.
-dryrun              Dry run, do not actually create the tickets
-help | -h           Prints this help text.

=head1 DESCRIPTION

Creates tickets for various core release tasks on JIRA at EBI.

=cut
}
EOF
    exit 0;
}

#TODO check date timeline
