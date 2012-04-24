# FASTA Pipeline

This is a re-implementation of an existing pipeline developed originally by
core and the webteam. The new version uses eHive, so familiarity with this
system is essential, and has been written to use as little memory as possible.

## The Registry File

This is the way we retrieve the database connections to work with. The
registry file should specify:

* The core (and any other) databases to dump from
* A production database
	* **species = multi**
	* **group = production**
	* Used to find which species require new DNA
* A web database
	* **species = multi**
	* **group = web**
	* Used to name BLAT index files

Here is an example of a file for v67 of Ensembl. Note the use of the
Registry object within a registry file and the scoping of the package. If
you omit the *-db_version* parameter and only use HEAD checkouts of Ensembl
then this will automatically select the latest version of the API. Any
change to version here must be reflected in the configuration file.

	package Reg;

	use Bio::EnsEMBL::Registry;
	use Bio::EnsEMBL::DBSQL::DBAdaptor;

	Bio::EnsEMBL::Registry->no_version_check(1);
	Bio::EnsEMBL::Registry->no_cache_warnings(1);
  
	{
	  my $version = 67;
	  Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
	    {
	      -host => "mydb-1",
	      -port => 3306,
	      -db_version => $version,
	      -user => "user",
	      -NO_CACHE => 1,
	    },
	    {    
	      -host => "mydb-2",
	      -port => 3306,
	      -db_version => $version,
	      -user => "user",
	      -NO_CACHE => 1,
	    },
	  );
  
	  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	    -HOST => 'mydb-2',
	    -PORT => 3306,
	    -USER => 'user',
	    -DBNAME => 'ensembl_website',
	    -SPECIES => 'multi',
	    -GROUP => 'web'
	  );
  
	  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	    -HOST => 'mydb-2',
	    -PORT => 3306,
	    -USER => 'user',
	    -DBNAME => 'ensembl_production',
	    -SPECIES => 'multi',
	    -GROUP => 'production'
	  );
	}

	1;

You give the registry to the **init_pipeline.pl** script via the **-registry** option

## Overriding Defaults Using a New Config File 

We recommend if you have a number of parameters which do not change
between releases to create a configuration file which inherits from the
root config file e.g.

	package MyCnf;

	use base qw/Bio::EnsEMBL::Pipeline::FASTA::FASTA_conf/;

	sub default_options {
	  my ($self) = @_;
	  return {
	    %{ $self->SUPER::default_options() },
    
	    #Override of options
	  };
	}

	1;

If you do override the config then you should use the package name for your overridden config in the upcoming example commands.

## Environment

### PERL5LIB

* ensembl
* ensembl-hive
* bioperl

### PATH

* ensembl-hive/scripts
* faToTwoBit (if not using a custom location)
* xdformat (if not using a custom location)

### ENSEMBL\_CVS\_ROOT\_DIR

Set to the base checkout of Ensembl. We should be able to add *ensembl-hive/sql* onto this path to find the SQL directory for hive e.g.

	export ENSEMBL_CVS_ROOT_DIR=$HOME/work/ensembl-checkouts

### ENSADMIN\_PSW

Give the password to use to log into a database server e.g.

	export ENSADMIN_PSW=wibble

## Example Commands

### To load use normally:

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -base_path /path/to/dumps -registry reg.pm

### Run a subset of species (no forcing & supports registry aliases):

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -species anolis -species celegans -species human \
	-base_path /path/to/dumps -registry reg.pm

### Specifying species to force (supports all registry aliases):

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -force_species anolis -force_species celegans -force_species human \
	-base_path /path/to/dumps -registry reg.pm

### Running & forcing a species:

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -species celegans -force_species celegans \
	-base_path /path/to/dumps -registry reg.pm

### Dumping just gene data (no DNA or ncRNA):

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -dump_type cdna \
	-base_path /path/to/dumps -registry reg.pm

### Using a different SCP user & identity:

	init_pipeline.pl Bio::EnsEMBL::Pipeline::PipeConfig:FASTA_conf \
	-pipeline_db -host=my-db-host -scp_user anotherusr -scp_identity /users/anotherusr/.pri/identity \
	-base_path /path/to/dumps -registry reg.pm

## Running the Pipeline

1. Start a screen session or get ready to run the beekeeper with a **nohup**
2. Choose a dump location
	* A fasta, blast and blat directory will be created 1 level below
3. Use an *init_pipeline.pl* configuration from above
	* Make sure to give it the **-base_path** parameter
4. Sync the database using one of the displayed from *init_pipeline.pl*
5. Run the pipeline in a loop with a good sleep between submissions and redirect log output (the following assumes you are using **bash**)
	* **2>&1** is important as this clobbers STDERR into STDOUT
	* **> my_run.log** then sends the output to this file. Use **tail -f** to track the pipeline
		beekeeper.pl -url mysql://usr:pass@server:port/db -reg_conf reg.pm -loop -sleep 5 2>&1 > my_run.log &

6. Wait

