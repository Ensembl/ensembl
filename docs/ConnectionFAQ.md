Ensembl API Connection FAQ
==========================

### MSG: human is not a valid species name (check DB and API version)

Usually this is caused by using the wrong Ensembl API version to access a server. The API can only find species from the same version as it. 

The master branch of the Ensembl git repository is typically one release ahead of the public servers, and will always fail to find a species by default. The master branch is in development, and is not guaranteed to work. To access an older release, the Registry option DB_VERSION can be set, but it is preferable to use the correct API version to avoid unintended consequences.

You can also try the systematic name of the species.

### MSG: (insert fungi/protist/etc here) is not a valid species name (check DB and API version)

Try running perl ensembl/misc-scripts/ping_ensembl.pl and check the output.

Ensembl Genomes provides these species, and releases roughly two weeks later than Ensembl. If you have just updated your API and Ensembl recently announced a release, your software may be too new for the Ensembl Genomes servers. You can wait until they release, or roll back your API version. This is easy if you installed from Github. 

```bash
VERSION=`perl -e 'use Bio::EnsEMBL::ApiVersion qw/software_version/; print software_version'`
git checkout release/`expr ${VERSION} - 1`
```

If you installed a downloaded package, then you will need to download an older Ensembl API release.

### error 2006 (MySQL server has gone away)

In a long-running process, it is possible to hit database server time limits for connections. Typically after 8 hours the server will close the connection, and your Perl code will die.

#### Strategies for avoiding timeouts

1 - Use the nearest database server to improve efficiency. Ensembl has mirrors in Asia and the USA as well as the main servers hosted in the UK. See the [Mirrors page](http://www.ensembl.org/info/about/mirrors.html) for specifics.

2 - For intense database access and high frequency querying, choose a good time to disconnect manually

```perl
...
# discrete work completed that takes an hour or two
$gene_adaptor->dbc->disconnect_if_idle;

# API re-opens connection automatically
$gene_adaptor->fetch_by_stable_id($stable_id);
```

3 - For scripts which occasionally consult Ensembl while working on a big problem for several hours

```perl
...
# For all code using Ensembl
Bio::Ensembl::Registry->set_disconnect_when_inactive(1);
# For just one adaptor
$adaptor->dbc->disconnect_when_inactive(1);
# For just one occasion
# This causes the connection to close whenever it is not being used. This is costly if there are very frequent database requests
$adaptor->dbc->disconnect_if_idle;


# In combination with disconnecting after every request, you can hold the connection open for the duration of a code block
my @gene_ids = ('ENSG0000001',...);
my %external_refs;
$gene_adaptor->dbc->prevent_disconnect(sub {
  while (my $id = shift @gene_ids ) {
    my $gene = $gene_adaptor->fetch_by_stable_id($id);
    my $xrefs = $gene->get_all_DBEntries;
    foreach my $xref (@$xrefs) {
      $external_refs{$id} = $xref->display_id;
    }
  }
});
# This will finish faster than if it continues to disconnect and reconnect all the time
```

4 - For scripts which access Ensembl a lot and have no easy opportunity to behave as in option 2 above.

```perl
# For all code using Ensembl
Bio::Ensembl::Registry->set_reconnect_when_lost(1);
# For one adaptor
$adaptor->dbc->reconnect_when_lost(1);
# This option adds an additional message to every call to the database, checking that the connection is still up
# It increases network traffic and latency of each request, but can restore a broken connection
# Not necessary if you call disconnect_if_idle
```

Option 2 is both fastest and makes best use of Ensembl servers. Option 4 is next quickest, and option 3 is slow for heavy access, but suitable for occasional requests.

### Error 99

The server can't host any more connections. Users are connecting too many times in a short time period. Contact Ensembl Helpdesk (helpdesk@ensembl.org) to let us know there is a problem, and try again later.

### DBI/DBD::mysql not in PERL5LIB

The Ensembl API requires both DBI and DBD::mysql packages, typically via cpan or cpanm. If you have installed these libraries but still have this problem, you will need to add them to your PERL5LIB environment variable.

```bash
echo $PERL5LIB
# Can I see where my libs are installed?
export PERL5LIB=$PERL5LIB:/path/to/perl/lib
```

### MSG: Cannot connect to the Ensembl MySQL server at ensembldb.ensembl.org:3306;

Firstly, check your connection parameters. Run perl ensembl/scripts/ping_ensembl.pl and see what it says. If both ping_ensembl and your script cannot connect, the most likely cause is that your local network prohibits this kind of traffic. Ask you sysadmins if they allow outbound database traffic on port 3306/5306. 
