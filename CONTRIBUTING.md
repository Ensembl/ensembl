# Contribution Guide

The Ensembl development team welcomes outside contributions, in fact we moved to Git to facilitate this. However, to ensure legibility for other users, we ask contributors to take a few moments to clean up their code, its comments, and its history before submitting a pull request. It takes a little bit of effort from everyone, but nobody likes to decipher cryptic comments, review commits overloaded with minor typesetting changes or re-trace the history of a file across fragmented commits. Together, let's keep Ensembl tidy! 

This guide covers how to contribute changes to an Ensembl project. Please do not create a pull request without reading this guide first.
We also invite you to read our code of conduct (http://www.ensembl.org/info/about/code_of_conduct.html) before continuing.

## Quick Guide - Using Forks & Pull Requests

1. Fork the ensembl repository
2. Switch to the branch you want to edit
  * Consider using feature branches over editing master
3. Hack away and commit with a useful commit message
  * First line is short but descriptive
  * Other lines are more explanatory
4. Make sure your forked master is up-to date with origin and rebase if needed
6. Push
7. Create a pull request
8. Communicate what the change is
9. Wait and celebrate as your change goes in

## Quick Guide - Using Patches

**Whilst patch submission is supported we really prefer users to submit pull requests. Patches also take longer to integrate**

1. Clone the repository
2. Switch to the branch you want to edit
3. Hack away and commit
4. Use `git format-patch` to create the patches
5. Send them to helpdesk for the attention of the team

# Why Could My Pull Request Be Rejected?

We attempt to integrate as many pull requests into Ensembl as possible but do reserve some rights to reject pull requests

* The pull request brought in a non-useful change
    - Single line comments which serve no help to the code
* The pull request removes essential code
* The pull request's history was too complex
    - We want to see as close to a linear set of commits resulting in your final change
    - We will ask for multiple internal merges to be squished into a sensible set of commits
* Poor commit messages
    - Do not repeat the same message in multiple commits
* Excessive whitespace changes
    - Do not reformat someone else's code; it's just bad manners
* The pull request modifies code without testing for regression
    - Please provide simple unit tests to test the changes
* The pull request fails unit testing
    - Please ensure the test suite runs successfully and update the test data if necessary


# License

Ensembl code is licensed under our Apache 2.0 license. Our expectation is that contributing code is made available under the same license. Any copyright assertion to other organistions should be declared in the modyfying file and in the root LICENSE section.

# Using Forks and Pull Requests

This is our preferred method of accepting code contributions. It is the cleanest way and allows us to review your changes quickly before integration.

## Fork and Clone The Repository

You must fork the repository before contributing changes. The big _Fork_ button at the top right of GitHub's interface will automate the process. There is more information from GitHub at https://help.github.com/articles/fork-a-repo. 

Once forked clone your new repository to your machine to start hacking. After cloning remember to add the upstream repository as a remote. For example to add the original Ensembl core project repo back do the following:

```
git remote add upstream https://github.com/Ensembl/ensembl.git
```

## Switching Branch

By default Ensembl projects have a default branch of the **latest stable release**. If you are contributing a fix for a specific release then please remain there; otherwise switch to master.

```
git checkout --track -b master origin/master
```

To help improve your hacking time consider developing on a branch from master. This will allow you to bring in changes from upstream and integrate them into your fork without fear of merge conflicts. The following prefixes are available to use:

* _feature/_ - A new feature going into a repository
* _hotfix/_ - Fixes to be integrated into a branch
* _experimental/_ - Experimental feature with no guarentee about hash stability

Switch to a new branch once you are on master:

```
git checkout -b hotfix/quickfixtoadaptor
```

## Hacking and Committing

Go forth and hack. Be aware that commit messages should attempt to summarise the changes in a single line and then describe the changes in more depth. For example the following commit message is bad:

```
Fixed insertion code
```

Compared to:

```
Fixed insertion code when storing a Gene.

GeneAdaptor's store method was not using the DBI insert id 
retrieval system but rather the MySQL last insert id variable.
```

Try also to minimise branches within your code base. If we see too many we will ask you to rebase/squish.

## Syncing master with upstream, rebasing your changes and pushing

First switch to master and pull in new changes from upstream held on master. This will bring those changes down and attempt to merge your local master with _upstream/master_. If you have changes on master be aware that this will probably require a merge commit. Staying away from master is a good idea.

```
git checkout master
git pull upstream master
```
Once the changes are down rebase your branch against master:

```
git checkout hotfix/quickfixtoadaptor
git rebase master
```
Now push to origin:

```
git push -u origin hotfix/quickfixtoadaptor
```

## Creating the Pull Request

https://help.github.com/articles/using-pull-requests

Go to your GitHub fork's page, switch to your branch and click on the _Compare and Review_ button. This will start the merge. Then click on the top left +- file icon and edit accordingly:

* Switch the base branch to _master_

This ensures you are submitting your change against the right branch in Ensembl. For more information see [GitHub's documentation on doing this](https://help.github.com/articles/using-pull-requests#changing-the-branch-range-and-destination-repository).

Now click "Start a discussion" and give a brief description of what your pull reuqest is. Also if you are not using an obvious username let us know who you are.

## Updating a Pull Request

Should you be told to perform some fixes for your pull request you should perform them on your local repo. To update the pull request run `git push -f origin hotfix/quickfixtoadaptor` but never do this on a branch you have shared with more than just us.


# Updating the schema

Any change that affects the underlying SQL schema needs to come with the following:
- a patch file, patch_oldrelease_newrelease_version.sql
See https://github.com/Ensembl/ensembl/blob/release/81/sql/patch_79_80_b.sql for an example
This patch describes the actual schema change and updates the meta table to record the change
- an update to the table.sql file
The change from the patch should also be included in the table.sql file
This is to ensure that a newly created database is identical to an existing database to which the patch was applied
- an update to the test databases
The databases used for testing are in https://github.com/Ensembl/ensembl/tree/release/81/modules/t/test-genome-DBs
These should be patched with the latest changes to ensure data consistency.
This can be done using the patching script provided in https://github.com/Ensembl/ensembl-test/blob/master/scripts/patch_test_databases.pl

# The test suite

The Ensembl code comes with a series of unit tests that check methods behave as expected.
This test suite is run as part of the TravisCi integration for each pull request.
See https://travis-ci.org/Ensembl/ensembl
This ensures that any change does not affect existing functionality.
A pull request can only be integrated if the test suite passes successfully.
If no tests are available for a new functionality, please provide some basic tests.
