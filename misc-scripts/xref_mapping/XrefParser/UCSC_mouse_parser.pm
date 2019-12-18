=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
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

# This module replicates the generic UCSC parser for mouse-specific data
# This prevents cross-mapping between species by treating each species as a separate source

package XrefParser::UCSC_mouse_parser;

use strict;
use warnings;

use parent qw( XrefParser::UCSCParser );


1;
