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

package SeqStoreConverter::vega::MusMusculus;

use strict;
use warnings;

use SeqStoreConverter::MusMusculus;
use SeqStoreConverter::vega::VBasicConverter;

use vars qw(@ISA);

@ISA = qw(SeqStoreConverter::MusMusculus SeqStoreConverter::vega::VBasicConverter);

1;
