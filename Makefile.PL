use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'Bio::Tools::CodonOptTable',
    AUTHOR              => 'Rakesh Kumar Shardiwal <rakesh.shardiwal@gmail.com>',
    VERSION_FROM        => 'lib/Bio/Tools/CodonOptTable.pm',
    ABSTRACT_FROM       => 'lib/Bio/Tools/CodonOptTable.pm',
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Bio-Tools-CodonOptTable-*' },
);
