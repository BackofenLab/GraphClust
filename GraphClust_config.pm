package GraphClust_config;

use strict;
use warnings;

require Exporter;

################################################################################
## default config

our %CONFIG = (

  PATH_LOCARNA       => "/usr/local/user/locarna-1.7.2/bin/",
  PATH_RNASHAPES     => "/usr/local/rnashapes/2.1.6/bin//",
  PATH_VRNA          => "/usr/local/user/ViennaRNA-2.0.7/bin/",
  PATH_INFERNAL_1_0  => "/usr/local/infernal/1.0.2/bin//",
  PATH_INFERNAL_1_1  => "/usr/local/infernal/1.1.1/bin//",
  PATH_R             => "/usr/bin/",
  PATH_RNAZ          => "/usr/local/rnaz/2.1-a0130d9/bin//",
  PATH_BLASTCLUST    => "/usr/local/ncbiblast/2.2.15/bin//",
  PATH_OCTAVE        => "/usr/bin/",
  PATH_CMFINDER      => "/usr/local/user/cmfinder-0.2/bin/",

  VERSION_INFO       =>  "GraphClust 0.7.6",


  ## paths not automatically configured
  ## please adapt to your own system
  PATH_TMP           => "/var/tmp/",
  ## path to qsub for SGE/OGE job submission (sun/oracle grid engine)
  PATH_SGE           => "/opt/sge-6.0/bin/lx24-amd64/",
  ## hostname which is allowed to submit SGE/OGE jobs, ssh is used to login
  SGE_HOSTNAME       => "biui",
);
